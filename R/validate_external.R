#' Validate a gene signature using an external GEO cohort
#'
#' @param signature A CancerSignature object from load_signature()
#' @param geo_id Character. GEO accession number (e.g., 'GSE68465')
#' @param surv_time_col Column name in pData for survival time
#' @param surv_event_col Column name in pData for survival event
#' @param stage_col Column name in pData for TNM stage (or NULL)
#' @param gene_col Column name in fData for gene symbols
#' @param event_value Value indicating event (default: 'deceased')
#' @param log2_transform Apply log2(x+1)? Default: TRUE
#' @param cache_dir Directory to cache GEO data
#' @return ExternalValidation object
#' @export
validate_external <- function(signature,
                               geo_id,
                               surv_time_col,
                               surv_event_col,
                               stage_col      = NULL,
                               gene_col       = 'Gene Symbol',
                               event_value    = 'deceased',
                               log2_transform = TRUE,
                               cache_dir      = tempdir()) {

  if (!inherits(signature, 'CancerSignature')) {
    stop('`signature` must be a CancerSignature object from load_signature()')
  }
  if (!requireNamespace('GEOquery', quietly = TRUE)) {
    stop('Package GEOquery is required. Install: BiocManager::install(GEOquery)')
  }

  message('Downloading GEO dataset: ', geo_id, ' ...')

  gse <- tryCatch(
    GEOquery::getGEO(geo_id, GSEMatrix = TRUE, destdir = cache_dir,
                     AnnotGPL = FALSE, getGPL = FALSE),
    error = function(e) stop('Failed to download ', geo_id, ': ', e$message)
  )
  if (is.list(gse)) gse <- gse[[1]]

  expr_raw <- Biobase::exprs(gse)
  pheno    <- Biobase::pData(gse)
  fdata    <- Biobase::fData(gse)

  message('  Matrix: ', nrow(expr_raw), ' probes x ', ncol(expr_raw), ' samples')

  # Gene symbol kolonu bul
  if (!gene_col %in% colnames(fdata)) {
    auto_col <- colnames(fdata)[grepl('symbol|gene_name|gene$',
                                       tolower(colnames(fdata)))][1]
    if (!is.na(auto_col)) {
      message('  Using column: ', auto_col)
      gene_col <- auto_col
    } else {
      stop('Gene symbol column not found. Available: ',
           paste(head(colnames(fdata), 10), collapse = ', '))
    }
  }

  gene_symbols <- fdata[[gene_col]]
  valid_probes <- !is.na(gene_symbols) & gene_symbols != '' & gene_symbols != '---'
  expr_raw     <- expr_raw[valid_probes, ]
  gene_symbols <- gene_symbols[valid_probes]
  gene_symbols <- trimws(sub('\\s*///.*', '', gene_symbols))

  # Duplicate probeler: en yuksek ortalama ekspresyonlu probeyi tut
  row_means  <- rowMeans(expr_raw, na.rm = TRUE)
  best_idx   <- tapply(seq_len(nrow(expr_raw)), gene_symbols,
                       function(idx) idx[which.max(row_means[idx])])
  best_idx   <- unlist(best_idx)
  expr_mat   <- expr_raw[best_idx, ]
  rownames(expr_mat) <- names(best_idx)

  # Log2 donusumu
  if (log2_transform && max(expr_mat, na.rm = TRUE) > 50) {
    expr_mat <- log2(expr_mat + 1)
    message('  Applied log2(x+1) transformation.')
  }

  message('  Gene matrix: ', nrow(expr_mat), ' x ', ncol(expr_mat))

  # Survival kolonu bul
  if (!surv_time_col %in% colnames(pheno)) {
    match_col <- colnames(pheno)[grepl(surv_time_col, colnames(pheno),
                                        ignore.case = TRUE)][1]
    if (!is.na(match_col)) {
      message('  Using time column: ', match_col)
      surv_time_col <- match_col
    } else {
      stop('surv_time_col not found. Available: ',
           paste(head(colnames(pheno), 20), collapse = ', '))
    }
  }

  if (!surv_event_col %in% colnames(pheno)) {
    match_col <- colnames(pheno)[grepl(surv_event_col, colnames(pheno),
                                        ignore.case = TRUE)][1]
    if (!is.na(match_col)) {
      message('  Using event column: ', match_col)
      surv_event_col <- match_col
    } else {
      stop('surv_event_col not found. Available: ',
           paste(head(colnames(pheno), 20), collapse = ', '))
    }
  }

  os_time  <- suppressWarnings(as.numeric(as.character(pheno[[surv_time_col]])))
  event_lower <- tolower(trimws(as.character(pheno[[surv_event_col]])))
  os_event <- as.integer(event_lower %in%
                c('1', 'true', 'deceased', 'dead', 'yes', tolower(event_value)))

  names(os_time)  <- colnames(expr_mat)
  names(os_event) <- colnames(expr_mat)

  valid    <- !is.na(os_time) & !is.na(os_event) & os_time > 0
  expr_mat <- expr_mat[, valid]
  os_time  <- os_time[valid]
  os_event <- os_event[valid]

  message('  Valid patients: ', sum(valid))

  # Stage verisi (opsiyonel)
  stage_vec <- NULL
  if (!is.null(stage_col) && stage_col %in% colnames(pheno)) {
    stage_vec        <- pheno[[stage_col]][valid]
    names(stage_vec) <- names(os_time)
    message('  Stage column found: ', stage_col)
  }

  # Imzayi skorla
  message('  Scoring: ', signature$name, ' ...')
  scored <- score_signature(signature = signature, expr_matrix = expr_mat)
  message('  Coverage: ', scored$coverage_pct, '%')

  if (scored$coverage_pct < 50) {
    warning('Low gene coverage (', scored$coverage_pct, '%). Missing: ',
            paste(scored$genes_missing, collapse = ', '))
  }

  # Gun -> ay donusumu
  common_ids <- intersect(names(scored$risk_scores), names(os_time))
  os_time_m  <- os_time[common_ids] / 30.44
  os_event_m <- os_event[common_ids]

  # C-indeks hesapla
  cfit      <- survival::concordance(
    survival::Surv(os_time_m, os_event_m) ~ scored$risk_scores[common_ids]
  )
  c_raw     <- cfit$concordance
  c_index   <- round(max(c_raw, 1 - c_raw), 4)
  c_se      <- round(sqrt(cfit$var), 4)
  km_fit    <- survival::survdiff(
    survival::Surv(os_time_m, os_event_m) ~ scored$risk_group[common_ids]
  )
  logrank_p <- round(1 - pchisq(km_fit$chisq, df = 1), 4)

  # Staging karsilastirmasi
  staging_result <- NULL
  if (!is.null(stage_vec)) {
    staging_result <- tryCatch(
      compare_to_staging(
        scored     = scored,
        surv_time  = os_time_m,
        surv_event = os_event_m,
        stage      = stage_vec[common_ids]
      ),
      error = function(e) { warning('Staging failed: ', e$message); NULL }
    )
  }

  # Sonuc
  result <- list(
    geo_id         = geo_id,
    signature_name = signature$name,
    cancer_type    = signature$cancer,
    n_patients     = length(common_ids),
    genes_used     = scored$genes_used,
    genes_missing  = scored$genes_missing,
    coverage_pct   = scored$coverage_pct,
    c_index        = c_index,
    c_index_se     = c_se,
    logrank_p      = logrank_p,
    scored         = scored,
    staging        = staging_result,
    surv_time      = os_time_m,
    surv_event     = os_event_m,
    expr_matrix    = expr_mat
  )

  class(result) <- 'ExternalValidation'
  message('Done. C-index = ', c_index, ' (p = ', logrank_p, ')')
  return(result)
}

#' Print method for ExternalValidation objects
#' @param x An ExternalValidation object
#' @param ... Further arguments (ignored)
#' @export
print.ExternalValidation <- function(x, ...) {
  cat('== ExternalValidation =======================\n')
  cat('  GEO dataset  :', x$geo_id, '\n')
  cat('  Signature    :', x$signature_name, '\n')
  cat('  Cancer type  :', x$cancer_type, '\n')
  cat('  Patients     :', x$n_patients, '\n')
  cat('  Gene coverage:', x$coverage_pct, '%',
      '(', length(x$genes_used), 'genes )\n')
  if (length(x$genes_missing) > 0) {
    cat('  Missing genes:', paste(x$genes_missing, collapse = ', '), '\n')
  }
  cat('\n  Performance:\n')
  cat(sprintf('    C-index   : %.4f (SE=%.4f)\n', x$c_index, x$c_index_se))
  cat(sprintf('    Log-rank p: %.4f\n', x$logrank_p))
  if (!is.null(x$staging)) {
    cat('\n  vs. TNM staging:\n')
    cat(sprintf('    Sig C     : %.4f\n', x$staging$sig_c_index))
    cat(sprintf('    TNM C     : %.4f\n', x$staging$stage_c_index))
    cat(sprintf('    Delta C   : %+.4f (p=%.4f)\n',
                x$staging$c_index_diff, x$staging$p_value_diff))
    cat('   ', x$staging$verdict, '\n')
  }
  cat('=============================================\n')
  invisible(x)
}
