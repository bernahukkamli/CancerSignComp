#' Compare multiple cancer prognostic gene signatures on the same dataset
#'
#' @param signatures A list of \code{CancerSignature} objects
#' @param expr_matrix A numeric matrix: rows = genes, columns = patients
#' @param surv_time Numeric vector of survival times (same order as matrix columns)
#' @param surv_event Numeric vector: 1 = event occurred, 0 = censored
#' @param cutoff Cutoff method passed to \code{score_signature()}:
#'   \code{"median"} (default), \code{"mean"}, or numeric quantile
#'
#' @return A list of class \code{SignatureComparison} with per-signature
#'   metrics and a summary table
#' @export
#'
#' @examples
#' set.seed(42)
#' mat <- matrix(rnorm(1000), nrow = 10,
#'               dimnames = list(
#'                 c("TP53","KRAS","MYC","CDKN2A","SMAD4",
#'                   "EGFR","PTEN","RB1","BRCA2","PIK3CA"),
#'                 paste0("Patient_", 1:100)))
#' sig1 <- load_signature(c("TP53","KRAS","MYC"), "Sig_A", "PAAD")
#' sig2 <- load_signature(c("EGFR","PTEN","RB1"), "Sig_B", "PAAD")
#' surv_t <- rexp(100, rate = 0.05)
#' surv_e <- rbinom(100, 1, 0.7)
#' result <- compare_signatures(list(sig1, sig2), mat, surv_t, surv_e)
#' print(result)
compare_signatures <- function(signatures,
                               expr_matrix,
                               surv_time,
                               surv_event,
                               cutoff = "median") {

  # --- Girdi kontrolleri ---
  if (!is.list(signatures) || length(signatures) < 2) {
    stop("`signatures` must be a list of at least 2 CancerSignature objects")
  }

  # ScoredSignature verilmisse CancerSignature listesine donustur
  signatures <- lapply(seq_along(signatures), function(i) {
    s <- signatures[[i]]
    if (inherits(s, "ScoredSignature")) {
      # ScoredSignature icerisinden CancerSignature olustur
      load_signature(
        genes       = s$genes_used,
        name        = s$signature_name,
        cancer_type = s$cancer_type
      )
    } else if (inherits(s, "CancerSignature")) {
      s
    } else {
      stop("Element ", i, " must be a CancerSignature or ScoredSignature object.\n",
           "  Got: ", class(s)[1], "\n",
           "  Use load_signature() to create a CancerSignature object.")
    }
  })

  n_patients <- ncol(expr_matrix)

  if (length(surv_time) != n_patients) {
    stop("`surv_time` length (", length(surv_time),
         ") must equal number of patients (", n_patients, ")")
  }

  if (length(surv_event) != n_patients) {
    stop("`surv_event` length (", length(surv_event),
         ") must equal number of patients (", n_patients, ")")
  }

  if (!all(surv_event %in% c(0, 1))) {
    stop("`surv_event` must contain only 0 (censored) and 1 (event)")
  }

  # --- survival paketi gerekli ---
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required. Install with: install.packages('survival')")
  }

  # --- Her imza için metrik hesapla ---
  metrics_list <- lapply(signatures, function(sig) {

    # Skoru hesapla
    scored <- score_signature(sig, expr_matrix, cutoff = cutoff)

    # C-index (concordance index)
    surv_obj  <- survival::Surv(surv_time, surv_event)
    cfit      <- survival::concordance(surv_obj ~ scored$risk_scores)
    c_raw     <- cfit$concordance
    c_index   <- round(max(c_raw, 1 - c_raw), 4)
    c_se      <- round(sqrt(cfit$var), 4)

    # Log-rank p-değeri (Low vs High risk grubu)
    group     <- scored$risk_group
    km_fit    <- survival::survdiff(surv_obj ~ group)
    logrank_p <- round(1 - pchisq(km_fit$chisq, df = 1), 4)

    # Medyan survival (Low vs High)
    km_surv   <- survival::survfit(surv_obj ~ group)
    med_surv  <- summary(km_surv)$table[, "median"]

    list(
      name          = sig$name,
      cancer_type   = sig$cancer_type,
      n_genes       = sig$n_genes,
      genes_used    = length(scored$genes_used),
      coverage_pct  = scored$coverage_pct,
      c_index       = c_index,
      c_index_se    = c_se,
      logrank_p     = logrank_p,
      n_low         = sum(scored$risk_group == "Low"),
      n_high        = sum(scored$risk_group == "High"),
      median_surv_low  = unname(med_surv[1]),
      median_surv_high = unname(med_surv[2]),
      scored        = scored
    )
  })

  names(metrics_list) <- sapply(signatures, function(s) s$name)

  # --- Özet tablo oluştur ---
  summary_df <- data.frame(
    signature    = sapply(metrics_list, `[[`, "name"),
    cancer_type  = sapply(metrics_list, `[[`, "cancer_type"),
    n_genes      = sapply(metrics_list, `[[`, "n_genes"),
    coverage_pct = sapply(metrics_list, `[[`, "coverage_pct"),
    c_index      = sapply(metrics_list, `[[`, "c_index"),
    c_index_se   = sapply(metrics_list, `[[`, "c_index_se"),
    logrank_p    = sapply(metrics_list, `[[`, "logrank_p"),
    n_low        = sapply(metrics_list, `[[`, "n_low"),
    n_high       = sapply(metrics_list, `[[`, "n_high"),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  # C-index'e göre sırala (büyükten küçüğe)
  summary_df <- summary_df[order(summary_df$c_index, decreasing = TRUE), ]

  # --- SignatureComparison objesi ---
  result <- list(
    summary_table = summary_df,
    metrics       = metrics_list,
    n_signatures  = length(signatures),
    n_patients    = n_patients,
    cutoff_method = cutoff
  )

  class(result) <- "SignatureComparison"
  return(result)
}


#' Print method for SignatureComparison objects
#'
#' @param x A \code{SignatureComparison} object
#' @param ... Further arguments (ignored)
#' @export
print.SignatureComparison <- function(x, ...) {
  cat("== SignatureComparison ======================\n")
  cat("  Signatures  :", x$n_signatures, "\n")
  cat("  Patients    :", x$n_patients, "\n")
  cat("  Cutoff      :", x$cutoff_method, "\n\n")
  cat("  Ranked by C-index:\n")
  df <- x$summary_table
  for (i in seq_len(nrow(df))) {
    cat(sprintf("  [%d] %-25s C-index=%.4f  p=%.4f  coverage=%.1f%%\n",
                i,
                df$signature[i],
                df$c_index[i],
                df$logrank_p[i],
                df$coverage_pct[i]))
  }
  cat("=============================================\n")
  invisible(x)
}
