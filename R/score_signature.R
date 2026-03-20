#' Score patients using a cancer prognostic gene signature
#'
#' @param signature  A \code{CancerSignature} object from \code{load_signature()}
#' @param expr_matrix A numeric matrix: rows = genes, columns = patients.
#'   Row names must be gene symbols (e.g., "TP53", "KRAS").
#' @param weights Optional named numeric vector of gene weights.
#'   Names must match gene symbols. If NULL, equal weights are used.
#' @param cutoff Method to split patients into risk groups:
#'   \code{"median"} (default), \code{"mean"}, or a numeric value (0-1 quantile).
#'
#' @return A list of class \code{ScoredSignature} containing risk scores
#'   and group assignments for each patient.
#' @export
#'
#' @examples
#' # Simulated expression matrix (100 patients, 5 genes)
#' set.seed(42)
#' mat <- matrix(rnorm(500), nrow = 5,
#'               dimnames = list(c("TP53","KRAS","MYC","CDKN2A","SMAD4"),
#'                               paste0("Patient_", 1:100)))
#' sig <- load_signature(c("TP53","KRAS","MYC"), "TestSig", "PAAD")
#' scored <- score_signature(sig, mat)
#' print(scored)
score_signature <- function(signature,
                            expr_matrix,
                            weights = NULL,
                            cutoff  = "median") {

  # --- Girdi kontrolleri ---
  if (!inherits(signature, "CancerSignature")) {
    stop("`signature` must be a CancerSignature object from load_signature()")
  }

  if (!is.matrix(expr_matrix) || !is.numeric(expr_matrix)) {
    stop("`expr_matrix` must be a numeric matrix (rows=genes, cols=patients)")
  }

  if (is.null(rownames(expr_matrix))) {
    stop("`expr_matrix` must have row names (gene symbols)")
  }

  if (is.null(colnames(expr_matrix))) {
    stop("`expr_matrix` must have column names (patient IDs)")
  }

  # --- İmza genleri matriste var mı? ---
  genes_found   <- intersect(signature$genes, rownames(expr_matrix))
  genes_missing <- setdiff(signature$genes, rownames(expr_matrix))

  if (length(genes_found) == 0) {
    stop("None of the signature genes found in expr_matrix row names.")
  }

  if (length(genes_missing) > 0) {
    warning(length(genes_missing), " signature gene(s) not found in matrix: ",
            paste(genes_missing, collapse = ", "))
  }

  coverage_pct <- round(length(genes_found) / length(signature$genes) * 100, 1)

  # --- Ağırlıkları hazırla ---
  if (is.null(weights)) {
    w <- rep(1, length(genes_found))
    names(w) <- genes_found
  } else {
    if (!is.numeric(weights) || is.null(names(weights))) {
      stop("`weights` must be a named numeric vector")
    }
    w <- weights[genes_found]
    if (any(is.na(w))) {
      missing_w <- genes_found[is.na(w)]
      warning("No weight provided for: ", paste(missing_w, collapse = ", "),
              ". Using weight = 1 for these genes.")
      w[is.na(w)] <- 1
    }
  }

  # --- Risk skoru hesapla (ağırlıklı ortalama) ---
  sub_matrix  <- expr_matrix[genes_found, , drop = FALSE]
  risk_scores <- apply(sub_matrix, 2, function(patient_expr) {
    weighted.mean(patient_expr, w = w)
  })

  # --- Risk gruplarını belirle ---
  threshold <- if (is.numeric(cutoff)) {
    quantile(risk_scores, probs = cutoff)
  } else if (cutoff == "median") {
    median(risk_scores)
  } else if (cutoff == "mean") {
    mean(risk_scores)
  } else {
    stop("`cutoff` must be 'median', 'mean', or a numeric value between 0 and 1")
  }

  risk_group <- ifelse(risk_scores >= threshold, "High", "Low")

  # --- ScoredSignature objesi oluştur ---
  result <- list(
    signature_name = signature$name,
    cancer_type    = signature$cancer_type,
    genes_used     = genes_found,
    genes_missing  = genes_missing,
    coverage_pct   = coverage_pct,
    n_patients     = ncol(expr_matrix),
    risk_scores    = risk_scores,
    risk_group     = factor(risk_group, levels = c("Low", "High")),
    threshold      = threshold,
    cutoff_method  = ifelse(is.numeric(cutoff), paste0("quantile_", cutoff), cutoff),
    weights_used   = w
  )

  class(result) <- "ScoredSignature"
  return(result)
}


#' Print method for ScoredSignature objects
#'
#' @param x A \code{ScoredSignature} object
#' @param ... Further arguments (ignored)
#' @export
print.ScoredSignature <- function(x, ...) {
  cat("== ScoredSignature ==========================\n")
  cat("  Signature   :", x$signature_name, "\n")
  cat("  Cancer type :", x$cancer_type, "\n")
  cat("  Patients    :", x$n_patients, "\n")
  cat("  Genes used  :", length(x$genes_used),
      paste0("(", x$coverage_pct, "% coverage)"), "\n")
  cat("  Cutoff      :", round(x$threshold, 4),
      paste0("[", x$cutoff_method, "]"), "\n")
  cat("  Risk groups : Low =",  sum(x$risk_group == "Low"),
      "| High =", sum(x$risk_group == "High"), "\n")
  if (length(x$genes_missing) > 0)
    cat("  Missing genes:", paste(x$genes_missing, collapse = ", "), "\n")
  cat("=============================================\n")
  invisible(x)
}
