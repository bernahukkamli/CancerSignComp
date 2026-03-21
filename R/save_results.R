#' Save CancerSignComp results to CSV or Excel
#'
#' @param x A \code{SignatureComparison} or \code{StagingComparison} object
#' @param file Output file path. Extension determines format:
#'   \code{.csv} for CSV, \code{.xlsx} for Excel.
#' @param overwrite Logical: overwrite existing file? (default: FALSE)
#'
#' @return Invisibly returns the saved data frame
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' mat <- matrix(rnorm(1000), nrow = 10,
#'               dimnames = list(
#'                 c("TP53","KRAS","MYC","CDKN2A","SMAD4",
#'                   "EGFR","PTEN","RB1","BRCA2","PIK3CA"),
#'                 paste0("Patient_", 1:100)))
#' sig1   <- load_signature(c("TP53","KRAS","MYC"),  "Sig_A", "PAAD")
#' sig2   <- load_signature(c("EGFR","PTEN","RB1"),  "Sig_B", "PAAD")
#' surv_t <- rexp(100, rate = 0.05)
#' surv_e <- rbinom(100, 1, 0.7)
#' comp   <- compare_signatures(list(sig1, sig2), mat, surv_t, surv_e)
#' save_results(comp, "results.csv")
#' }
save_results <- function(x, file, overwrite = FALSE) {

  if (!inherits(x, c("SignatureComparison", "StagingComparison"))) {
    stop("`x` must be a SignatureComparison or StagingComparison object")
  }

  if (file.exists(file) && !overwrite) {
    stop("File already exists. Use overwrite = TRUE to replace.")
  }

  ext <- tolower(tools::file_ext(file))

  if (!ext %in% c("csv", "xlsx")) {
    stop("`file` must end in .csv or .xlsx")
  }

  # Veriyi hazirla
  if (inherits(x, "SignatureComparison")) {
    df <- x$summary_table
  } else {
    df <- data.frame(
      signature      = x$signature_name,
      cancer_type    = x$cancer_type,
      n_patients     = x$n_patients_total,
      n_staged       = x$n_patients_staged,
      sig_c_index    = x$sig_c_index,
      sig_c_index_se = x$sig_c_index_se,
      sig_logrank_p  = x$sig_logrank_p,
      stage_c_index  = x$stage_c_index,
      stage_c_index_se = x$stage_c_index_se,
      stage_logrank_p  = x$stage_logrank_p,
      c_index_diff   = x$c_index_diff,
      z_score        = x$z_score,
      p_value_diff   = x$p_value_diff,
      verdict        = x$verdict,
      nri            = if (!is.null(x$nri)) x$nri else NA,
      nri_p          = if (!is.null(x$nri_p)) x$nri_p else NA,
      idi            = if (!is.null(x$idi)) x$idi else NA,
      idi_p          = if (!is.null(x$idi_p)) x$idi_p else NA,
      stringsAsFactors = FALSE
    )
  }

  # Kaydet
  if (ext == "csv") {
    utils::write.csv(df, file, row.names = FALSE)
  } else {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      stop("Package 'openxlsx' is required for Excel export. ",
           "Install with: install.packages('openxlsx')")
    }
    openxlsx::write.xlsx(df, file, overwrite = overwrite)
  }

  message("Results saved to: ", file)
  invisible(df)
}
