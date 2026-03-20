#' Load and validate a cancer prognostic gene signature
#'
#' @param genes Character vector of gene symbols (e.g., c("TP53", "BRCA1"))
#' @param name  A short name for this signature (e.g., "Smith2023_PAAD")
#' @param cancer_type Cancer type abbreviation (e.g., "PAAD", "COAD", "KIRC")
#' @param description Optional free-text description of the signature
#'
#' @return A list of class \code{CancerSignature} with validated gene symbols
#' @export
#'
#' @examples
#' sig <- load_signature(
#'   genes       = c("TP53", "KRAS", "CDKN2A"),
#'   name        = "MySignature_2024",
#'   cancer_type = "PAAD"
#' )
#' print(sig)
load_signature <- function(genes,
                           name,
                           cancer_type,
                           description = "") {

  # --- Girdi kontrolleri ---
  if (!is.character(genes)) {
    stop("`genes` must be a character vector (e.g., c('TP53', 'KRAS'))")
  }

  if (length(genes) < 2) {
    stop("`genes` must contain at least 2 gene symbols")
  }

  if (!is.character(name) || nchar(trimws(name)) == 0) {
    stop("`name` must be a non-empty string")
  }

  if (!is.character(cancer_type) || nchar(trimws(cancer_type)) == 0) {
    stop("`cancer_type` must be a non-empty string (e.g., 'PAAD')")
  }

  # --- Gen sembollerini temizle ---
  genes_clean <- trimws(toupper(unique(genes)))

  # Geçersiz karakterler var mı? (gen sembolleri sadece harf, rakam, tire içerir)
  invalid <- genes_clean[!grepl("^[A-Z0-9\\-\\.]+$", genes_clean)]
  if (length(invalid) > 0) {
    warning("Possibly invalid gene symbols removed: ",
            paste(invalid, collapse = ", "))
    genes_clean <- setdiff(genes_clean, invalid)
  }

  if (length(genes_clean) < 2) {
    stop("Fewer than 2 valid gene symbols remain after cleaning.")
  }

  # --- CancerSignature objesi oluştur ---
  signature <- list(
    name        = trimws(name),
    cancer_type = trimws(toupper(cancer_type)),
    genes       = genes_clean,
    n_genes     = length(genes_clean),
    description = description,
    created_at  = Sys.time()
  )

  class(signature) <- "CancerSignature"
  return(signature)
}


#' Print method for CancerSignature objects
#'
#' @param x A \code{CancerSignature} object
#' @param ... Further arguments (ignored)
#' @export
print.CancerSignature <- function(x, ...) {
  cat("== CancerSignature ==========================\n")
  cat("  Name        :", x$name, "\n")
  cat("  Cancer type :", x$cancer_type, "\n")
  cat("  Genes (n =", x$n_genes, "):",
      paste(head(x$genes, 6), collapse = ", "),
      if (x$n_genes > 6) "..." else "", "\n")
  if (nchar(x$description) > 0)
    cat("  Description :", x$description, "\n")
  cat("=============================================\n")
  invisible(x)
}
