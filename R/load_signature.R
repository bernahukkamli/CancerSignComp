#' Load and validate a cancer prognostic gene signature
#'
#' @param genes Character vector of gene symbols
#' @param name A short name for this signature
#' @param cancer_type Cancer type abbreviation (e.g., "PAAD")
#' @param cancer Alias for cancer_type
#' @param description Optional free-text description
#' @param source Optional citation or reference
#' @param weights Optional named numeric vector of gene weights
#' @return A list of class \code{CancerSignature}
#' @export
#' @examples
#' sig <- load_signature(genes = c("TP53", "KRAS"), name = "Test", cancer_type = "PAAD")
load_signature <- function(genes,
                           name,
                           cancer_type = NULL,
                           cancer      = NULL,
                           description = "",
                           source      = "",
                           weights     = NULL) {
  # cancer alias
  if (is.null(cancer_type) && !is.null(cancer)) cancer_type <- cancer
  if (is.null(cancer_type)) stop("`cancer_type` must be provided")
  if (!is.character(genes)) stop("`genes` must be a character vector")
  if (length(genes) < 2) stop("`genes` must contain at least 2 gene symbols")
  if (!is.character(name) || nchar(trimws(name)) == 0) stop("`name` must be a non-empty string")
  if (!is.null(weights) && !is.numeric(weights)) stop("`weights` must be a named numeric vector")
  if (!is.null(weights) && is.null(names(weights))) stop("`weights` must have names (gene symbols)")
  genes_clean <- trimws(toupper(unique(genes)))
  invalid <- genes_clean[!grepl("^[A-Z0-9\\-\\.]+$", genes_clean)]
  if (length(invalid) > 0) {
    warning("Possibly invalid gene symbols removed: ", paste(invalid, collapse = ", "))
    genes_clean <- setdiff(genes_clean, invalid)
  }
  if (length(genes_clean) < 2) stop("Fewer than 2 valid gene symbols remain after cleaning.")
  signature <- list(
    name        = trimws(name),
    cancer_type = trimws(toupper(cancer_type)),
    genes       = genes_clean,
    n_genes     = length(genes_clean),
    description = description,
    source      = source,
    weights     = weights,
    created_at  = Sys.time()
  )
  class(signature) <- "CancerSignature"
  return(signature)
}

#' @export
print.CancerSignature <- function(x, ...) {
  cat("== CancerSignature ==========================\n")
  cat("  Name        :", x$name, "\n")
  cat("  Cancer type :", x$cancer_type, "\n")
  cat("  Genes (n =", x$n_genes, "):",
      paste(head(x$genes, 6), collapse = ", "),
      if (x$n_genes > 6) "..." else "", "\n")
  if (!is.null(x$source) && nchar(x$source) > 0) cat("  Source      :", x$source, "\n")
  if (nchar(x$description) > 0) cat("  Description :", x$description, "\n")
  if (!is.null(x$weights)) cat("  Weights     : provided (n =", length(x$weights), ")\n")
  cat("=============================================\n")
  invisible(x)
}
