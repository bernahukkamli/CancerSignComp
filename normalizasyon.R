# compare_normalization.R dosyasını oluştur
cat('
#\' Compare normalization methods for signature scoring
#\'
#\' @param signature A SignatureSet object from load_signature()
#\' @param expr_matrix Gene expression matrix (genes x samples)
#\' @param methods Character vector of normalization methods to compare
#\' @param weights Optional named numeric vector of gene weights
#\' @param cutoff Cutoff method for high/low groups
#\' @return A NormComparison S3 object
#\' @export
compare_normalization <- function(signature, expr_matrix,
                                   methods = c("none", "log2", "zscore", "quantile"),
                                   weights = NULL,
                                   cutoff = "median") {
  results <- list()
  for (m in methods) {
    scored <- tryCatch(
      score_signature(signature, expr_matrix, weights = weights,
                      cutoff = cutoff, norm_method = m),
      error = function(e) NULL
    )
    if (!is.null(scored)) {
      results[[m]] <- scored
    }
  }
  structure(
    list(results = results, methods = names(results), signature = signature),
    class = "NormComparison"
  )
}

#\' @export
print.NormComparison <- function(x, ...) {
  cat("NormComparison object\\n")
  cat("Signature:", x$signature$name, "\\n")
  cat("Methods tested:", paste(x$methods, collapse = ", "), "\\n")
  for (m in x$methods) {
    r <- x$results[[m]]
    cat(sprintf("  %s: C-index = %.3f\\n", m, r$cindex))
  }
  invisible(x)
}

#\' @export
plot.NormComparison <- function(x, ...) {
  cindex_vals <- sapply(x$methods, function(m) x$results[[m]]$cindex)
  df <- data.frame(method = names(cindex_vals), cindex = cindex_vals)
  df$method <- factor(df$method, levels = df$method[order(df$cindex)])
  ggplot2::ggplot(df, ggplot2::aes(x = method, y = cindex)) +
    ggplot2::geom_col(fill = "steelblue") +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
    ggplot2::coord_flip() +
    ggplot2::labs(title = paste("Normalization comparison:", x$signature$name),
                  x = "Method", y = "C-index") +
    ggplot2::theme_minimal()
}
', file = "C:/Users/berna/CancerSignComp/R/compare_normalization.R")

cat("✓ Dosya oluşturuldu\n")
