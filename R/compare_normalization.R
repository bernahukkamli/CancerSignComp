
#' Compare normalization methods for signature scoring
#'
#' @param signature A SignatureSet object from load_signature()
#' @param expr_matrix Gene expression matrix (genes x samples)
#' @param methods Character vector of normalization methods to compare
#' @param weights Optional named numeric vector of gene weights
#' @param cutoff Cutoff method for high/low groups
#' @param survival_time Numeric vector of survival times
#' @param survival_event Numeric vector of survival events (0/1)
#' @return A NormComparison S3 object
#' @export
compare_normalization <- function(signature, expr_matrix,
                                   methods = c("none", "log2", "zscore", "quantile"),
                                   weights = NULL,
                                   cutoff = "median",
                                   survival_time = NULL,
                                   survival_event = NULL) {
  results <- list()
  cindex_vals <- c()
  
  for (m in methods) {
    scored <- tryCatch(
      score_signature(signature, expr_matrix, weights = weights,
                      cutoff = cutoff, norm_method = m),
      error = function(e) NULL
    )
    if (!is.null(scored)) {
      results[[m]] <- scored
      # C-index hesapla (survival bilgisi varsa)
      if (!is.null(survival_time) && !is.null(survival_event)) {
        common <- intersect(names(scored$risk_scores), names(survival_time))
        if (length(common) > 10) {
          ci <- tryCatch({
            surv_obj <- survival::Surv(survival_time[common], survival_event[common])
            survival::concordance(surv_obj ~ scored$risk_scores[common])$concordance
          }, error = function(e) NA)
          cindex_vals[m] <- ci
        } else {
          cindex_vals[m] <- NA
        }
      }
    }
  }
  
  structure(
    list(results = results, methods = names(results), 
         cindex = cindex_vals, signature = signature),
    class = "NormComparison"
  )
}

#' @export
print.NormComparison <- function(x, ...) {
  cat("NormComparison object\n")
  cat("Signature:", x$signature$name, "\n")
  cat("Methods tested:", paste(x$methods, collapse = ", "), "\n")
  if (length(x$cindex) > 0) {
    cat("C-index values:\n")
    for (m in names(x$cindex)) {
      cat(sprintf("  %-10s: %.3f\n", m, x$cindex[m]))
    }
  } else {
    cat("(Provide survival_time and survival_event for C-index calculation)\n")
    for (m in x$methods) {
      r <- x$results[[m]]
      cat(sprintf("  %-10s: n=%d, norm=%s\n", m, r$n_patients, r$norm_method))
    }
  }
  invisible(x)
}

#' @export
plot.NormComparison <- function(x, ...) {
  if (length(x$cindex) == 0 || all(is.na(x$cindex))) {
    stop("No C-index values. Provide survival_time and survival_event.")
  }
  df <- data.frame(
    method = names(x$cindex), 
    cindex = as.numeric(x$cindex),
    stringsAsFactors = FALSE
  )
  df <- df[!is.na(df$cindex), ]
  df$method <- factor(df$method, levels = df$method[order(df$cindex)])
  ggplot2::ggplot(df, ggplot2::aes(x = method, y = cindex)) +
    ggplot2::geom_col(fill = "steelblue") +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
    ggplot2::coord_flip() +
    ggplot2::labs(title = paste("Normalization comparison:", x$signature$name),
                  x = "Method", y = "C-index") +
    ggplot2::theme_minimal()
}
