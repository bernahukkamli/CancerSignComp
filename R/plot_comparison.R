#' Plot C-index comparison across multiple signatures
#'
#' @param comparison A \code{SignatureComparison} object from
#'   \code{compare_signatures()}
#' @param show_se Logical: show standard error bars (default: TRUE)
#' @param title Optional plot title
#' @param color Bar color (default: "#2E86AB")
#' @param ref_line Numeric: draw a reference line at this C-index value.
#'   Default is 0.5 (random classifier). Set to NULL to hide.
#'
#' @return A \code{ggplot} object
#' @export
#'
#' @examples
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
#' plot_comparison(comp)
plot_comparison <- function(comparison,
                            show_se  = TRUE,
                            title    = NULL,
                            color    = "#2E86AB",
                            ref_line = 0.5) {

  # --- Girdi kontrolleri ---
  if (!inherits(comparison, "SignatureComparison")) {
    stop("`comparison` must be a SignatureComparison object from compare_signatures()")
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Install with: install.packages('ggplot2')")
  }

  df <- comparison$summary_table

  # Significance yıldızları ekle
  df$sig_label <- ifelse(df$logrank_p < 0.001, "***",
                         ifelse(df$logrank_p < 0.01,  "**",
                                ifelse(df$logrank_p < 0.05,  "*", "ns")))

  df$label <- paste0(df$signature, "\n(n genes = ", df$n_genes, ")")

  # C-index'e göre sırala
  df$label <- factor(df$label,
                     levels = df$label[order(df$c_index, decreasing = FALSE)])

  # --- Plot ---
  plot_title <- if (is.null(title)) {
    paste0("C-index Comparison (n = ", comparison$n_patients, " patients)")
  } else {
    title
  }

  p <- ggplot2::ggplot(df,
                       ggplot2::aes(x = .data$c_index, y = .data$label)) +
    ggplot2::geom_col(fill = color, alpha = 0.85, width = 0.6) +
    ggplot2::geom_text(
      ggplot2::aes(label = paste0(round(.data$c_index, 3),
                                  " ", .data$sig_label)),
      hjust = -0.1, size = 3.5
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, max(df$c_index) * 1.25),
      breaks = seq(0, 1, by = 0.1)
    ) +
    ggplot2::labs(
      title = plot_title,
      x     = "C-index (higher = better)",
      y     = NULL
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title  = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 10)
    )

  # Referans çizgisi
  if (!is.null(ref_line)) {
    p <- p + ggplot2::geom_vline(
      xintercept = ref_line,
      linetype   = "dashed",
      color      = "red",
      linewidth  = 0.5
    )
  }

  # SE hata çubukları
  if (show_se && "c_index_se" %in% names(df)) {
    p <- p + ggplot2::geom_errorbarh(
      ggplot2::aes(
        xmin = .data$c_index - .data$c_index_se,
        xmax = .data$c_index + .data$c_index_se
      ),
      height = 0.2, color = "gray40"
    )
  }

  return(p)
}
