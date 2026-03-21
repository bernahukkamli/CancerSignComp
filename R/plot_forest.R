#' Forest plot for subgroup survival analysis
#'
#' @param expr_matrix Numeric matrix: rows = genes, columns = patients
#' @param surv_time Numeric vector of survival times
#' @param surv_event Numeric vector: 1 = event, 0 = censored
#' @param signature A \code{CancerSignature} object
#' @param subgroups Named list of patient index vectors defining subgroups.
#'   Example: list(Stage_I_II = c(1,5,10), Stage_III_IV = c(2,3,8))
#' @param subgroup_var Optional character vector (same length as patients)
#'   for automatic subgroup splitting. If provided, \code{subgroups} is ignored.
#' @param title Optional plot title
#'
#' @return A \code{ggplot} object
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#'   set.seed(42)
#'   mat <- matrix(rnorm(500), nrow = 5,
#'                 dimnames = list(c("TP53","KRAS","MYC","CDKN2A","SMAD4"),
#'                                 paste0("Patient_", 1:100)))
#'   sig       <- load_signature(c("TP53","KRAS","MYC"), "TestSig", "PAAD")
#'   surv_t    <- rexp(100, rate = 0.05)
#'   surv_e    <- rbinom(100, 1, 0.7)
#'   age_group <- sample(c("Young","Old"), 100, replace = TRUE)
#'   plot_forest(mat, surv_t, surv_e, sig, subgroup_var = age_group)
#' }
plot_forest <- function(expr_matrix,
                        surv_time,
                        surv_event,
                        signature,
                        subgroups    = NULL,
                        subgroup_var = NULL,
                        title        = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }
  if (!inherits(signature, "CancerSignature")) {
    stop("`signature` must be a CancerSignature object")
  }

  n <- ncol(expr_matrix)

  if (!is.null(subgroup_var)) {
    if (length(subgroup_var) != n) {
      stop("`subgroup_var` length must equal number of patients")
    }
    subgroups <- split(seq_len(n), subgroup_var)
  }

  if (is.null(subgroups)) {
    stop("Provide either `subgroups` or `subgroup_var`")
  }

  # Her alt grup icin HR hesapla
  results <- lapply(names(subgroups), function(sg_name) {
    idx <- subgroups[[sg_name]]
    if (length(idx) < 10) {
      warning("Subgroup '", sg_name, "' has fewer than 10 patients - skipping")
      return(NULL)
    }
    sub_mat   <- expr_matrix[, idx, drop = FALSE]
    sub_time  <- surv_time[idx]
    sub_event <- surv_event[idx]
    scored    <- tryCatch(score_signature(signature, sub_mat),
                          error = function(e) NULL)
    if (is.null(scored)) return(NULL)
    surv_obj <- survival::Surv(sub_time, sub_event)
    cox_fit  <- survival::coxph(surv_obj ~ scored$risk_scores)
    cox_sum  <- summary(cox_fit)
    data.frame(
      subgroup = sg_name,
      n        = length(idx),
      hr       = round(cox_sum$conf.int[1, "exp(coef)"],  3),
      hr_lo    = round(cox_sum$conf.int[1, "lower .95"],  3),
      hr_hi    = round(cox_sum$conf.int[1, "upper .95"],  3),
      p_val    = round(cox_sum$coefficients[1, "Pr(>|z|)"], 4),
      stringsAsFactors = FALSE
    )
  })

  # Overall
  scored_all <- score_signature(signature, expr_matrix)
  surv_all   <- survival::Surv(surv_time, surv_event)
  cox_all    <- survival::coxph(surv_all ~ scored_all$risk_scores)
  cox_all_s  <- summary(cox_all)

  overall <- data.frame(
    subgroup = "Overall",
    n        = n,
    hr       = round(cox_all_s$conf.int[1, "exp(coef)"],  3),
    hr_lo    = round(cox_all_s$conf.int[1, "lower .95"],  3),
    hr_hi    = round(cox_all_s$conf.int[1, "upper .95"],  3),
    p_val    = round(cox_all_s$coefficients[1, "Pr(>|z|)"], 4),
    stringsAsFactors = FALSE
  )

  results_clean <- Filter(Negate(is.null), results)
  df <- do.call(rbind, c(results_clean, list(overall)))

  df$sig_label <- ifelse(df$p_val < 0.001, "***",
                         ifelse(df$p_val < 0.01,  "**",
                                ifelse(df$p_val < 0.05,  "*", "ns")))
  df$hr_label  <- paste0("HR=", sprintf("%.2f", df$hr), " ", df$sig_label)
  df$label     <- paste0(df$subgroup, " (n=", df$n, ")")
  df$label     <- factor(df$label, levels = rev(df$label))

  x_max <- max(df$hr_hi, na.rm = TRUE)
  x_min <- min(df$hr_lo, na.rm = TRUE)

  plot_title <- if (is.null(title)) {
    paste0("Forest Plot: ", signature$name)
  } else {
    title
  }

  ggplot2::ggplot(df, ggplot2::aes(x = .data$hr, y = .data$label)) +
    ggplot2::geom_vline(
      xintercept = 1, linetype = "dashed",
      color = "red", linewidth = 0.5
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = .data$hr_lo, xmax = .data$hr_hi),
      width = 0.25, color = "gray40", linewidth = 0.7,
      orientation = "y"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(size = .data$n),
      color = "#2E86AB", alpha = 0.85
    ) +
    ggplot2::annotate(
      "text",
      x     = x_max * 1.08,
      y     = df$label,
      label = df$hr_label,
      hjust = 0,
      size  = 3.2
    ) +
    ggplot2::scale_size_continuous(range = c(3, 8), guide = "none") +
    ggplot2::scale_x_continuous(
      limits = c(x_min * 0.85, x_max * 1.5),
      breaks = seq(round(x_min * 0.85, 1), round(x_max * 1.5, 1), by = 0.1)
    ) +
    ggplot2::labs(
      title = plot_title,
      x     = "Hazard Ratio (95% CI)",
      y     = NULL
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title  = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 10)
    )
}
