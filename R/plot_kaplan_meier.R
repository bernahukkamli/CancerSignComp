#' Plot Kaplan-Meier survival curves for a scored signature
#'
#' @param scored A \code{ScoredSignature} object from \code{score_signature()}
#' @param surv_time Numeric vector of survival times
#' @param surv_event Numeric vector: 1 = event, 0 = censored
#' @param time_unit Label for the x-axis time unit (default: "Months")
#' @param title Optional plot title. If NULL, uses signature name.
#' @param palette Color palette: "jco" (default), "lancet", "nejm", "aaas"
#' @param risk_table Logical: show risk table below plot (default: TRUE)
#' @param conf_int Logical: show confidence intervals (default: TRUE)
#'
#' @return A \code{ggsurvplot} object (printable and saveable)
#' @export
#'
#' @examples
#' \dontrun{
#'   set.seed(42)
#'   mat <- matrix(rnorm(500), nrow = 5,
#'                 dimnames = list(c("TP53","KRAS","MYC","CDKN2A","SMAD4"),
#'                                 paste0("Patient_", 1:100)))
#'   sig    <- load_signature(c("TP53","KRAS","MYC"), "TestSig", "PAAD")
#'   scored <- score_signature(sig, mat)
#'   surv_t <- rexp(100, rate = 0.05)
#'   surv_e <- rbinom(100, 1, 0.7)
#'   p <- plot_kaplan_meier(scored, surv_t, surv_e)
#'   print(p)
#' }
plot_kaplan_meier <- function(scored,
                              surv_time,
                              surv_event,
                              time_unit  = "Months",
                              title      = NULL,
                              palette    = "jco",
                              risk_table = TRUE,
                              conf_int   = TRUE) {
  if (!inherits(scored, "ScoredSignature")) {
    stop("`scored` must be a ScoredSignature object from score_signature()")
  }
  if (!requireNamespace("survminer", quietly = TRUE)) {
    stop("Package 'survminer' is required.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }
  n <- scored$n_patients
  if (length(surv_time) != n) stop("`surv_time` length must equal n_patients")
  if (length(surv_event) != n) stop("`surv_event` length must equal n_patients")
  df_surv <- data.frame(
    surv_time  = surv_time,
    surv_event = surv_event,
    risk_group = scored$risk_group
  )
  km_fit <- survival::survfit(
    survival::Surv(surv_time, surv_event) ~ risk_group,
    data = df_surv
  )
  km_diff <- survival::survdiff(
    survival::Surv(surv_time, surv_event) ~ risk_group,
    data = df_surv
  )
  logrank_p <- 1 - pchisq(km_diff$chisq, df = 1)
  plot_title <- if (is.null(title)) {
    paste0(scored$signature_name, " (", scored$cancer_type, ")")
  } else {
    title
  }
  p_label <- if (logrank_p < 0.001) {
    "p < 0.001"
  } else {
    paste0("p = ", round(logrank_p, 3))
  }
  p <- survminer::ggsurvplot(
    fit               = km_fit,
    data              = df_surv,
    risk.table        = risk_table,
    conf.int          = conf_int,
    palette           = palette,
    legend.labs       = c("Low risk", "High risk"),
    legend.title      = "Risk group",
    title             = plot_title,
    xlab              = paste0("Time (", time_unit, ")"),
    ylab              = "Survival probability",
    pval              = p_label,
    pval.coord        = c(max(surv_time) * 0.6, 0.9),
    risk.table.height = if (risk_table) 0.30 else 0,
    risk.table.y.text = FALSE,
    risk.table.col    = "strata",
    fontsize          = 3.5,
    ggtheme           = ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        plot.title      = ggplot2::element_text(hjust = 0.5, face = "bold"),
        legend.position = "top"
      )
  )
  return(p)
}
