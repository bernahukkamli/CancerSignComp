#' Compare a gene signature against clinical staging (TNM/AJCC)
#'
#' @param scored A \code{ScoredSignature} object from \code{score_signature()}
#' @param surv_time Numeric vector of survival times
#' @param surv_event Numeric vector: 1 = event, 0 = censored
#' @param stage Factor or character vector of clinical stages for each patient.
#'   Accepted formats: "I","II","III","IV" or "Stage I","Stage II", etc.
#'   Numeric values 1-4 are also accepted.
#' @param stage_grouping How to group stages for survival analysis:
#'   \code{"auto"} (default): I+II = Early, III+IV = Late
#'   \code{"each"}: treat each stage separately
#'
#' @return A list of class \code{StagingComparison}
#' @export
#'
#' @examples
#' set.seed(42)
#' mat <- matrix(rnorm(500), nrow = 5,
#'               dimnames = list(c("TP53","KRAS","MYC","CDKN2A","SMAD4"),
#'                               paste0("Patient_", 1:100)))
#' sig    <- load_signature(c("TP53","KRAS","MYC"), "TestSig", "PAAD")
#' scored <- score_signature(sig, mat)
#' surv_t <- rexp(100, rate = 0.05)
#' surv_e <- rbinom(100, 1, 0.7)
#' stages <- sample(c("I","II","III","IV"), 100, replace = TRUE)
#' result <- compare_to_staging(scored, surv_t, surv_e, stages)
#' print(result)
compare_to_staging <- function(scored,
                               surv_time,
                               surv_event,
                               stage,
                               stage_grouping = "auto") {

  if (!inherits(scored, "ScoredSignature")) {
    stop("`scored` must be a ScoredSignature object from score_signature()")
  }
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required.")
  }

  n <- scored$n_patients

  # --- Sample eslestirme (isimli vektorler icin) ---
  patient_ids <- names(scored$risk_scores)

  # surv_time eslestir
  if (!is.null(names(surv_time))) {
    common_t <- intersect(patient_ids, names(surv_time))
    if (length(common_t) < 10) stop("Fewer than 10 common samples between scored and surv_time")
    surv_time  <- surv_time[common_t]
    surv_event <- if (!is.null(names(surv_event))) surv_event[common_t] else surv_event[seq_along(common_t)]
    scored$risk_scores <- scored$risk_scores[common_t]
    scored$risk_group  <- scored$risk_group[common_t]
    n <- length(common_t)
  }

  # stage eslestir
  if (!is.null(names(stage))) {
    common_s <- intersect(names(scored$risk_scores), names(stage))
    if (length(common_s) < 10) stop("Fewer than 10 common samples between scored and stage")
    stage      <- stage[common_s]
    surv_time  <- surv_time[common_s]
    surv_event <- surv_event[common_s]
    scored$risk_scores <- scored$risk_scores[common_s]
    scored$risk_group  <- scored$risk_group[common_s]
    n <- length(common_s)
  }

  if (length(surv_time)  != n) stop("`surv_time` length must equal n_patients")
  if (length(surv_event) != n) stop("`surv_event` length must equal n_patients")
  if (!all(surv_event %in% c(0, 1))) stop("`surv_event` must contain only 0 and 1")
  if (length(stage) != n) stop("`stage` length must equal n_patients")

  # --- Stage temizle ---
  stage_clean <- toupper(trimws(as.character(stage)))
  stage_clean <- gsub("STAGE\\s*", "", stage_clean)
  stage_clean <- gsub("^1$", "I",   stage_clean)
  stage_clean <- gsub("^2$", "II",  stage_clean)
  stage_clean <- gsub("^3$", "III", stage_clean)
  stage_clean <- gsub("^4$", "IV",  stage_clean)
  stage_clean <- gsub("^(I{1,3}V?).*", "\\1", stage_clean)

  valid_stages <- c("I", "II", "III", "IV")
  unknown      <- !stage_clean %in% valid_stages

  if (any(unknown)) {
    warning(sum(unknown), " patient(s) have unrecognized stage values: ",
            paste(unique(stage_clean[unknown]), collapse = ", "),
            ". These will be excluded from staging analysis.")
  }

  keep       <- !unknown
  stage_use  <- stage_clean[keep]
  surv_t_use <- surv_time[keep]
  surv_e_use <- surv_event[keep]
  scores_use <- scored$risk_scores[keep]

  surv_obj_full    <- survival::Surv(surv_time, surv_event)
  surv_obj_staging <- survival::Surv(surv_t_use, surv_e_use)

  # --- İmza metrikleri ---
  cfit_sig   <- survival::concordance(surv_obj_full ~ scored$risk_scores)
  c_raw      <- cfit_sig$concordance
  c_sig      <- round(max(c_raw, 1 - c_raw), 4)
  c_sig_se   <- round(sqrt(cfit_sig$var), 4)

  km_sig     <- survival::survdiff(surv_obj_full ~ scored$risk_group)
  p_sig      <- round(1 - pchisq(km_sig$chisq, df = 1), 4)

  # --- Staging metrikleri ---
  stage_numeric <- as.numeric(factor(stage_use, levels = c("I","II","III","IV")))

  if (stage_grouping == "auto") {
    stage_group <- ifelse(stage_use %in% c("I", "II"), "Early", "Late")
    stage_group <- factor(stage_group, levels = c("Early", "Late"))
  } else {
    stage_group <- factor(stage_use, levels = valid_stages)
  }

  cfit_stage <- survival::concordance(surv_obj_staging ~ stage_numeric)
  c_raw_s    <- cfit_stage$concordance
  c_stage    <- round(max(c_raw_s, 1 - c_raw_s), 4)
  c_stage_se <- round(sqrt(cfit_stage$var), 4)

  km_stage   <- survival::survdiff(surv_obj_staging ~ stage_group)
  df_stage   <- length(levels(stage_group)) - 1
  p_stage    <- round(1 - pchisq(km_stage$chisq, df = df_stage), 4)

  # --- C-index farki ---
  c_diff    <- round(c_sig - c_stage, 4)
  c_diff_se <- round(sqrt(cfit_sig$var + cfit_stage$var), 4)
  z_score   <- round(c_diff / c_diff_se, 3)
  p_diff    <- round(2 * pnorm(-abs(z_score)), 4)

  verdict <- if (c_diff > 0.05 && p_diff < 0.05) {
    "Signature significantly OUTPERFORMS clinical staging"
  } else if (c_diff < -0.05 && p_diff < 0.05) {
    "Signature significantly UNDERPERFORMS clinical staging"
  } else {
    "Signature performs SIMILARLY to clinical staging"
  }

  # --- NRI / IDI ---
  nri_result <- tryCatch({
    sig_prob   <- (scores_use - min(scores_use)) /
      (max(scores_use) - min(scores_use))
    stage_prob <- (stage_numeric - min(stage_numeric)) /
      (max(stage_numeric) - min(stage_numeric))

    df_nri <- data.frame(
      surv_t_use  = surv_t_use,
      surv_e_use  = surv_e_use,
      stage_prob  = stage_prob,
      sig_prob    = sig_prob
    )

    nri_obj <- nricens::nribin(
      mdl.std = survival::coxph(
        survival::Surv(surv_t_use, surv_e_use) ~ stage_prob,
        data = df_nri
      ),
      mdl.new = survival::coxph(
        survival::Surv(surv_t_use, surv_e_use) ~ sig_prob,
        data = df_nri
      ),
      updown = "category",
      cut    = median(sig_prob),
      niter  = 100,
      alpha  = 0.05
    )
    list(
      nri       = round(nri_obj$nri["NRI", "Estimate"], 4),
      nri_p     = round(nri_obj$nri["NRI", "p-value"], 4),
      idi       = round(nri_obj$idi["IDI", "Estimate"], 4),
      idi_p     = round(nri_obj$idi["IDI", "p-value"], 4),
      available = TRUE
    )
  }, error = function(e) {
    list(nri = NA, nri_p = NA, idi = NA, idi_p = NA, available = FALSE)
  })

  # --- Sonuc ---
  result <- list(
    signature_name     = scored$signature_name,
    cancer_type        = scored$cancer_type,
    n_patients_total   = n,
    n_patients_staged  = sum(keep),
    stage_grouping     = stage_grouping,
    sig_c_index        = c_sig,
    sig_c_index_se     = c_sig_se,
    sig_logrank_p      = p_sig,
    stage_c_index      = c_stage,
    stage_c_index_se   = c_stage_se,
    stage_logrank_p    = p_stage,
    c_index_diff       = c_diff,
    c_index_diff_se    = c_diff_se,
    z_score            = z_score,
    p_value_diff       = p_diff,
    verdict            = verdict,
    nri                = nri_result$nri,
    nri_p              = nri_result$nri_p,
    idi                = nri_result$idi,
    idi_p              = nri_result$idi_p,
    nri_available      = nri_result$available,
    stage_distribution = table(stage_use)
  )

  class(result) <- "StagingComparison"
  return(result)
}


#' Print method for StagingComparison objects
#'
#' @param x A \code{StagingComparison} object
#' @param ... Further arguments (ignored)
#' @export
print.StagingComparison <- function(x, ...) {
  cat("== StagingComparison ========================\n")
  cat("  Signature   :", x$signature_name, "\n")
  cat("  Cancer type :", x$cancer_type, "\n")
  cat("  Patients    :", x$n_patients_total,
      "(staged:", x$n_patients_staged, ")\n\n")
  cat("  Gene signature:\n")
  cat(sprintf("    C-index   : %.4f (SE=%.4f)\n",
              x$sig_c_index, x$sig_c_index_se))
  cat(sprintf("    Log-rank p: %.4f\n", x$sig_logrank_p))
  cat("\n  Clinical staging (TNM/AJCC):\n")
  cat(sprintf("    C-index   : %.4f (SE=%.4f)\n",
              x$stage_c_index, x$stage_c_index_se))
  cat(sprintf("    Log-rank p: %.4f\n", x$stage_logrank_p))
  cat("\n  Comparison:\n")
  cat(sprintf("    Delta C   : %+.4f (z=%.3f, p=%.4f)\n",
              x$c_index_diff, x$z_score, x$p_value_diff))
  cat("   ", x$verdict, "\n")
  if (!is.null(x$nri_available) && x$nri_available) {
    cat(sprintf("    NRI       : %.4f (p=%.4f)\n", x$nri, x$nri_p))
    cat(sprintf("    IDI       : %.4f (p=%.4f)\n", x$idi, x$idi_p))
  }
  cat("\n  Stage distribution:\n  ")
  print(x$stage_distribution)
  cat("=============================================\n")
  invisible(x)
}
