
#' Validate a signature in an external cohort
#'
#' @param signature A CancerSignature object from load_signature()
#' @param expr_matrix A numeric matrix: rows = genes, columns = patients
#' @param surv_time Numeric vector of survival times
#' @param surv_event Numeric vector of survival events (1=event, 0=censored)
#' @param stage Optional numeric vector of AJCC/TNM stage
#' @param cohort_name Character string identifying the external cohort (e.g. "GSE17536")
#' @return A list of class ExternalValidation with results
#' @export
validate_external <- function(signature, expr_matrix, surv_time, surv_event,
                               stage = NULL, cohort_name = "External") {

    scored <- score_signature(signature, expr_matrix)

    if (!is.null(stage)) {
        result     <- compare_to_staging(
            scored     = scored,
            surv_time  = surv_time,
            surv_event = surv_event,
            stage      = stage
        )
        c_index    <- result$sig_c_index
        logrank_p  <- result$sig_logrank_p
        delta_c    <- result$c_index_diff
        verdict    <- result$verdict
        tnm_cindex <- result$stage_c_index
    } else {
        surv_obj   <- survival::Surv(surv_time, surv_event)
        cox_fit    <- survival::coxph(surv_obj ~ scored$risk_scores)
        c_index    <- summary(cox_fit)$concordance[1]
        lr         <- survival::survdiff(surv_obj ~ scored$risk_group)
        logrank_p  <- 1 - stats::pchisq(lr$chisq, df = 1)
        delta_c    <- NA
        verdict    <- "No staging available"
        tnm_cindex <- NA
    }

    out <- list(
        signature_name = signature$name,
        cancer_type    = signature$cancer_type,
        cohort         = cohort_name,
        n_patients     = ncol(expr_matrix),
        c_index        = c_index,
        tnm_cindex     = tnm_cindex,
        logrank_p      = logrank_p,
        delta_c        = delta_c,
        verdict        = verdict,
        scored         = scored
    )
    class(out) <- "ExternalValidation"

    cat("== ExternalValidation ==========================
")
    cat("  Signature  :", signature$name, "
")
    cat("  Cohort     :", cohort_name, "
")
    cat("  Patients   :", ncol(expr_matrix), "
")
    cat("  C-index    :", round(c_index, 4), "
")
    cat("  TNM C-index:", ifelse(is.na(tnm_cindex), "N/A",
                                 as.character(round(tnm_cindex, 4))), "
")
    cat("  Log-rank p :", round(logrank_p, 4), "
")
    if (!is.na(delta_c))
        cat("  Delta C    :", round(delta_c, 4), "
")
    cat("  Verdict    :", verdict, "
")
    cat("================================================
")

    invisible(out)
}

