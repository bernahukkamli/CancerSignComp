# Sheffer_2009 - kolorektal kanser için iyi bilinen imza
sig_sheffer <- load_signature(
  genes  = c("CXCL10", "MMP3", "CXCL9", "IDO1", "GBP1", "STAT1", "IRF1"),
  name   = "Sheffer_2009",
  cancer = "COAD",
  source = "Sheffer et al. PNAS 2009"
)

cat("Genes found:", sum(sig_sheffer$genes %in% rownames(expr_sym_coad)),
    "/", sig_sheffer$n_genes, "\n")

scored_sheffer <- score_signature(
  sig_sheffer,
  expr_sym_coad,
  cutoff      = "median",
  norm_method = "none"
)

result_sheffer <- compare_to_staging(
  scored     = scored_sheffer,
  surv_time  = surv_time_coad,
  surv_event = surv_event_coad,
  stage      = stage_coad
)
print(result_sheffer)
