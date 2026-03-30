# Guinney_2015 CMS (Consensus Molecular Subtypes) - en iyi bilinen COAD imzası
# Ama daha basit: Kemper_2012 - 7-gen COAD prognostik imzası
sig_kemper <- load_signature(
  genes  = c("MACC1", "MET", "EGFR", "VEGFA", "VEGFB", "VEGFC", "VEGFD"),
  name   = "Kemper_2012",
  cancer = "COAD",
  source = "Kemper et al. Ann Surg 2012"
)

# Coverage kontrol
cat("Genes found:", sum(sig_kemper$genes %in% rownames(expr_sym_coad)),
    "/", sig_kemper$n_genes, "\n")

# Score et
scored_kemper <- score_signature(
  sig_kemper,
  expr_sym_coad,
  cutoff      = "median",
  norm_method = "none"
)

cat("Coverage:", scored_kemper$coverage_pct, "%\n")

# C-index
common <- intersect(names(scored_kemper$risk_scores), names(surv_time_coad))
ci_kemper <- survival::concordance(
  survival::Surv(surv_time_coad[common], surv_event_coad[common]) ~
    scored_kemper$risk_scores[common]
)
cat("Kemper_2012 C-index:", ci_kemper$concordance, "\n")
