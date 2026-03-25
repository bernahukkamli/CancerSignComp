lines <- readLines("C:/Users/berna/CancerSignComp/R/compare_normalization.R")
# survcomp satırını bul ve değiştir
lines <- gsub(
  "survcomp::concordance.index\\(scored\\$risk_scores\\[common\\], surv_obj\\)\\$c.index",
  "survival::concordance(surv_obj ~ scored$risk_scores[common])$concordance",
  lines
)
writeLines(lines, "C:/Users/berna/CancerSignComp/R/compare_normalization.R")

devtools::load_all()

# Tekrar test et
norm_comparison <- compare_normalization(
  signature      = sig_zhou2022,
  expr_matrix    = expr_sym_luad,
  weights        = weights_zhou2022,
  methods        = c("none", "log2", "zscore", "quantile"),
  cutoff         = "median",
  survival_time  = surv_time_luad,
  survival_event = surv_event_luad
)

print(norm_comparison)
plot(norm_comparison)
