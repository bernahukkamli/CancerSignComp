# Zhou et al. 2022 imzası
sig_zhou2022 <- load_signature(
  genes       = c("PLEK2", "COL1A1", "GPX3"),
  name        = "Zhou_2022",
  cancer_type = "LUAD",
  description = "Zhou et al. Transl Lung Cancer Res 2022 - 3-gene signature"
)

weights_zhou2022 <- c(
  COL1A1 =  0.1446053,
  GPX3   = -0.2426827,
  PLEK2  =  0.2697514
)

# Orijinal katsayılarla skorla
scored_zhou <- score_signature(
  sig_zhou2022,
  expr_final_luad,
  weights = weights_zhou2022,
  method  = "weighted_mean"
)
print(scored_zhou)

# C-index
surv_obj_luad <- Surv(surv_t_luad, surv_e_luad)
c_zhou <- concordance(surv_obj_luad ~ scored_zhou$risk_scores)
cat("\nZhou_2022 C-index (bizim)    :",
    round(max(c_zhou$concordance, 1 - c_zhou$concordance), 4), "\n")
cat("Zhou_2022 C-index (makalede) : 0.6375\n")

# Staging ile karşılaştır
result_zhou <- compare_to_staging(
  scored_zhou, surv_t_luad, surv_e_luad, stage_luad_clean
)
print(result_zhou)

# KM grafiği
km_zhou <- plot_kaplan_meier(
  scored_zhou, surv_t_luad, surv_e_luad,
  time_unit = "Days",
  title     = "Zhou_2022 (PLEK2+COL1A1+GPX3) — TCGA-LUAD"
)
print(km_zhou)
