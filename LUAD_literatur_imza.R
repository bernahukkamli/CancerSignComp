devtools::load_all("C:/Users/berna/CancerSignComp")

# Stage temizle
stage_luad_clean <- gsub("Stage ", "", stage_luad)
stage_luad_clean <- gsub("IA|IB", "I", stage_luad_clean)
stage_luad_clean <- gsub("IIA|IIB", "II", stage_luad_clean)
stage_luad_clean <- gsub("IIIA|IIIB", "III", stage_luad_clean)
table(stage_luad_clean)

# --- İmzalar ---
sig_chen2007 <- load_signature(
  genes       = c("DUSP6", "MMD", "STAT1", "ERBB3", "LCK"),
  name        = "Chen_2007",
  cancer_type = "LUAD",
  description = "Chen et al. NEJM 2007 - 5-gen klinik imza"
)

sig_li2023 <- load_signature(
  genes       = c("BTK", "FGFR2", "PIM2", "CHEK1", "CDK1"),
  name        = "Li_2023",
  cancer_type = "LUAD",
  description = "Li et al. Clin Respir J 2023 - LASSO Cox"
)

sig_xia2023 <- load_signature(
  genes       = c("TCN1", "CENPF", "MAOB", "CRTAC1", "PLEK2"),
  name        = "Xia_2023",
  cancer_type = "LUAD",
  description = "Xia et al. World J Clin Oncol 2023"
)

# --- Benchmark ---
comp_luad <- compare_signatures(
  signatures  = list(sig_chen2007, sig_li2023, sig_xia2023),
  expr_matrix = expr_final_luad,
  surv_time   = surv_t_luad,
  surv_event  = surv_e_luad
)
print(comp_luad)

# --- Staging karşılaştırması ---
scored_best_luad <- score_signature(
  sig_li2023, expr_final_luad
)

result_luad <- compare_to_staging(
  scored_best_luad, surv_t_luad, surv_e_luad, stage_luad_clean
)
print(result_luad)
