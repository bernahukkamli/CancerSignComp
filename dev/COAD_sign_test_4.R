# Zheng_2022 - 4-gen metabolik COAD imzası (Frontiers in Public Health 2022)
# DFS risk score formülü makaleden (sayfa 5)
sig_zheng2022 <- load_signature(
  genes  = c("CCND1", "EDAR", "FUT1", "PPAT"),
  name   = "Zheng_2022",
  cancer = "COAD",
  source = "Zheng et al. Front Public Health 2022"
)

weights_zheng2022 <- c(
  CCND1 =  0.411,
  EDAR  =  0.547,
  FUT1  =  0.302,
  PPAT  = -0.584
)

# Coverage kontrol
cat("Genes found:", sum(sig_zheng2022$genes %in% rownames(expr_sym_coad)),
    "/", sig_zheng2022$n_genes, "\n")

# Score et
scored_zheng <- score_signature(
  sig_zheng2022,
  expr_sym_coad,
  weights     = weights_zheng2022,
  cutoff      = "median",
  norm_method = "none"
)
cat("Coverage:", scored_zheng$coverage_pct, "%\n")

# Staging karşılaştırması
result_zheng <- compare_to_staging(
  scored     = scored_zheng,
  surv_time  = surv_time_coad,
  surv_event = surv_event_coad,
  stage      = stage_coad
)
print(result_zheng)
