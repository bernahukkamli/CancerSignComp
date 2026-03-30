sig_aldherasi2021 <- load_signature(
  genes  = c("GRIA1", "UCN2", "PKHD1L1", "RIMS2", "PGM5", "CLIC6", "CAVIN2"),
  name   = "AlDherasi_2021",
  cancer = "LUAD",
  source = "Al-Dherasi et al. Cancer Cell Int 2021"
)

weights_aldherasi2021 <- c(
  GRIA1   = -0.3658,
  UCN2    =  0.5701,
  PKHD1L1 = -0.6010,
  RIMS2   =  0.2192,
  PGM5    = -0.3617,
  CLIC6   = -0.6036,
  CAVIN2  =  1.1686
)

cat("Genes found:", sum(sig_aldherasi2021$genes %in% rownames(expr_sym_luad)),
    "/", sig_aldherasi2021$n_genes, "\n")

scored_aldherasi <- score_signature(
  sig_aldherasi2021,
  expr_sym_luad,
  weights     = weights_aldherasi2021,
  cutoff      = "median",
  norm_method = "none"
)
cat("Coverage:", scored_aldherasi$coverage_pct, "%\n")

result_aldherasi <- compare_to_staging(
  scored     = scored_aldherasi,
  surv_time  = surv_time_luad,
  surv_event = surv_event_luad,
  stage      = stage_named
)
print(result_aldherasi)
