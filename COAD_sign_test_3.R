# Cai_2021 - 11-gen COAD imzası (Gastroenterology Report 2021)
# Katsayılar Table 2'den (multivariable Cox)
sig_cai2021 <- load_signature(
  genes  = c("NDRG1", "FLT1", "LBP", "FABP4", "ADIPOQ",
             "AGT", "ACVRL1", "CCL11", "CDC42", "TRAV9_2", "POMC"),
  name   = "Cai_2021",
  cancer = "COAD",
  source = "Cai et al. Gastroenterology Report 2021"
)

weights_cai2021 <- c(
  POMC    = -5.752,
  AGT     = -1.497,
  LBP     =  4.162,
  CCL11   = -0.929,
  ACVRL1  = -2.563,
  TRAV9_2 = -3.124,
  CDC42   = -2.958,
  NDRG1   =  0.828,
  FABP4   =  0.710,
  ADIPOQ  = -0.356,
  FLT1    =  0.118
)

# Coverage kontrol
cat("Genes found:", sum(sig_cai2021$genes %in% rownames(expr_sym_coad)),
    "/", sig_cai2021$n_genes, "\n")

# Score et
scored_cai <- score_signature(
  sig_cai2021,
  expr_sym_coad,
  weights     = weights_cai2021,
  cutoff      = "median",
  norm_method = "none"
)
cat("Coverage:", scored_cai$coverage_pct, "%\n")
