sig_zhang2022 <- load_signature(
  name   = "Zhang_2022",
  genes  = c("SCAF11", "NOD1", "NLRP2", "NLRP1", "GPX4", "CASP8", "AIM2"),
  cancer = "LUAD"
)

weights_zhang2022 <- c(
  SCAF11 = -0.280,
  NOD1   = -0.261,
  NLRP2  = -0.076,
  NLRP1  = -0.141,
  GPX4   = -0.301,
  CASP8  =  0.293,
  AIM2   =  0.131
)

# Gene coverage kontrol et
sum(sig_zhang2022$genes %in% rownames(expr_sym_luad))

# Score et
scored_zhang <- score_signature(
  sig_zhang2022,
  expr_sym_luad,
  weights = weights_zhang2022,
  cutoff  = "median",
  norm_method = "none"
)

cat("Genes used:", paste(scored_zhang$genes_used, collapse=", "), "\n")
cat("Coverage:", scored_zhang$coverage_pct, "%\n")
