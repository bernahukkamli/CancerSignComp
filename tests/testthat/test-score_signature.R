set.seed(42)
mat <- matrix(
  rnorm(500), nrow = 5,
  dimnames = list(
    c("TP53", "KRAS", "MYC", "CDKN2A", "SMAD4"),
    paste0("Patient_", 1:100)
  )
)
sig <- load_signature(c("TP53", "KRAS", "MYC"), "TestSig", "PAAD")

test_that("score_signature returns ScoredSignature object", {
  scored <- score_signature(sig, mat)
  expect_s3_class(scored, "ScoredSignature")
})

test_that("score_signature returns correct number of patients", {
  scored <- score_signature(sig, mat)
  expect_equal(scored$n_patients, 100)
})

test_that("score_signature splits into equal groups with median cutoff", {
  scored <- score_signature(sig, mat, cutoff = "median")
  expect_equal(sum(scored$risk_group == "Low"),
               sum(scored$risk_group == "High"))
})

test_that("score_signature warns on missing genes", {
  sig2 <- load_signature(c("TP53", "BRCA1"), "TestSig2", "PAAD")
  expect_warning(score_signature(sig2, mat), "not found")
})

test_that("score_signature errors on non-matrix input", {
  expect_error(
    score_signature(sig, as.data.frame(mat)),
    "numeric matrix"
  )
})

test_that("score_signature errors on invalid method", {
  expect_error(
    score_signature(sig, mat, method = "invalid"),
    "must be one of"
  )
})

test_that("score_signature weighted_mean and ssgsea give different scores", {
  scored_wm <- score_signature(sig, mat, method = "weighted_mean")
  scored_ss <- score_signature(sig, mat, method = "ssgsea")
  expect_false(identical(scored_wm$risk_scores, scored_ss$risk_scores))
})
