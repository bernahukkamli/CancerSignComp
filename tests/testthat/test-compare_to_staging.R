set.seed(42)
mat <- matrix(
  rnorm(500), nrow = 5,
  dimnames = list(c("TP53","KRAS","MYC","CDKN2A","SMAD4"),
                  paste0("Patient_", 1:100))
)
sig    <- load_signature(c("TP53","KRAS","MYC"), "TestSig", "PAAD")
scored <- score_signature(sig, mat)
surv_t <- rexp(100, rate = 0.05)
surv_e <- rbinom(100, 1, 0.7)
stages <- sample(c("I","II","III","IV"), 100, replace = TRUE)

test_that("compare_to_staging returns StagingComparison object", {
  result <- compare_to_staging(scored, surv_t, surv_e, stages)
  expect_s3_class(result, "StagingComparison")
})

test_that("compare_to_staging has required fields", {
  result <- compare_to_staging(scored, surv_t, surv_e, stages)
  expect_true("sig_c_index"   %in% names(result))
  expect_true("stage_c_index" %in% names(result))
  expect_true("c_index_diff"  %in% names(result))
  expect_true("verdict"       %in% names(result))
})

test_that("compare_to_staging normalizes stage formats", {
  stages2 <- sample(c("Stage I","Stage II","stage iii","4"), 100, replace = TRUE)
  result  <- compare_to_staging(scored, surv_t, surv_e, stages2)
  expect_s3_class(result, "StagingComparison")
})

test_that("compare_to_staging warns on unknown stages", {
  stages3 <- c(stages[1:95], rep("Unknown", 5))
  expect_warning(
    compare_to_staging(scored, surv_t, surv_e, stages3),
    "unrecognized"
  )
})

test_that("compare_to_staging c_index values are between 0.5 and 1", {
  result <- compare_to_staging(scored, surv_t, surv_e, stages)
  expect_gte(result$sig_c_index,   0.5)
  expect_lte(result$sig_c_index,   1.0)
  expect_gte(result$stage_c_index, 0.5)
  expect_lte(result$stage_c_index, 1.0)
})

test_that("compare_to_staging verdict is one of three options", {
  result <- compare_to_staging(scored, surv_t, surv_e, stages)
  expect_true(result$verdict %in% c(
    "Signature significantly OUTPERFORMS clinical staging",
    "Signature significantly UNDERPERFORMS clinical staging",
    "Signature performs SIMILARLY to clinical staging"
  ))
})
