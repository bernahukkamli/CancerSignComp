set.seed(42)
mat <- matrix(
  rnorm(1000), nrow = 10,
  dimnames = list(
    c("TP53","KRAS","MYC","CDKN2A","SMAD4",
      "EGFR","PTEN","RB1","BRCA2","PIK3CA"),
    paste0("Patient_", 1:100)
  )
)
sig1   <- load_signature(c("TP53","KRAS","MYC"),  "Sig_A", "PAAD")
sig2   <- load_signature(c("EGFR","PTEN","RB1"),  "Sig_B", "PAAD")
surv_t <- rexp(100, rate = 0.05)
surv_e <- rbinom(100, 1, 0.7)

test_that("compare_signatures returns SignatureComparison object", {
  comp <- compare_signatures(list(sig1, sig2), mat, surv_t, surv_e)
  expect_s3_class(comp, "SignatureComparison")
})

test_that("compare_signatures summary_table has correct dimensions", {
  comp <- compare_signatures(list(sig1, sig2), mat, surv_t, surv_e)
  expect_equal(nrow(comp$summary_table), 2)
  expect_true("c_index" %in% names(comp$summary_table))
  expect_true("logrank_p" %in% names(comp$summary_table))
})

test_that("compare_signatures ranks by c_index descending", {
  comp <- compare_signatures(list(sig1, sig2), mat, surv_t, surv_e)
  ci <- comp$summary_table$c_index
  expect_true(all(diff(ci) <= 0))
})

test_that("compare_signatures errors on fewer than 2 signatures", {
  expect_error(
    compare_signatures(list(sig1), mat, surv_t, surv_e),
    "at least 2"
  )
})

test_that("compare_signatures errors on wrong surv_time length", {
  expect_error(
    compare_signatures(list(sig1, sig2), mat, surv_t[1:50], surv_e),
    "number of patients"
  )
})

test_that("compare_signatures c_index values are between 0.5 and 1", {
  comp <- compare_signatures(list(sig1, sig2), mat, surv_t, surv_e)
  expect_true(all(comp$summary_table$c_index >= 0.5))
  expect_true(all(comp$summary_table$c_index <= 1))
})
