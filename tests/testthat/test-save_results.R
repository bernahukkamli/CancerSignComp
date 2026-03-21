set.seed(42)
mat <- matrix(
  rnorm(1000), nrow = 10,
  dimnames = list(
    c("TP53","KRAS","MYC","CDKN2A","SMAD4",
      "EGFR","PTEN","RB1","BRCA2","PIK3CA"),
    paste0("Patient_", 1:100)
  )
)
sig1   <- load_signature(c("TP53","KRAS","MYC"), "Sig_A", "PAAD")
sig2   <- load_signature(c("EGFR","PTEN","RB1"), "Sig_B", "PAAD")
surv_t <- rexp(100, rate = 0.05)
surv_e <- rbinom(100, 1, 0.7)
comp   <- compare_signatures(list(sig1, sig2), mat, surv_t, surv_e)

test_that("save_results creates CSV file", {
  tmp <- tempfile(fileext = ".csv")
  save_results(comp, tmp)
  expect_true(file.exists(tmp))
  df <- read.csv(tmp)
  expect_equal(nrow(df), 2)
  unlink(tmp)
})

test_that("save_results errors on existing file without overwrite", {
  tmp <- tempfile(fileext = ".csv")
  save_results(comp, tmp)
  expect_error(save_results(comp, tmp), "already exists")
  unlink(tmp)
})

test_that("save_results overwrites with overwrite = TRUE", {
  tmp <- tempfile(fileext = ".csv")
  save_results(comp, tmp)
  expect_no_error(save_results(comp, tmp, overwrite = TRUE))
  unlink(tmp)
})

test_that("save_results errors on wrong extension", {
  expect_error(save_results(comp, "file.txt"), ".csv or .xlsx")
})
