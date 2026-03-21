test_that("load_signature returns CancerSignature object", {
  sig <- load_signature(c("TP53", "KRAS", "MYC"), "TestSig", "PAAD")
  expect_s3_class(sig, "CancerSignature")
  expect_equal(sig$name, "TestSig")
  expect_equal(sig$cancer_type, "PAAD")
  expect_equal(sig$n_genes, 3)
})

test_that("load_signature converts to uppercase", {
  sig <- load_signature(c("tp53", "kras"), "TestSig", "PAAD")
  expect_equal(sig$genes, c("TP53", "KRAS"))
})

test_that("load_signature removes duplicates", {
  sig <- load_signature(c("TP53", "TP53", "KRAS"), "TestSig", "PAAD")
  expect_equal(sig$n_genes, 2)
})

test_that("load_signature errors on single gene", {
  expect_error(
    load_signature("TP53", "TestSig", "PAAD"),
    "at least 2"
  )
})

test_that("load_signature errors on non-character input", {
  expect_error(
    load_signature(c(1, 2, 3), "TestSig", "PAAD"),
    "character vector"
  )
})

test_that("load_signature warns on invalid gene symbols", {
  expect_warning(
    load_signature(c("TP53", "KRAS", "123INVALID!!!"), "TestSig", "PAAD"),
    "invalid"
  )
})

test_that("load_signature errors on empty name", {
  expect_error(
    load_signature(c("TP53", "KRAS"), "", "PAAD"),
    "non-empty"
  )
})
