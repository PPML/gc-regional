test_that("the order of the names in theta matches the priors data frame.", {
  expect_true(all(sapply(strsplit(names(gc_env$theta), "\\."), function(x)
    paste0(x[2:length(x)], collapse = ".")) == gc_env$priors$parameter))
})
