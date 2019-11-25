
test_that("the gc_env environment is created when the gcRegional package is loaded", {
    expect_true(exists('gc_env', as.environment('package:gcRegional'))) })


