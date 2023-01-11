## ------------------------------------------------------------------------- ##
## Checks, testdeJUNKERenv
## ------------------------------------------------------------------------- ##
test_that("testdeJUNKERenv works", {
    res <- testdeJUNKERenv()
    
    expect_identical(res, list(keras_available = TRUE, tf_version = "2.10.0"))
})
