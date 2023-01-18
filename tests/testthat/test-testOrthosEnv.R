## ------------------------------------------------------------------------- ##
## Checks, testOrthosEnv
## ------------------------------------------------------------------------- ##
test_that("testOrthosEnv works", {
    res <- testOrthosEnv()
    
    expect_identical(res, list(keras_available = TRUE, tf_version = "2.10.0"))
})
