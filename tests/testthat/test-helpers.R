## ------------------------------------------------------------------------- ##
## Checks, M2Mcor
## ------------------------------------------------------------------------- ##
test_that(".M2Mcor works", {
    set.seed(1)
    M1 <- matrix(stats::rnorm(30), nrow = 10)
    M2 <- matrix(stats::rnorm(30), nrow = 10)
    expect_equal(.M2Mcor(M1, M2),
                 c(stats::cor(M1[, 1], M2[, 1]),
                   stats::cor(M1[, 2], M2[, 2]),
                   stats::cor(M1[, 3], M2[, 3])))
})

## ------------------------------------------------------------------------- ##
## Checks, .grid_cor_wNAs
## ------------------------------------------------------------------------- ##
test_that(".grid_cor_wNAs works", {
    # create synthetic data
    nr <- 100L
    ncQuery <- 20L
    ncDBase <- 50L
    thr <- 3.0
    set.seed(56L)
    mQuery <- mQueryThr <- matrix(runif(nr * ncQuery) * 12, nrow = nr, ncol = ncQuery)
    mQueryThr[mQuery < thr] <- NA
    mDBase <- mDBaseThr <- matrix(runif(nr * ncDBase) * 12, nrow = nr, ncol = ncDBase)
    mDBaseThr[mDBase < thr] <- NA
    tf <- tempfile(fileext = ".h5")
    mDBaseHDF5 <- HDF5Array::writeHDF5Array(x = mDBase, filepath = tf)
    tfthr <- tempfile(fileext = ".h5")
    mDBaseThrHDF5 <- HDF5Array::writeHDF5Array(x = mDBaseThr, filepath = tfthr)

    # misspecified arguments
    expect_error(.grid_cor_wNAs(query = NULL, hdf5 = mDBaseHDF5), "NULL")
    expect_error(.grid_cor_wNAs(query = "error", hdf5 = mDBaseHDF5), "matrix")
    expect_error(.grid_cor_wNAs(query = mQuery, hdf5 = NULL), "NULL")
    expect_error(.grid_cor_wNAs(query = mQuery, hdf5 = "error"), "HDF5Matrix")
    expect_error(.grid_cor_wNAs(query = mQuery, hdf5 = mDBaseHDF5, chunk_size = "error"), "numeric")
    expect_error(.grid_cor_wNAs(query = mQuery, hdf5 = mDBaseHDF5, workers = "error"), "numeric")
    expect_error(.grid_cor_wNAs(query = mQuery, hdf5 = mDBaseHDF5, workers = -1), "within")
    expect_error(.grid_cor_wNAs(query = mQuery, hdf5 = mDBaseHDF5, thr = "error"), "numeric")
    
    # correct results
    res0 <- stats::cor(mQuery, mDBase)
    res0thr <- stats::cor(mQueryThr, mDBaseThr, use = "pairwise.complete")
    res1 <- .grid_cor_wNAs(query = mQuery, hdf5 = mDBaseHDF5, chunk_size = 10, workers = 1, thr = 0.0)
    res1thr <- .grid_cor_wNAs(query = mQueryThr, hdf5 = mDBaseThrHDF5, chunk_size = 10, workers = 1, thr = thr)
    res2 <- .grid_cor_wNAs(query = mQuery, hdf5 = mDBaseHDF5, chunk_size = 10, workers = 2, thr = 0.0)
    res2thr <- .grid_cor_wNAs(query = mQueryThr, hdf5 = mDBaseThrHDF5, chunk_size = 10, workers = 2, thr = thr)
    
    expect_type(res1, "double")
    expect_type(res1thr, "double")
    expect_type(res2, "double")
    expect_type(res2thr, "double")
    
    expect_identical(dim(res1), c(ncQuery, ncDBase))
    expect_identical(dim(res1thr), c(ncQuery, ncDBase))
    expect_identical(dim(res2), c(ncQuery, ncDBase))
    expect_identical(dim(res2thr), c(ncQuery, ncDBase))
    
    expect_equal(res0, res1)
    expect_equal(res0thr, res1thr)
    expect_equal(res0, res2)
    expect_equal(res0thr, res2thr)
    
    # clean up
    unlink(c(tf, tfthr))
})


## ------------------------------------------------------------------------- ##
## Checks, .grid_cor_woNAs
## ------------------------------------------------------------------------- ##
test_that(".grid_cor_woNAs works", {
    # create synthetic data
    nr <- 100L
    ncQuery <- 20L
    ncDBase <- 50L
    set.seed(57L)
    mQuery <- mQueryThr <- matrix(runif(nr * ncQuery) * 12, nrow = nr, ncol = ncQuery)
    mDBase <- mDBaseThr <- matrix(runif(nr * ncDBase) * 12, nrow = nr, ncol = ncDBase)
    tf <- tempfile(fileext = ".h5")
    mDBaseHDF5 <- HDF5Array::writeHDF5Array(x = mDBase, filepath = tf)

    # misspecified arguments
    expect_error(.grid_cor_woNAs(query = NULL, hdf5 = mDBaseHDF5), "NULL")
    expect_error(.grid_cor_woNAs(query = "error", hdf5 = mDBaseHDF5), "matrix")
    expect_error(.grid_cor_woNAs(query = mQuery, hdf5 = NULL), "NULL")
    expect_error(.grid_cor_woNAs(query = mQuery, hdf5 = "error"), "HDF5Matrix")
    expect_error(.grid_cor_woNAs(query = mQuery, hdf5 = mDBaseHDF5, chunk_size = "error"), "numeric")
    expect_error(.grid_cor_woNAs(query = mQuery, hdf5 = mDBaseHDF5, workers = "error"), "numeric")
    expect_error(.grid_cor_woNAs(query = mQuery, hdf5 = mDBaseHDF5, workers = -1), "within")

    # correct results
    res0 <- stats::cor(mQuery, mDBase)
    res1 <- .grid_cor_woNAs(query = mQuery, hdf5 = mDBaseHDF5, chunk_size = 10, workers = 1)
    res2 <- .grid_cor_woNAs(query = mQuery, hdf5 = mDBaseHDF5, chunk_size = 10, workers = 2)
    res3 <- .grid_cor_wNAs(query = mQuery, hdf5 = mDBaseHDF5, chunk_size = 10, workers = 1, thr = 0.0)
    
    expect_type(res1, "double")
    expect_type(res2, "double")

    expect_identical(dim(res1), c(ncQuery, ncDBase))
    expect_identical(dim(res2), c(ncQuery, ncDBase))

    expect_equal(res0, res1)
    expect_equal(res0, res2)
    expect_equal(res1, res3)
    
    # clean up
    unlink(tf)
})

## ------------------------------------------------------------------------- ##
## Checks, .assertScalar
## ------------------------------------------------------------------------- ##
test_that(".assertScalar works", {
    expect_error(.assertScalar(1, type = TRUE))
    expect_error(.assertScalar(1, type = 1))
    expect_error(.assertScalar(1, type = c("numeric", "character")))
    expect_error(.assertScalar(1, type = "numeric", rngIncl = TRUE))
    expect_error(.assertScalar(1, type = "numeric", rngIncl = "rng"))
    expect_error(.assertScalar(1, type = "numeric", rngIncl = 1))
    expect_error(.assertScalar(1, type = "numeric", rngIncl = 1:3))
    expect_error(.assertScalar(1, type = "numeric", rngExcl = TRUE))
    expect_error(.assertScalar(1, type = "numeric", rngExcl = "rng"))
    expect_error(.assertScalar(1, type = "numeric", rngExcl = 1))
    expect_error(.assertScalar(1, type = "numeric", rngExcl = 1:3))
    expect_error(.assertScalar(1, type = "numeric", rngIncl = c(0, 2), rngExcl = c(0, 2)))
    expect_error(.assertScalar(1, type = "numeric", allowNULL = 1))
    expect_error(.assertScalar(1, type = "numeric", allowNULL = "rng"))
    expect_error(.assertScalar(1, type = "numeric", allowNULL = NULL))
    expect_error(.assertScalar(1, type = "numeric", allowNULL = c(TRUE, FALSE)))
    
    expect_true(.assertScalar(1, type = "numeric", rngIncl = c(1, 3)))
    expect_error(.assertScalar(1, type = "numeric", rngExcl = c(1, 3)))
    expect_true(.assertScalar(1, type = "numeric", rngExcl = c(1, 3), validValues = 1))
    expect_true(.assertScalar(-1, type = "numeric", rngIncl = c(1, 3), validValues = c(-1, 0)))
    expect_error(.assertScalar(-1, type = "numeric", rngIncl = c(1, 3), validValues = 0))
    expect_true(.assertScalar(-1, type = "numeric", validValues = c(-1, 0)))
    expect_error(.assertScalar(-1, type = "numeric", validValues = c(-2, 0)))
    expect_true(.assertScalar(NA_real_, type = "numeric", rngIncl = c(1, 2), validValues = NA_real_))
    expect_error(.assertScalar(NA, type = "numeric", rngIncl = c(1, 2), validValues = NA_real_))
    expect_true(.assertScalar(NA_real_, type = "numeric", rngIncl = c(1, 2), validValues = NA))
    expect_true(.assertScalar(1, type = "numeric", rngIncl = c(0, 3), validValues = 3))
    expect_true(.assertScalar(1, rngIncl = c(0, 3), validValues = 3))
    expect_true(.assertScalar(1, type = "numeric", rngIncl = c(0, 1)))
    expect_error(.assertScalar(1, type = "numeric", rngExcl = c(0, 1)))
    expect_true(.assertScalar(1, type = "numeric", rngExcl = c(0, 1), validValues = 1))
    expect_error(.assertScalar(1, type = "numeric", rngExcl = c(0, 1), validValues = 3:4))
    expect_true(.assertScalar(NULL, type = "numeric", allowNULL = TRUE))
    expect_error(.assertScalar(NULL, type = "numeric", allowNULL = FALSE))
    expect_error(.assertScalar(1, type = "character"))
    expect_error(.assertScalar("x", type = "numeric"))
    expect_error(.assertScalar(FALSE, type = "character"))
    expect_error(.assertScalar(c(1, 2), type = "numeric"))
    test <- "text"
    expect_error(.assertScalar(x = test, type = "numeric"),
                 "'test' must be of class 'numeric")
})

## ------------------------------------------------------------------------- ##
## Checks, .assertVector
## ------------------------------------------------------------------------- ##
test_that(".assertVector works", {
    expect_error(.assertVector(1, type = TRUE))
    expect_error(.assertVector(1, type = 1))
    expect_error(.assertVector(1, type = c("numeric", "character")))
    expect_error(.assertVector(1, type = "numeric", rngIncl = TRUE))
    expect_error(.assertVector(1, type = "numeric", rngIncl = "rng"))
    expect_error(.assertVector(1, type = "numeric", rngIncl = 1))
    expect_error(.assertVector(1, type = "numeric", rngIncl = 1:3))
    expect_error(.assertVector(1, type = "numeric", rngExcl = TRUE))
    expect_error(.assertVector(1, type = "numeric", rngExcl = "rng"))
    expect_error(.assertVector(1, type = "numeric", rngExcl = 1))
    expect_error(.assertVector(1, type = "numeric", rngExcl = 1:3))
    expect_error(.assertVector(1, type = "numeric", rngIncl = c(0, 2), rngExcl = c(0, 2)))
    expect_error(.assertVector(1, type = "numeric", allowNULL = 1))
    expect_error(.assertVector(1, type = "numeric", allowNULL = "rng"))
    expect_error(.assertVector(1, type = "numeric", allowNULL = NULL))
    expect_error(.assertVector(1, type = "numeric", allowNULL = c(TRUE, FALSE)))
    expect_error(.assertVector(1, type = "numeric", len = TRUE))
    expect_error(.assertVector(1, type = "numeric", len = "rng"))
    expect_error(.assertVector(1, type = "numeric", len = 1:3))
    expect_error(.assertVector(1, type = "numeric", rngLen = TRUE))
    expect_error(.assertVector(1, type = "numeric", rngLen = "rng"))
    expect_error(.assertVector(1, type = "numeric", rngLen = 1))
    expect_error(.assertVector(1, type = "numeric", rngLen = 1:3))
    
    expect_true(.assertVector(c(1, 2), type = "numeric", rngIncl = c(1, 3)))
    expect_error(.assertVector(c(1, 2), type = "numeric", rngIncl = c(1, 1.5)))
    expect_error(.assertVector(c(1, 2), type = "numeric", rngExcl = c(1, 3)))
    expect_true(.assertVector(c(1, 2), type = "numeric", rngExcl = c(1, 3), validValues = 1))
    expect_error(.assertVector(c(1, 2), type = "numeric", validValues = c(1, 3)))
    expect_true(.assertVector(c(1, 2), type = "numeric", validValues = c(1, 2)))
    expect_error(.assertVector(c(1, 2), type = "numeric", len = 1))
    expect_true(.assertVector(c(1, 2), type = "numeric", len = 2))
    expect_error(.assertVector(c(1, 2), type = "numeric", rngLen = c(3, 5)))
    expect_true(.assertVector(c(1, 2), type = "numeric", rngLen = c(2, 5)))
    expect_true(.assertVector(c(1, 2), type = "numeric", rngLen = c(1, 2)))
    expect_error(.assertVector(c("a", "b"), type = "character", validValues = c("A", "B")))
    expect_true(.assertVector(LETTERS[1:2], type = "character", validValues = LETTERS))
    test <- "text"
    expect_error(.assertVector(x = test, type = "numeric"),
                 "'test' must be of class 'numeric")
})

## ------------------------------------------------------------------------- ##
## Checks, .assertPackagesAvailable
## ------------------------------------------------------------------------- ##
test_that(".assertPackagesAvailable works", {
    testfunc <- function(...) .assertPackagesAvailable(...)
    expect_error(testfunc(1L))
    expect_error(testfunc("test", "error"))
    expect_error(testfunc("test", c(TRUE, FALSE)))

    expect_true(testfunc("base"))
    expect_true(testfunc("githubuser/base"))
    expect_true(testfunc(c("base", "methods")))
    expect_error(testfunc(c("error", "error2")), "BiocManager")
    expect_error(testfunc("error1", suggestInstallation = FALSE), "installed.\n$")
    rm(testfunc)
})
