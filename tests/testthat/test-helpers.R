test_that("M2Mcor works", {
    set.seed(1)
    M1 <- matrix(stats::rnorm(30), nrow = 10)
    M2 <- matrix(stats::rnorm(30), nrow = 10)
    expect_equal(M2Mcor(M1, M2),
                 c(stats::cor(M1[, 1], M2[, 1]),
                   stats::cor(M1[, 2], M2[, 2]),
                   stats::cor(M1[, 3], M2[, 3])))
})
