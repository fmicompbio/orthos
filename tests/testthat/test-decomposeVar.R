## ------------------------------------------------------------------------- ##
## Checks, .preprocessInput
## ------------------------------------------------------------------------- ##
test_that(".preprocessInput works", {
    # create synthetic data
    set.seed(1)
    M <- matrix(stats::rpois(30, 20), nrow = 10)
    rownames(M) <- as.character(seq.int(nrow(M)))
    M2 <- M
    M2[M < 20] <- NA
    
    # checks without NAs
    N0 <- log2(sweep(M, 2, colSums(M), "/") * 1e6 + 4)
    expect_message({
        N1 <- .preprocessInput(M = M, ids = rownames(M),
                               verbose = TRUE, pseudocount = 4)
    })
    N2 <- .preprocessInput(M = M, ids = as.character(1:5),
                           verbose = FALSE, pseudocount = 4)
    
    expect_type(N1, "double")
    expect_type(N2, "double")
    expect_identical(dim(M), dim(N1))
    expect_identical(c(5L, ncol(M)), dim(N2))
    expect_equal(N0, N1)
    expect_equal(N0[rownames(N2), ], N2)
    
    # checks with NAs
    N0 <- log2(sweep(M2, 2, colSums(M2, na.rm = TRUE), "/") * 1e6 + 4)
    N0[is.na(N0)] <- log2(4)
    expect_message(
        expect_message({
            N1 <- .preprocessInput(M = M2, ids = rownames(M2),
                                   verbose = TRUE, pseudocount = 4)
        })
    )
    N2 <- .preprocessInput(M = M2, ids = as.character(1:5),
                           verbose = FALSE, pseudocount = 4)
    
    expect_type(N1, "double")
    expect_type(N2, "double")
    expect_identical(dim(M2), dim(N1))
    expect_identical(c(5L, ncol(M2)), dim(N2))
    expect_equal(N0, N1)
    expect_equal(N0[rownames(N2), ], N2)
})
