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

## ------------------------------------------------------------------------- ##
## Checks, .readGeneInformation
## ------------------------------------------------------------------------- ##
test_that(".readGeneInformation works", {
    expect_error(.readGeneInformation("error", mustSucceed = FALSE))
    expect_error(.readGeneInformation("error", mustSucceed = TRUE))
    
    genesMouse <- .readGeneInformation("mouse", mustSucceed = FALSE)
    genesHuman <- .readGeneInformation("human", mustSucceed = FALSE)
    
    if (nrow(genesMouse) == 0 && nrow(genesHuman) == 0) {
        expect_error(.readGeneInformation("mouse", mustSucceed = TRUE))
        expect_error(.readGeneInformation("human", mustSucceed = TRUE))
        
    } else {
        idTypes <- c("ENSEMBL_GENE_ID", "GENE_SYMBOL",
                     "ENTREZ_GENE_ID", "ARCHS4_ID")
        
        expect_true(all(idTypes %in% colnames(genesMouse)))
        expect_true(all(idTypes %in% colnames(genesHuman)))
    }
})

## ------------------------------------------------------------------------- ##
## Checks, .detectFeatureIdType
## ------------------------------------------------------------------------- ##
test_that(".detectFeatureIdType works", {
    # load annotation and create synthetic data
    genesMouse <- .readGeneInformation("mouse", mustSucceed = FALSE)
    genesHuman <- .readGeneInformation("human", mustSucceed = FALSE)
    skip_if(nrow(genesMouse) == 0 || nrow(genesHuman) == 0,
            message = paste0("skipping .detectFeatureIdType tests - ",
                             "gene information not vailable"))
    
    idTypes <- c("ENSEMBL_GENE_ID", "GENE_SYMBOL",
                 "ENTREZ_GENE_ID", "ARCHS4_ID")
    nr <- 1000
    set.seed(926L)
    mMouse <- lapply(structure(idTypes, names = idTypes), function(idType) {
        ids <- sample(genesMouse[,idType], nr)
        ids[is.na(ids)] <- paste0("id_", seq.int(sum(is.na(ids))))
        matrix(runif(2 * nr, 0, 100), nrow = nr, ncol = 2,
               dimnames = list(ids, 1:2))
    })
    mHuman <- lapply(structure(idTypes, names = idTypes), function(idType) {
        ids <- sample(genesHuman[,idType], nr)
        ids[is.na(ids)] <- paste0("id_", seq.int(sum(is.na(ids))))
        matrix(runif(2 * nr, 0, 100), nrow = nr, ncol = 2,
               dimnames = list(ids, 1:2))
    })

    # verbose
    expect_message(
        expect_message(
            expect_message(.detectFeatureIdType(featureType = "AUTO",
                                                genes = genesMouse,
                                                M = mMouse[["GENE_SYMBOL"]],
                                                maxMissing = nrow(genesMouse),
                                                verbose = TRUE)
            )
        )
    )
    
    # too few matched genes
    expect_error(.detectFeatureIdType(featureType = "AUTO",
                                      genes = genesMouse,
                                      M = mMouse[["GENE_SYMBOL"]],
                                      maxMissing = 1000,
                                      verbose = FALSE)
    )

    # result checks
    resMouseL <- lapply(structure(c("AUTO",idTypes), names = c("AUTO",idTypes)),
                        function(idType) {
                            .detectFeatureIdType(featureType = idType,
                                                 genes = genesMouse,
                                                 M = mMouse[[idType]],
                                                 maxMissing = nrow(genesMouse),
                                                 verbose = FALSE)
                   })
    resHumanL <- lapply(structure(c("AUTO",idTypes), names = c("AUTO",idTypes)),
                        function(idType) {
                            .detectFeatureIdType(featureType = idType,
                                                 genes = genesHuman,
                                                 M = mHuman[[idType]],
                                                 maxMissing = nrow(genesHuman),
                                                 verbose = FALSE)
                        })
    resExpected <- list(AUTO = "ENSEMBL_GENE_ID", ENSEMBL_GENE_ID = "ENSEMBL_GENE_ID", 
                        GENE_SYMBOL = "GENE_SYMBOL", ENTREZ_GENE_ID = "ENTREZ_GENE_ID", 
                        ARCHS4_ID = "ARCHS4_ID")
    expect_identical(resMouseL, resExpected)
    expect_identical(resHumanL, resExpected)
})
