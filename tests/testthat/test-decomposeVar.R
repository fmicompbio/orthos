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
    expect_error(.readGeneInformation("error", mustWork = FALSE))
    expect_error(.readGeneInformation("error", mustWork = TRUE))
    
    genesMouse <- .readGeneInformation("Mouse", mustWork = FALSE)
    genesHuman <- .readGeneInformation("Human", mustWork = FALSE)
    
    expect_s4_class(genesMouse, "DFrame")
    expect_s4_class(genesHuman, "DFrame")
    
    idTypes <- c("ENSEMBL_GENE_ID", "GENE_SYMBOL",
                 "ENTREZ_GENE_ID", "ARCHS4_ID")
    expect_true(all(idTypes %in% colnames(genesMouse)))
    expect_true(all(idTypes %in% colnames(genesHuman)))
    
    skip_if(nrow(genesMouse) > 0 || nrow(genesHuman) > 0,
            message = "cannot test `mustWork` when data is available")
    
    expect_error(.readGeneInformation("Mouse", mustWork = TRUE))
    expect_error(.readGeneInformation("Human", mustWork = TRUE))
})

## ------------------------------------------------------------------------- ##
## Checks, .detectFeatureIdType
## ------------------------------------------------------------------------- ##
test_that(".detectFeatureIdType works", {
    # load annotation and create synthetic data
    genesMouse <- .readGeneInformation("Mouse", mustWork = FALSE)
    genesHuman <- .readGeneInformation("Human", mustWork = FALSE)

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


## ------------------------------------------------------------------------- ##
## Checks, decomposeVar
## ------------------------------------------------------------------------- ##
test_that("decomposeVar works", {
    genesMouse <- .readGeneInformation("Mouse", mustWork = FALSE)
    genesHuman <- .readGeneInformation("Human", mustWork = FALSE)
    
    skip_if(nrow(genesMouse) == 0 || nrow(genesHuman) == 0,
            message = paste0("skipping decomposeVar tests - ",
                             "gene information not vailable"))

    fnameHuman <- system.file("extdata", "GSE215150_MKL1_Human.rds",
                              package = "orthos", mustWork = TRUE)
    fnameMouse <- system.file("extdata", "GSE215150_MKL1_Mouse.rds",
                              package = "orthos", mustWork = TRUE)

    countsHuman <- readRDS(fnameHuman)
    countsMouse <- readRDS(fnameMouse)
    
    decHuman <- decomposeVar(M = countsHuman, treatm = c(2, 3), cntr = c(1, 1), 
                             organism = "Human", verbose = FALSE)
    decMouse <- decomposeVar(M = countsMouse, treatm = c(2, 3), cntr = c(1, 1),
                             organism = "Mouse", verbose = FALSE)
    
    expect_s4_class(decHuman, "SummarizedExperiment")
    expect_s4_class(decMouse, "SummarizedExperiment")
    
    expect_identical(colnames(decHuman), c("MKL1", "caMKL1"))
    expect_identical(colnames(decMouse), c("MKL1", "caMKL1"))
    
    anms <- c("INPUT_CONTRASTS", "DECODED_CONTRASTS", "RESIDUAL_CONTRASTS", 
              "CONTEXT")
    expect_identical(SummarizedExperiment::assayNames(decHuman), anms)
    expect_identical(SummarizedExperiment::assayNames(decMouse), anms)
    
    for (anm1 in anms) {
        expect_type(assay(decHuman, anm1), "double")
        expect_type(assay(decMouse, anm1), "double")
    }
    
    expect_identical(dim(decHuman), c(20411L, 2L))
    expect_identical(dim(decMouse), c(20339L, 2L))
    
    expect_equal(colSums(assay(decHuman, "INPUT_CONTRASTS")),
                 c(MKL1 = -287.736062968947, caMKL1 = -2875.33388549955))
    expect_equal(colSums(assay(decHuman, "DECODED_CONTRASTS")),
                 c(MKL1 = -345.53792401588, caMKL1 = -2228.12848719998))
    expect_equal(colSums(assay(decHuman, "RESIDUAL_CONTRASTS")),
                 c(MKL1 = 57.8018610469331, caMKL1 = -647.205398299573))
    expect_equal(colSums(assay(decHuman, "CONTEXT")),
                 c(MKL1 = 78030.7107896163, caMKL1 = 78030.7107896163))
    
    expect_equal(colSums(assay(decMouse, "INPUT_CONTRASTS")),
                 c(MKL1 = -205.650429365804, caMKL1 = -637.007214334137))
    expect_equal(colSums(assay(decMouse, "DECODED_CONTRASTS")),
                 c(MKL1 = -99.0546978751445, caMKL1 = -586.600493863487))
    expect_equal(colSums(assay(decMouse, "RESIDUAL_CONTRASTS")),
                 c(MKL1 = -106.595731490659, caMKL1 = -50.40672047065))
    expect_equal(colSums(assay(decMouse, "CONTEXT")),
                 c(MKL1 = 79313.4220740248, caMKL1 = 79313.4220740248))
    
    idsHuman <- intersect(rownames(countsHuman), rownames(decHuman))
    idsMouse <- intersect(rownames(countsMouse), rownames(decMouse))
    
    expect_length(idsHuman, 18051L)
    expect_length(idsMouse, 19776L)
    
    decHuman2 <- decomposeVar(M = countsHuman[idsHuman, 2:3],
                              MD = assay(decHuman, "INPUT_CONTRASTS")[idsHuman, ], 
                              organism = "Human", verbose = TRUE)
    decMouse2 <- decomposeVar(M = countsMouse[idsMouse, 2:3],
                              MD = assay(decMouse, "INPUT_CONTRASTS")[idsMouse, ],
                              organism = "Mouse", verbose = TRUE)
    
    expect_identical(dim(decHuman), dim(decHuman2))
    expect_identical(dim(decMouse), dim(decMouse2))
    
    expect_true(all(
        diag(cor(assay(decHuman, "RESIDUAL_CONTRASTS"),
                 assay(decHuman2, "RESIDUAL_CONTRASTS"))) > 0.99
    ))
    expect_true(all(
        diag(cor(assay(decMouse, "RESIDUAL_CONTRASTS"),
                 assay(decMouse2, "RESIDUAL_CONTRASTS"))) > 0.64
    ))
})
