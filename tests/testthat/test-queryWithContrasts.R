## ------------------------------------------------------------------------- ##
## Checks, loadContrastDatabase
## ------------------------------------------------------------------------- ##
test_that("loadContrastDatabase works", {
    expect_error(loadContrastDatabase(organism = "error"))
    expect_error(loadContrastDatabase(organism = "Human", mode = "error"))
    
    useMode <- "DEMO"
    seHuman <- loadContrastDatabase("Human", mode = useMode, mustWork = FALSE)
    seMouse <- loadContrastDatabase("Mouse", mode = useMode, mustWork = FALSE)
    
    expect_s4_class(seHuman, "SummarizedExperiment")
    expect_s4_class(seMouse, "SummarizedExperiment")
    
    skip_if(identical(dim(seHuman), c(0L, 0L)) ||
                identical(dim(seMouse), c(0L, 0L)),
            message = "contrast database not available")
    
    assayNms <- c("INPUT_CONTRASTS", "DECODED_CONTRASTS",
                  "RESIDUAL_CONTRASTS", "CONTEXT")
    expect_identical(SummarizedExperiment::assayNames(seHuman),
                     assayNms)
    expect_identical(SummarizedExperiment::assayNames(seMouse),
                     assayNms)

    if (identical(useMode, "ANALYSIS")) {
        # mode = "ANALYSIS"
        expect_identical(dim(seHuman), c(20411L, 74731L))
        expect_identical(dim(seMouse), c(20339L, 58532L))
    } else {
        # mode = "DEMO"
        expect_identical(dim(seHuman), c(20411L, 1000L))
        expect_identical(dim(seMouse), c(20339L,  988L))
    }

    for (assayNm1 in assayNms) {
        expect_s4_class(SummarizedExperiment::assay(seHuman, assayNm1),
                        "DelayedArray")
        expect_s4_class(SummarizedExperiment::assay(seMouse, assayNm1),
                        "DelayedArray")
    }
})

## ------------------------------------------------------------------------- ##
## Checks, queryWithContrasts
## ------------------------------------------------------------------------- ##
test_that("queryWithContrasts works", {
    genesHuman <- .readGeneInformation("Human", mustWork = FALSE)
    
    skip_if(nrow(genesHuman) == 0,
            message = paste0("skipping queryWithContrasts tests - ",
                             "gene information not vailable"))
    
    fnameHuman <- system.file("extdata", "GSE215150_MKL1_Human.rds",
                              package = "orthos", mustWork = TRUE)
    countsHuman <- readRDS(fnameHuman)

    decHuman <- decomposeVar(M = countsHuman, treatm = c(2, 3), cntr = c(1, 1), 
                             organism = "Human", verbose = FALSE)
    
    resHuman <- queryWithContrasts(decHuman, organism = "Human",
                                   detailTopn = 10L,
                                   BPPARAM = BiocParallel::SerialParam(),
                                   verbose = FALSE, plotType = "none",
                                   mode = "DEMO")
    resHuman2 <- queryWithContrasts(decHuman, organism = "Human",
                                    detailTopn = 10L,
                                    BPPARAM = BiocParallel::MulticoreParam(2L),
                                    verbose = TRUE, plotType = "none",
                                    mode = "DEMO")
    
    anms <- c("INPUT_CONTRASTS", "DECODED_CONTRASTS", "RESIDUAL_CONTRASTS")
    cnms <- colnames(decHuman)
    
    # serial vs parallel
    expect_identical(resHuman, resHuman2)
    
    # result type
    expect_type(resHuman, "list")
    expect_length(resHuman, 3L)
    expect_identical(names(resHuman), c("pearson.rhos", "zscores", "TopHits"))
    
    # pearson.rhos
    expect_type(resHuman$pearson.rhos, "list")
    expect_identical(names(resHuman$pearson.rhos), anms)
    expect_identical(dim(resHuman$pearson.rhos$INPUT_CONTRASTS),
                     dim(resHuman$pearson.rhos$DECODED_CONTRASTS))
    expect_identical(dim(resHuman$pearson.rhos$INPUT_CONTRASTS),
                     dim(resHuman$pearson.rhos$RESIDUAL_CONTRASTS))
    expect_equal(rowSums(resHuman$pearson.rhos$INPUT_CONTRASTS),
                 c(MKL1 = 11.8424116074152, caMKL1 = 38.7952653127559), tolerance = 1e-3)
    expect_equal(rowSums(resHuman$pearson.rhos$DECODED_CONTRASTS),
                 c(MKL1 = 20.7882406720692, caMKL1 = 103.779724182865), tolerance = 1e-3)
    expect_equal(rowSums(resHuman$pearson.rhos$RESIDUAL_CONTRASTS),
                 c(MKL1 = 4.64374405966305, caMKL1 = 4.33762865756607), tolerance = 1e-3)

    # zscores
    expect_type(resHuman$zscores, "list")
    expect_identical(lapply(resHuman$pearson.rhos, dimnames),
                     lapply(resHuman$zscores, dimnames))
    
    # TopHits
    expect_type(resHuman$TopHits, "list")
    expect_length(resHuman$TopHits, 3L)
    expect_identical(names(resHuman$TopHits), anms)
    expect_type(resHuman$TopHits$RESIDUAL_CONTRASTS, "list")
    expect_length(resHuman$TopHits$RESIDUAL_CONTRASTS, 2L)
    expect_identical(names(resHuman$TopHits$RESIDUAL_CONTRASTS), cnms)
    expect_s4_class(resHuman$TopHits$RESIDUAL_CONTRASTS[[cnms[1]]],
                    "DataFrame")
    expect_true(all(unlist(lapply(resHuman$TopHits, function(L) {
        unlist(lapply(L, nrow))
    })) == 10L))
})

