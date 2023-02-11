## ------------------------------------------------------------------------- ##
## Checks, .loadContrastDatabase
## ------------------------------------------------------------------------- ##
test_that(".loadContrastDatabase works", {
    expect_error(.loadContrastDatabase("error"))
    
    seHuman <- .loadContrastDatabase("Human")
    seMouse <- .loadContrastDatabase("Mouse")
    
    expect_s4_class(seHuman, "SummarizedExperiment")
    expect_s4_class(seMouse, "SummarizedExperiment")
    
    assayNms <- c("INPUT_CONTRASTS", "DECODED_CONTRASTS",
                  "RESIDUAL_CONTRASTS", "CONTEXT")
    expect_identical(SummarizedExperiment::assayNames(seHuman),
                     c("INPUT_CONTRASTS", "DECODED_CONTRASTS",
                       "RESIDUAL_CONTRASTS", "CONTEXT"))
    expect_identical(SummarizedExperiment::assayNames(seMouse),
                     c("INPUT_CONTRASTS", "DECODED_CONTRASTS",
                       "RESIDUAL_CONTRASTS", "CONTEXT"))
    
    expect_identical(dim(seHuman), c(20411L, 74731L))
    expect_identical(dim(seMouse), c(20339L, 58532L))

    for (assayNm1 in assayNms) {
        expect_s4_class(SummarizedExperiment::assay(seHuman, assayNm1),
                        "DelayedArray")
        expect_s4_class(SummarizedExperiment::assay(seMouse, assayNm1),
                        "DelayedArray")
    }
})
