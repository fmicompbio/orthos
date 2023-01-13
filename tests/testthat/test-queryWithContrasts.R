## ------------------------------------------------------------------------- ##
## Checks, .loadContrastDatabase
## ------------------------------------------------------------------------- ##
test_that(".loadContrastDatabase works", {
    expect_error(.loadContrastDatabase("error"))
    
    seHuman <- .loadContrastDatabase("Human")
    seMouse <- .loadContrastDatabase("Mouse")
    
    expect_s4_class(seHuman, "SummarizedExperiment")
    expect_s4_class(seMouse, "SummarizedExperiment")
    
    expect_identical(SummarizedExperiment::assayNames(seHuman),
                     c("INPUT_CONTRASTS", "DECODED_CONTRASTS",
                       "RESIDUAL_CONTRASTS", "CONTEXT"))
    expect_identical(SummarizedExperiment::assayNames(seMouse),
                     c("INPUT_CONTRASTS", "DECODED_CONTRASTS",
                       "RESIDUAL_CONTRASTS", "CONTEXT"))
    
    expect_identical(dim(seHuman), c(20411L, 77016L))
    expect_identical(dim(seMouse), c(20339L, 59552L))

    expect_true(all(unlist(lapply(SummarizedExperiment::assays(seHuman),
                                  function(x) inherits(x, "DelayedArray")))))
    expect_true(all(unlist(lapply(SummarizedExperiment::assays(seMouse),
                                  function(x) inherits(x, "DelayedArray")))))
})
