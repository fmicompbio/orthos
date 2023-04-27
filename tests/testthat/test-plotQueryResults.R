.createSyntheticQuery <- function(n = 100L,
                                  ntop = 10L,
                                  ids = sprintf("sample%03d", seq.int(n))) {
    tmpI <- matrix(rnorm(2 * n, mean = 0, sd = 2), nrow = 2L,
                   dimnames = list(c("treat1", "treat2"), ids))
    tmpD <- matrix(rnorm(2 * n, mean = 0, sd = 1), nrow = 2L,
                   dimnames = list(c("treat1", "treat2"), ids))
    tmpR <- tmpI - tmpD
    idxtopI <- sample(x = n, size = ntop)
    idxtopD <- sample(x = n, size = ntop)
    idxtopR <- sample(x = n, size = ntop)
    tmpTopI <- S4Vectors::DataFrame(
        characteristics_ch1 = rep("characteristics", ntop),
        geo_accession = ids[idxtopI],
        series_id = sprintf("GSE%05d", round(seq(1, ntop) / 4) + 1),
        source_name_ch1 = rep("source", ntop),
        title = rep("title", ntop),
        Cor2CNT = runif(ntop, 0.6, 1.0),
        CNTname = ids[idxtopI])
    tmpTopD <- tmpTopR <- tmpTopI
    rownames(tmpTopI) <- tmpTopI$geo_accession
    tmpTopD$geo_accession <- tmpTopD$CNTname <- ids[idxtopD]
    rownames(tmpTopD) <- tmpTopD$geo_accession
    tmpTopR$geo_accession <- tmpTopR$CNTname <- ids[idxtopR]
    rownames(tmpTopR) <- tmpTopR$geo_accession
    qres <- list(
        pearson.rhos = list(INPUT_CONTRASTS = tmpI,
                            DECODED_CONTRASTS = tmpD,
                            RESIDUAL_CONTRASTS = tmpR),
        zscores = list(INPUT_CONTRASTS = t(scale(t(tmpI))),
                       DECODED_CONTRASTS = t(scale(t(tmpD))),
                       RESIDUAL_CONTRASTS = t(scale(t(tmpR)))),
        TopHits = list(INPUT_CONTRASTS = list(treat1 = tmpTopI,
                                              treat2 = tmpTopI),
                       DECODED_CONTRASTS = list(treat1 = tmpTopD,
                                                treat2 = tmpTopD),
                       RESIDUAL_CONTRASTS = list(treat1 = tmpTopR,
                                                 treat2 = tmpTopR)))
    return(qres)
}


## ------------------------------------------------------------------------- ##
## Checks, plotQueryResultsManh
## ------------------------------------------------------------------------- ##
test_that("plotQueryResultsManh works", {
    # create synthetic query result
    set.seed(42L)
    qres <- .createSyntheticQuery(n = 50L, ntop = 5L)
    
    # arguments
    expect_error(plotQueryResultsManh(queryResults = "error"))
    expect_error(plotQueryResultsManh(queryResults = qres, doPlot = "error"))
    
    # results
    tf <- tempfile(fileext = ".png")
    grDevices::png(filename = tf, width = 1000, height = 1000, pointsize = 16)
    p <- plotQueryResultsManh(queryResults = qres, doPlot = TRUE)
    grDevices::dev.off()
    
    expect_type(p, "list")
    expect_length(p, nrow(qres$pearson.rhos$INPUT_CONTRASTS))
    expect_named(p, rownames(qres$pearson.rhos$INPUT_CONTRASTS))
    for (i in seq_along(p)) {
        expect_s3_class(p[[i]], "ggplot")
    }
    expect_true(file.info(tf)['size'] > 5000)
    
    # clean up
    unlink(tf)
})

## ------------------------------------------------------------------------- ##
## Checks, plotQueryResultsViolin
## ------------------------------------------------------------------------- ##
test_that("plotQueryResultsViolin works", {
    # create synthetic query result
    set.seed(41L)
    qres <- .createSyntheticQuery(n = 50L, ntop = 5L)
    
    # arguments
    expect_error(plotQueryResultsViolin(queryResults = "error"))
    expect_error(plotQueryResultsViolin(queryResults = qres, doPlot = "error"))
    
    # results
    tf <- tempfile(fileext = ".png")
    grDevices::png(filename = tf)
    p <- plotQueryResultsViolin(queryResults = qres, doPlot = TRUE)
    grDevices::dev.off()
    
    expect_type(p, "list")
    expect_length(p, nrow(qres$pearson.rhos$INPUT_CONTRASTS))
    expect_named(p, rownames(qres$pearson.rhos$INPUT_CONTRASTS))
    for (i in seq_along(p)) {
        expect_s3_class(p[[i]], "ggplot")
    }
    expect_true(file.info(tf)['size'] > 20000)
    
    # clean up
    unlink(tf)
})
