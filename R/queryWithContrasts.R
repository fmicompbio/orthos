#' Load contrast database
#'
#' Load a pre-calculated, organism-specific contrast database and return it
#' as a \code{SummarizedExperiment}.
#'
#' @param organism Character scalar selecting the organism for which to load the
#'     contrast database. One of \code{"Human"} or \code{"Mouse"}.
#' @param mode When in "ANALYSIS" mode (default) the complete contrast DB is
#'     queried. "DEMO" mode employs a small "toy" database for the queries.
#'     "DEMO" should only be used for testing/demonstration purposes
#'     and never for actual analysis purposes.
#' @param mustWork Logical scalar. If \code{FALSE} and the contrast database
#'     is not available, return an empty \code{SummarizedExperiment} object.
#'     If \code{TRUE} (the default) and the contrast database is
#'     not available, \code{.loadContrastDatabase} throws an error.
#'
#' @return A \code{SummarizedExperiment} with pre-calculated contrasts as
#'     assays.
#'
#' @author Panagiotis Papasaikas, Michael Stadler
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom HDF5Array loadHDF5SummarizedExperiment
# #' @importFrom digest digest
#'
#' @keywords internal
#' @noRd
.loadContrastDatabase <- function(organism = c("Human", "Mouse"),
                                  mode = c("ANALYSIS", "DEMO"),
                                  mustWork = TRUE) {
    organism <- match.arg(organism)
    mode <- match.arg(mode)

    # load SummarizedExperiment
    # currently, this loads a local file
    # in the future, this will obtain the database using BiocFileCache or
    # ExperimentHub
    
    # dataDir <- "/tungstenfs/groups/gbioinfo/papapana/DEEP_LEARNING/Autoencoders/ARCHS4/Rdata/DECOMPOSED_CONTRASTS_HDF5"

    dataDir <- orthosData::GetorthosContrastDB(organism=organism, mode=mode)
    if (identical(mode, "DEMO")) {
        prefix <- paste0(tolower(organism), "_v212_NDF_c100_DEMO")
        seFile <- file.path(dataDir, paste0(prefix, "se.rds"))
    } else {
        prefix <- paste0(tolower(organism), "_v212_NDF_c100")
        seFile <- file.path(dataDir, paste0(prefix, "se.rds"))
    }
    
    if (file.exists(seFile)) {
        se <- HDF5Array::loadHDF5SummarizedExperiment(dir = dataDir,
                                                      prefix = prefix)
    } else {
        if (mustWork) {
            stop("contrast database for '", organism, "' is not available")
        }
        se <- SummarizedExperiment::SummarizedExperiment()
    }
    
    # check validity
    # DBhash <- digest::digest(se, algo = "xxhash64")
    # hashvals <- list(Human = "87e231a6567c61e0", Mouse = "800d2113e4a41175" )
    #    stopifnot(paste0("The contrast DB contained for ", organism,
    #                     " has not been correctly loaded. ",
    #                     "Please remove it and try again.") =
    #                  DBhash == hashvals[[organism]])

    # return
    return(se)
}


#' Query the contrast database with a set of contrasts
#'
#' @export
#' @author Panagiotis Papasaikas
#'
#' @param contrasts A \code{SummarizedExperiment} object with assays containing
#'     contrasts named INPUT_CONTRASTS, DECODED_CONTRASTS and RESIDUAL_CONTRASTS
#'     (at least one should be present) and context information in an assay
#'     named CONTEXT. The latter is only required when use="expressed.in.both".
#'     This is typically generated using
#'     \code{decomposeVar}.
#' @param use Determines if all.genes or genes expressed in both query and
#'     target context will be used. Note that "expressed.in.both", though more
#'     accurate, is slower.
#' @param exprThr is the quantile in the provided context that determines the
#'     expression value above which a gene is considered to be expressed. This
#'     same value is then used for thresholding the contrast database. Only
#'     applies when use="expressed.in.both"
#' @param organism Uses the `orthosData` contrast Database  from this
#'     species. One of \code{"Human"} or \code{"Mouse"}.
#' @param plotType Select the type of visualization for the query results
#'     \code{"violin", "manh"} or \code{"none"} to suppress
#'     the plotting.
#' @param detailTopn specifies the number of top hits for which metadata will
#'     be returned in the TopHits slot of the results.
#' @param verbose Logical scalar indicating whether to print messages along
#'     the way.
#' @param BPPARAM BiocParallelParam object specifying how parallelization is to
#'     be performed using e.g. \code{\link[BiocParallel]{MulticoreParam}})
#'     or \code{\link[BiocParallel]{SnowParam}})
#' @param chunk_size Column dimension for the grid used to read blocks from the
#'     HDF5 Matrix. Sizes between 250 and 1000 are recommended. Smaller sizes
#'     reduce memory usage.
#' @param mode When in "ANALYSIS" mode (default) the complete contrast DB is
#'     queried. "DEMO" mode employs a small "toy" database for the queries.
#'     "DEMO" should only be used for testing/demonstration purposes
#'     and never for actual analysis purposes.
#'
#' @return A list with three elements called "pearson.rhos", "zscores" and
#'     "TopHitsof", containing raw and z-scored Pearson's rho correlation
#'     coefficients between the query contrast(s) and the contrasts in the
#'     database, as well as detailed metadata for the \code{detailTopn} best
#'     hits.
#'
#' @importFrom SummarizedExperiment assays colData
#' @importFrom parallel detectCores
#' @importFrom BiocParallel bpparam bpisup bpstart bpstop bpprogressbar
#'     bpworkers bptasks
#' @importFrom cowplot plot_grid
#' @importFrom stats quantile
#'
queryWithContrasts <- function(contrasts,
                               use = c("expressed.in.both", "all.genes"),
                               exprThr = 0.25,
                               organism = c("Human", "Mouse"),
                               plotType = c("violin", "manh", "none"),
                               detailTopn = 10,
                               verbose = TRUE,
                               BPPARAM = BiocParallel::bpparam(),
                               chunk_size = 500,
                               mode = c("ANALYSIS", "DEMO")) {

    ## -------------------------------------------------------------------------
    ## Check inputs
    ## -------------------------------------------------------------------------
    .assertVector(x = contrasts, type = "SummarizedExperiment")
    validContrasts <- c("INPUT_CONTRASTS", "DECODED_CONTRASTS",
                        "RESIDUAL_CONTRASTS")
    presentContrasts <- intersect(validContrasts, names(assays(contrasts)))
    stopifnot("The assays slot in the provided SummarizedExperiment does not contain valid contrast names " =
                  length(presentContrasts) > 0)

    if (verbose) {
        message(paste("provided contrast: ", presentContrasts,
                      collapse = "\n"))
    }
    use <- match.arg(use)
    organism <- match.arg(organism)
    plotType <- match.arg(plotType)
    mode <- match.arg(mode)
    .assertScalar(x = exprThr, type = "numeric", rngIncl = c(0, 1))
    .assertScalar(x = detailTopn, type = "numeric", rngExcl = c(0, Inf))
    .assertScalar(x = verbose, type = "logical")
    .assertScalar(x = chunk_size, type = "numeric", rngIncl = c(100, Inf))
    .assertScalar(x = BPPARAM, type = "BiocParallelParam")
    
    ## ------------------------------------------------------------------------
    ## Setup bpprogressbar, initialize cluster
    ## ------------------------------------------------------------------------
    if (verbose) {
        BiocParallel::bpprogressbar(BPPARAM) <- TRUE
        if (BiocParallel::bpworkers(BPPARAM) < 10) {
            BiocParallel::bptasks(BPPARAM) <- ifelse(
                identical(class(BPPARAM)[[1]],"SerialParam"),
                10, min(2 * BiocParallel::bpworkers(BPPARAM), 10))
        }
    }
    if (!BiocParallel::bpisup(BPPARAM) && !is(BPPARAM, "MulticoreParam")) {
        BiocParallel::bpstart(BPPARAM)
        on.exit(BiocParallel::bpstop(BPPARAM), add = TRUE)
    }

    ## -------------------------------------------------------------------------
    ## Load contrast database
    ## -------------------------------------------------------------------------
    if (verbose) {
        message("Loading contrast database...")
    }
    targetContrasts <- .loadContrastDatabase(organism = organism, mode)

    stopifnot( "Incompatible rownames in the provided SummarizedExperiment.
Rownames should be the same as in the contrast database.
You can make sure by generating your SE generated using `decomposeVar`" =
                   identical(rownames(contrasts), rownames(targetContrasts)))
    context <- SummarizedExperiment::assays(contrasts)[["CONTEXT"]]
    contrasts <- SummarizedExperiment::assays(contrasts)[presentContrasts]

    ## -------------------------------------------------------------------------
    ## Calculate correlations
    ## -------------------------------------------------------------------------
    if (use == "expressed.in.both") {
        # Set a global expression threshold according to a quantile in the
        # query data context
        if (verbose) {
            message("Thresholding genes...")
        }
        thr <- stats::quantile(context, exprThr)

        pearson.rhos <- sapply(presentContrasts, function(x) {
            if (verbose) {
                message("Querying contrast database with ", x, "...")
            }
            query <- contrasts[[x]]
            query[context <= thr] <- NA
            .grid_cor_wNAs(
                query,
                hdf5 = SummarizedExperiment::assays(targetContrasts)[[x]],
                hdf5_ctx = SummarizedExperiment::assays(targetContrasts)[["CONTEXT"]],
                thr = thr, BPPARAM = BPPARAM, chunk_size = chunk_size)
        }, simplify = FALSE, USE.NAMES = TRUE
        )
    } else if (use == "all.genes") {
        pearson.rhos <- sapply(presentContrasts, function(x) {
            if (verbose) {
                message("Querying contrast database with ", x, "...")
            }
            query <- contrasts[[x]]
            .grid_cor_woNAs(
                query,
                hdf5 = SummarizedExperiment::assays(targetContrasts)[[x]],
                BPPARAM = BPPARAM, chunk_size = chunk_size)
        }, simplify = FALSE, USE.NAMES = TRUE
        )
    }

    ## -------------------------------------------------------------------------
    ## Calculate z-scores for each decomposed component query
    ## -------------------------------------------------------------------------
    if (verbose) {
        message("Compiling query statistics...")
    }
    zscores <- sapply(presentContrasts, function(x) {
        t(scale(t(pearson.rhos[[x]])))
    }, simplify = FALSE, USE.NAMES = TRUE
    )

    ## -------------------------------------------------------------------------
    ## Get top hits
    ## -------------------------------------------------------------------------
    TopHits <- sapply(presentContrasts, function(contr) {
        apply(zscores[[contr]], 1, function(x) {
            #N <- names(sort(x, decreasing = TRUE)[seq_len(detailTopn)])
            Zscore <- sort(x, decreasing = TRUE)[seq_len(detailTopn)]
            N <- names(Zscore)
            DBinfo <- SummarizedExperiment::colData(targetContrasts)[N, c(12, 29, 3, 22, 24,
                                                                          33, 31)]
            colnames(DBinfo)[colnames(DBinfo)=="CNTname"] <- "CNT_geo_accession"
            colnames(DBinfo)[colnames(DBinfo)=="geo_accession"] <- "TREATM_geo_accession"
            colnames(DBinfo)[colnames(DBinfo)=="Cor2CNT"] <- "corr_TREATM_CNT"
            cbind(Zscore, DBinfo)
            
        })
    }, simplify = FALSE, USE.NAMES = TRUE
    )

    ## -------------------------------------------------------------------------
    ## Gather results
    ## -------------------------------------------------------------------------
    RESULTS <- list(pearson.rhos = pearson.rhos,
                    zscores = zscores,
                    TopHits = TopHits)

    ## -------------------------------------------------------------------------
    ## Plot
    ## -------------------------------------------------------------------------
    if (plotType != "none") {
        if (verbose) {
            message("Generating plots...")
        }

        suppressWarnings({
            if (plotType == "violin") {
                PLOTS <- plotQueryResultsViolin(RESULTS)
            } else if (plotType == "manh") {
                PLOTS <- plotQueryResultsManh(RESULTS)
            }
        })
    }

    if (verbose) {
        message("Done!")
    }

    return(RESULTS)
}
