#' Load contrast database
#' 
#' Load a pre-calculated, organism-specific contrast database and return it
#' as a \code{SummarizedExperiment}.
#'
#' @param organism Character scalar selecting the organism for which to load the
#'     contrast database. One of \code{"Human"} or \code{"Mouse"}.
#'
#' @return A \code{SummarizedExperiment} with pre-calculate contrasts as
#'     assays.
#'
#' @author Panagiotis Papasaikas, Michael Stadler
#'
#' @importFrom HDF5Array loadHDF5SummarizedExperiment
# #' @importFrom digest digest
#' 
#' @keywords internal
.loadContrastDatabase <- function(organism = c("Human", "Mouse")) {
    organism <- match.arg(organism)

    # load SummarizedExperiment    
    # currently, this loads a local file
    # in the future, this will obtain the database using BiocFileCache or ExperimentHub
    se <- HDF5Array::loadHDF5SummarizedExperiment(dir = "/tungstenfs/groups/gbioinfo/papapana/DEEP_LEARNING/Autoencoders/ARCHS4/Rdata/DECOMPOSED_CONTRASTS_HDF5",
                                       prefix = paste0(tolower(organism),"_v212_c100"))
    
    # check validity
    # DBhash <- digest::digest(se, algo = "xxhash64")
    # hashvals <- list(Human = "87e231a6567c61e0", Mouse = "800d2113e4a41175" )
    #    stopifnot(paste0("The contrast DB contained for ", organism,
    #                     " has not been correctly loaded. Please remove it and try again.") = 
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
#' @param organism Selects the autoencoder model trained on data from this
#'     species. One of \code{"Human"} or \code{"Mouse"}.
#' @param plotContrast Select a contrast to be plotted, one of
#'     \code{"RESIDUAL", "INPUT", "DECODED"} or \code{"NONE"} to suppress
#'     the plotting.
#' @param detailTopn specifies the number of top hits for which metadata will 
#'     be returned in the TopHits slot of the results.
#' @param verbose Logical scalar indicating whether to print messages along 
#'     the way.
#' @param BPPARAM BiocParallelParam object specifying how parallelization is to be performed using for e.g
#'     \code{\link[BiocParallel]{MulticoreParam}}) or \code{\link[BiocParallel]{SnowParam}})
#' @param chunk_size Column dimension for the grid used to read blocks from the
#'     HDF5 Matrix. Sizes between 250 and 1000 are recommended. Smaller sizes
#'     reduce memory usage.
#'
#' @return A list of PearsonRhos, Zscores against the datbase as well as 
#'     detailed Metadata for the detailTopn hits.
#' 
#' @importFrom SummarizedExperiment assays colData
#' @importFrom parallel detectCores
#' @importFrom BiocParallel bpparam bpisup bpstart bpstop
#' @importFrom cowplot plot_grid
#' @importFrom stats quantile
#' 
queryWithContrasts <- function(contrasts,
                               use = c("expressed.in.both", "all.genes"),
                               exprThr = 0.25, organism = c("Human","Mouse"),
                               plotContrast = c("RESIDUAL", "INPUT", "DECODED", "NONE"),
                               detailTopn = 10, verbose = TRUE,
                               BPPARAM = BiocParallel::bpparam(),
                               chunk_size = 500 ) {
    
    ## -------------------------------------------------------------------------
    ## Check inputs
    ## -------------------------------------------------------------------------
    .assertVector(x = contrasts, type = "SummarizedExperiment")
    valid.contrasts <- c("INPUT_CONTRASTS", "DECODED_CONTRASTS", "RESIDUAL_CONTRASTS")
    present.contrasts <- intersect(valid.contrasts, names(assays(contrasts)))
    stopifnot("The assays slot in the provided SummarizedExperiment does not contain valid contrast names " = 
                  length(present.contrasts) > 0)
    message(paste("provided contrast: ", present.contrasts, collapse = "\n"))
    
    use <- match.arg(use)
    .assertScalar(x = exprThr, type = "numeric", rngIncl = c(0, 1))
    organism <- match.arg(organism)
    plotContrast <- match.arg(plotContrast)
    .assertScalar(x = detailTopn, type = "numeric", rngExcl = c(0, Inf))
    .assertScalar(x = verbose, type = "logical")
    .assertScalar(x = chunk_size, type = "numeric", rngIncl = c(100, Inf))
    stopifnot("BPPARAM must be an instance of a BiocParallelParam class e.g., MulticoreParam, SnowParam or SerialParam" =
                identical(attr(class(BPPARAM), "package"), "BiocParallel") )
    if ( !BiocParallel::bpisup(BPPARAM) && !is(BPPARAM, "MulticoreParam")     ) {
      BiocParallel::bpstart(BPPARAM)
      on.exit(BiocParallel::bpstop(BPPARAM), add=TRUE)
    }
    
    ## -------------------------------------------------------------------------
    ## Load contrast database
    ## -------------------------------------------------------------------------
    if (verbose) {
        message("Loading contrast database...")
    }
    target.contrasts <- .loadContrastDatabase(organism = organism)
    
    stopifnot( "Incompatible rownames in the provided SummarizedExperiment.
Rownames should be the same as in the contrast database.
You can make sure by generating your SE generated using `decomposeVar`" = 
                   identical(rownames(contrasts), rownames(target.contrasts)))
    contrasts <- SummarizedExperiment::assays(contrasts)[present.contrasts]
    
    ## -------------------------------------------------------------------------
    ## Calculate correlations
    ## -------------------------------------------------------------------------
    if (use == "expressed.in.both") {
        # Set a global expression threshold according to a quantile in the query data context
        message("Thresholding genes...")
        thr <- quantile(contrasts[["CONTEXT"]], exprThr)
        
        pearson.rhos <- sapply(present.contrasts, function(x) {
            message(paste0("Querying contrast database with ", x, "..."))
            query <- contrasts[[x]]
            query[query <= thr] <- NA
            .grid_cor_wNAs(query, hdf5 = SummarizedExperiment::assays(target.contrasts)[[x]], 
                           thr = thr, BPPARAM = BPPARAM, chunk_size = chunk_size) 
        }, simplify = FALSE, USE.NAMES = TRUE
        )
        
    } else if (use == "all.genes") {
        pearson.rhos <- sapply(present.contrasts, function(x) {
            message(paste0("Querying contrast database with ", x, "..."))
            query <- contrasts[[x]]
            .grid_cor_woNAs(query, hdf5 = SummarizedExperiment::assays(target.contrasts)[[x]], 
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
    zscores <- sapply(present.contrasts, function(x) {
        t(scale(t(pearson.rhos[[x]])))
    }, simplify = FALSE, USE.NAMES = TRUE
    )
    
    ## -------------------------------------------------------------------------
    ## Get top hits
    ## -------------------------------------------------------------------------
    TopHits <- sapply(present.contrasts, function(contr) {
        apply(zscores[[contr]], 1, function(x) { 
            N <- names(sort(x, decreasing = TRUE)[seq_len(detailTopn)])
            SummarizedExperiment::colData(target.contrasts)[N, c(3, 12, 22, 24, 29, 31, 33)]
        })
    }, simplify = FALSE, USE.NAMES = TRUE
    )
    
    ## -------------------------------------------------------------------------
    ## Plot
    ## -------------------------------------------------------------------------
    if (plotContrast != "NONE" && 
        paste0(plotContrast, "_CONTRASTS") %in% present.contrasts) {
        message("Generating plots...")
        plot.data <- zscores[[paste0(plotContrast, "_CONTRASTS")]]
        PLOTLIST <- sapply(seq_len(nrow(plot.data)), function(i) {
            plotQueryResults(plot.data[i, ], annot = target.contrasts$series_id,
                             topn = detailTopn)
        }, simplify = FALSE)
        suppressWarnings({
            print(
                cowplot::plot_grid(plotlist = PLOTLIST,
                                   labels = rownames(plot.data))
            )
        })
    }
    
    if (verbose) {
        message("Done!")
    }
    
    return(list(pearson.rhos = pearson.rhos, zscores = zscores,
                TopHits = TopHits))
}
