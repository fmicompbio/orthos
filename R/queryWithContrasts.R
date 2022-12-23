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
#'     accurate, is much slower.
#' @param exprThr is the quantile in the provided context that determines the 
#'     expression value above which a gene is considered to be expressed. This 
#'     same value is then used for thresholding the contrast database. Only 
#'     applies when use="expressed.in.both"
#' @param organism Selects the autoencoder model trained on data from this
#'     species. One of \code{"Human"} or \code{"Mouse"}.
#' @param preserveInGlobalEnv specifies whether the contrast database (stored 
#'     as a \code{SummarizedExperiment} object in \code{target.contrasts})
#'     should be preserved in the \code{.GlobalEnv} environment either for
#'     future queries or for accessing its metadata
#' @param plotContrast Select a contrast to be plotted, one of
#'     \code{"RESIDUAL", "INPUT", "DECODED"} or \code{"NONE"} to suppress
#'     the plotting.
#' @param detailTopn specifies the number of top hits for which metadata will 
#'     be returned in the TopHits slot of the results.
#' @param verbose Logical scalar indicating whether to print messages along 
#'     the way.
#'
#' @return A list of PearsonRhos, Zscores against the datbase as well as 
#'     detailed Metadata for the detailTopn hits.
#' 
#' @importFrom digest digest
#' @importFrom SummarizedExperiment assays
#' @importFrom cowplot plot_grid
#' 
#' 
queryWithContrasts <- function(contrasts, 
                               use = c("expressed.in.both", "all.genes"),
                               exprThr = 0.25, organism = c("Human","Mouse"), 
                               preserveInGlobalEnv = TRUE,
                               plotContrast = c("RESIDUAL", "INPUT", "DECODED", "NONE"),
                               detailTopn = 10, verbose = TRUE) {
    
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
    .assertScalar(x = preserveInGlobalEnv, type = "logical")
    plotContrast <- match.arg(plotContrast)
    .assertScalar(x = detailTopn, type = "numeric", rngExcl = c(0, Inf))
    .assertScalar(x = verbose, type = "logical")
    
    ## -------------------------------------------------------------------------
    ## Load contrast database
    ## -------------------------------------------------------------------------
    if (!exists("target.contrasts", envir = .GlobalEnv)) {
        if (verbose) {
            message("Loading contrast database...")
        }
        target.contrasts <- readRDS(paste0("/tungstenfs/groups/gbioinfo/papapana/DEEP_LEARNING/Autoencoders/ARCHS4/Rdata/DECOMPOSED_CONTRASTS_ARCHS4_v212_",organism,".rds") )
        # target.contrasts <- readRDS(system.file("data", paste0("DECOMPOSED_CONTRASTS_ARCHS4_v212_",organism,"_smpl.rds"), package = "deJUNKER"))
        
        ### Probably better to move to HDF5 based implementation as this will 
        ### be both a significant speed-up in terms of loading
        ### and much more lean on memory requirements. However this needs a 
        ### reworked implementation using DelayedMatrices.
        #   target.contrasts <- loadHDF5SummarizedExperiment(dir="/tungstenfs/groups/gbioinfo/papapana/DEEP_LEARNING/Autoencoders/ARCHS4/Rdata/DECOMPOSED_CONTRASTS_HDF5/",
        #                                               prefix="human_v212")
        if (preserveInGlobalEnv) {
            assign("target.contrasts", target.contrasts, envir = .GlobalEnv)
        }
    }
    # A sample of col indices to quickly check data integrity
    smpl.col <- if (ncol(target.contrasts) == 72) 1:72 else seq(1, 50000, 700)
    DBhash <- digest::digest(target.contrasts[, smpl.col], algo = "xxhash64")
    hashvals <- list(Human = "803df1fc71c30ba6", Mouse = "ffd499bd8c1060ec")
    stopifnot("The contrast DB contained in the `target.contrasts` object has not been correctly loaded.
Please remove `target.contrasts` and try again." = 
                  DBhash == hashvals[[organism]])
    
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
        set.to.NA <- SummarizedExperiment::assays(target.contrasts)[["CONTEXT"]] <= thr
        
        pearson.rhos <- sapply(present.contrasts, function(x) {
            message(paste0("Querying contrast database with ", x, "..."))
            Thresholded.Contrast <- SummarizedExperiment::assays(target.contrasts)[[x]]
            Thresholded.Contrast[ set.to.NA ] <- NA
            stats::cor(contrasts[[x]], Thresholded.Contrast, 
                       use = "pairwise.complete.obs") 
        }, simplify = FALSE, USE.NAMES = TRUE
        )
    } else if (use == "all.genes") {
        pearson.rhos <- sapply(present.contrasts, function(x) {
            message(paste0("Querying contrast database with ", x, "..."))
            stats::cor(contrasts[[x]], SummarizedExperiment::assays(target.contrasts)[[x]])  
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
            colData(target.contrasts)[N, c(3, 12, 22, 24, 29, 31, 33)]
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
