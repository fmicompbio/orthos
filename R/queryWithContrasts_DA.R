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
#' @param preserveInGlobalEnv specifies whether the contrast database (stored 
#'     as a \code{SummarizedExperiment} object in \code{target.contrasts})
#'     should be preserved in the GlobalEnv either for future queries or for 
#'     accessing its metadata
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
#' @importFrom DelayedMatrixStats rowSds rowMeans2
#' 
queryWithContrasts_DA <- function(contrasts = NULL, 
                               use = c("expressed.in.both", "all.genes"),
                               exprThr = 0.25, organism = c("Human","Mouse"), 
                               preserveInGlobalEnv = TRUE,
                               plotContrast = c("RESIDUAL", "INPUT", "DECODED", "NONE"),
                               detailTopn = 10, verbose = TRUE) {
    
    ## -------------------------------------------------------------------------
    ## Check inputs
    ## -------------------------------------------------------------------------
    stopifnot("`contrasts` should be a valid SummarizedExperiment" = 
                  is(contrasts, "SummarizedExperiment"))
    valid.contrasts <- c("INPUT_CONTRASTS", "DECODED_CONTRASTS", "RESIDUAL_CONTRASTS")
    present.contrasts <- intersect(valid.contrasts, names(assays(contrasts)))
    stopifnot("The assays slot in the provided SummarizedExperiment does not contain valid contrast names " = 
                  length(present.contrasts) > 0)
    message(paste ("provided contrast: ", present.contrasts, collapse = "\n"))
    
    use <- match.arg(use)
    organism <- match.arg(organism)
    plotContrast <- match.arg(plotContrast)
    
    ## -------------------------------------------------------------------------
    ## Load contrast database
    ## -------------------------------------------------------------------------
    if (!exists ("target.contrasts")) {
        if (verbose) {
            message("Loading contrast database...")
        }
        
        ### reworked implementation using DelayedMatrices.
        #target.contrasts <- loadHDF5SummarizedExperiment(dir="/tungstenfs/groups/gbioinfo/papapana/DEEP_LEARNING/Autoencoders/ARCHS4/Rdata/DECOMPOSED_CONTRASTS_HDF5/",prefix="human_v212")
        target.contrasts <- readRDS(system.file("data", paste0("DECOMPOSED_CONTRASTS_ARCHS4_v212_",organism,"_smpl.rds"), package = "deJUNKER"))
        
        if (preserveInGlobalEnv) {
            assign("target.contrasts", target.contrasts, envir = .GlobalEnv)
        }
    }
    
    stopifnot( "Incompatible rownames in the provided SummarizedExperiment.
Rownames should be the same as in the contrast database.
You can make sure by generating your SE generated using `decomposeVar`" = 
                   identical(rownames(contrasts), rownames(target.contrasts)))
    contrasts <- SummarizedExperiment::assays(contrasts)[present.contrasts]
    
    ## -------------------------------------------------------------------------
    ## Calculate correlations
    ## -------------------------------------------------------------------------
    pearson.rhos <- list()
    ## Calculate correlation
    for (ctr in present.contrasts){
        query <- contrasts[[ctr]]
        targt <- SummarizedExperiment::assays(target.contrasts)[[ctr]]
        pearson.rhos[[ctr]] <- do.call(cbind, mclapply(seq_len(ncol(query)), function(i) {
            idx <- which(!is.na(query[, i]))
            targtsub <- targt[idx, ]
            querysub <- query[idx, i, drop = FALSE]
            n_not_na <- t(!is.na(targtsub)) %*% matrix(1, nrow = length(idx), 
                                                       ncol = 1)
            query_sum_squared <- t(!is.na(targtsub)) %*% querysub ^ 2
            query_sum <- t(!is.na(targtsub)) %*% querysub
            
            targtsub <- scale(targtsub, center = TRUE, scale = TRUE)
            
            ## Set NA values to 0
            targtsub[is.na(targtsub)] <- 0
            
            ## Calculate correlation
            ((t(targtsub) %*% querysub) / (n_not_na - 1)) / 
                (sqrt(n_not_na * query_sum_squared - (query_sum) ^ 2) / 
                     sqrt(n_not_na * (n_not_na - 1)))
        }, mc.cores = 2L))
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
