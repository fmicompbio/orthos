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
#'
#' @return A list of PearsonRhos, Zscores against the datbase as well as 
#'     detailed Metadata for the detailTopn hits.
#' 
queryWithContrasts <- function(contrasts = NULL, 
                               use = c("expressed.in.both", "all.genes"),
                               exprThr = 0.25, organism = "Human", 
                               preserveInGlobalEnv = TRUE,
                               plotContrast = c("RESIDUAL", "INPUT", "DECODED", "NONE"),
                               detailTopn = 10) {
    
    stopifnot("`contrasts` should be a valid SummarizedExperiment" = is(contrasts, "SummarizedExperiment"))
    valid.contrasts <- c("INPUT_CONTRASTS", "DECODED_CONTRASTS", "RESIDUAL_CONTRASTS")
    present.contrasts <- intersect(valid.contrasts, names(assays(contrasts)))
    stopifnot("The assays slot in the provided SummarizedExperiment does not contain valid contrast names " = length(present.contrasts) > 0)
    message(paste ("provided contrast: ", present.contrasts, collapse = "\n"))
    
    plotContrast <-  match.arg(plotContrast)
    
    smpl.col <- seq(1, 75000, 100) # A sample of row/col indexes to quickly check data integrity
    smpl.row <- seq(1, 20000, 100)
    if (!exists ("target.contrasts")) {
        message("Loading contrast database...")
        target.contrasts <- readRDS("/tungstenfs/groups/gbioinfo/papapana/DEEP_LEARNING/Autoencoders/ARCHS4/Rdata/DECOMPOSED_CONTRASTS_human_v212_uncompressed.rds")
        ### Probably better to move to HDF5 based implementation as this will 
        ### be both a significant speed-up in terms of loading
        ### and much more lean on memory requirements. However this needs a 
        ### reworked implementation using DelayedMatrices.
        #   target.contrasts <- loadHDF5SummarizedExperiment(dir="/tungstenfs/groups/gbioinfo/papapana/DEEP_LEARNING/Autoencoders/ARCHS4/Rdata/DECOMPOSED_CONTRASTS_HDF5/",
        #                                               prefix="human_v212")
        if (preserveInGlobalEnv){
            assign("target.contrasts", target.contrasts, envir = .GlobalEnv)
        }
    }
    
    DBhash <- digest::digest(target.contrasts[smpl.row,smpl.col], algo = "xxhash64")
    
    stopifnot("The contrast DB contained in the `target.contrasts` object has not been correctly loaded.
Please remove `target.contrasts` and try again." = DBhash == "f9abc421d4e93c08")
    
    stopifnot( "Incompatible rownames in the provided SummarizedExperiment.
Rownames should be the same as in the contrast database.
You can make sure by generating your SE generated using `decomposeVar`"
               = identical(rownames(contrasts), rownames(target.contrasts)))
    contrasts <- assays(contrasts)[present.contrasts]
    
    
    if (use == "expressed.in.both") {
        #Set a global expression threshold according to a quantile in the query data context
        message("Thresholding genes...")
        thr <- quantile(contrasts[["CONTEXT"]], exprThr)
        set.to.NA <- assays(target.contrasts)[["CONTEXT"]] <= thr
        
        pearson.rhos <- sapply(present.contrasts, function(x) {
            message(paste0("Querying contrast database with ", x, "..."))
            Thresholded.Contrast <- assays(target.contrasts)[[x]]
            Thresholded.Contrast[ set.to.NA ] <- NA
            cor(contrasts[[x]], Thresholded.Contrast, 
                use = "pairwise.complete.obs") 
            }, simplify = FALSE, USE.NAMES = TRUE
        )
    } else if (use == "all.genes") {
        pearson.rhos <- sapply(present.contrasts, function(x) {
            message(paste0("Querying contrast database with ", x, "..."))
            cor(contrasts[[x]], assays(target.contrasts)[[x]])  
            }, simplify = FALSE, USE.NAMES = TRUE
        )
    }
    
    ### Z-score and pval calculation for each decomposed component query:
    message("Compiling query statistics...")
    zscores <- sapply(present.contrasts, function(x) {
        t(scale(t(query.res$pearson.rhos[[x]])))
    }, simplify = FALSE, USE.NAMES = TRUE
    )
    
    TopHits <- sapply(present.contrasts, function(contr){
        apply(zscores[[contr]], 1, function(x) { 
            N <- names(sort(x, decreasing = TRUE)[seq_len(detailTopn)])
        colData(target.contrasts)[N, c(3, 12, 22, 24, 29, 31, 33)]
        })
    }, simplify = FALSE, USE.NAMES = TRUE
    )
    
    if (plot.contrast != "NONE" && 
        paste0(plot.contrast, "_CONTRASTS") %in% present.contrasts) {
        message("Generating plots...")
        plot.data <- zscores[[paste0(plot.contrast, "_CONTRASTS")]]
        PLOTLIST <- sapply(seq_len(nrow(plot.data)), function(i) {
            plot.query_results(plot.data[i, ])
        }, simplify = FALSE)
        suppressWarnings({
            print(
                cowplot::plot_grid(plotlist = PLOTLIST,
                                   labels = rownames(plot.data))
            )
        })
    }
    
    message("Done!")
    return(list(pearson.rhos = pearson.rhos, zscores = zscores,
                TopHits = TopHits))
}
