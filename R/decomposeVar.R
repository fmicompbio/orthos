#' Decompose input contrasts to decoded and residual fractions
#' 
#' Decompose input contrasts to decode and residual fractions according to a 
#' trained contrast encoder-decoder.
#' 
#' @author Panagiotis Papasaikas
#' @export
#' 
#' @param M Matrix of raw gene counts.
#' @param treatm,cntr Vectors indicating column indices in M corresponding 
#'     to treatments and controls. If treatm and cntr are specified MD has to be 
#'     \code{NULL}. Alternatively (if MD is specified), M is assumed to be a 
#'     raw gene count matrix specifying context for contrasts specified in MD.
#'     MD is then a matrix of gene deltas with the same dimensions as M. 
#'     If MD is specified treatm and cntr have to be \code{NULL}.
#' @param processInput If set to \code{TRUE} (default) the count matrix will 
#'     be preprocessed (library normalized, log-transformed, NA values will 
#'     be set to 0).
#' @param featureType: Set to AUTO for automatic feature id-type detection. 
#'     Alternatively specify the type of supplied id features.
#'
#' @return A \code{SummarizedExperiment} object with the decomposed contrasts 
#'     in the assays and the decomposed variance as the \code{colData}.
#' 
#' @examples
#' RES <- decomposeVar(M = ContextM, MD = DeltaM, processInput = FALSE)
#' RES <- decomposeVar(M = CountM, treatm = c(3, 4, 5, 6), 
#'                     cntr = c(1, 1, 2, 2), processInput = FALSE)
#' 
decomposeVar <- function(M, MD = NULL, treatm = NULL, cntr = NULL, 
                         processInput = TRUE, organism = "Human",
                         featureType = c("AUTO", "ENSEMBL_GENE_ID",
                                         "GENE_SYMBOL", "ENTREZ_GENE_ID",
                                         "ARCHS4_ID")) {
    
    featureType <- match.arg(featureType)
    
    message("Checking input...")
    stopifnot("`M` must be a (raw) count matrix with features (genes) in the rows and conditions in the columns" = is.matrix(M) & all(M >= 0, na.rm = TRUE))
    stopifnot("`specify either both `treatm` and `cntr` OR  `MD`" = !is.null(MD) | (!is.null(treatm) & !is.null(cntr)))
    
    if (is.null(MD)) {
        stopifnot("`treatm` and `cntr` vectors must have the same length" = identical(length(treatm), length(cntr)))
        stopifnot("`treatm` and `cntr` vectors must be numeric vectors indicating `M` column indices  " = all(c(treatm, cntr) %in% 1:ncol(M)) & !is.null(treatm))
    }
    
    if (!is.null(MD)) {
        stopifnot("`treatm` and `cntr` arguments must be NULL when `MD` is specified" = is.null(treatm) & is.null(cntr))
        stopifnot("`M` and `MD` matrices have to have the same dimensions and dimnames" = identical( dimnames(M), dimnames(MD)))
    }
    
    #### Feature ID type detection:
    geneIDs <- rownames(M)
    rownames(M) <- toupper(gsub("\\.\\d.+","",rownames(M)))
    if (featureType == "AUTO") {
        message("Detecting feature ids-type...")
        IDtypes <- c("ENSEMBL_GENE_ID", "GENE_SYMBOL", "ENTREZ_GENE_ID", "ARCHS4_ID")
        ID_OLaps <- sapply(IDtypes, function(x) {  
            sum(rownames(M) %in% toupper(genes[,x])) 
        })
        featureType <- IDtypes[which.max(ID_OLaps)]
        message("Feature ids-type detected: ", featureType)
    } else {
        message("Feature ids-type specified by user: ", featureType)
        ID_OLaps <- sum(rownames(M) %in% toupper(genes[,featureType]))
    }
    
    message(max(ID_OLaps), "/", nrow(M) , 
            " provided input features mapped against a total of ", 
            nrow(genes), " model features.")
    message(nrow(genes) - max(ID_OLaps), " missing features will be set to 0.")
    message("--> Missing features corresponding to non/lowly expressed genes ", 
            "in your context(s) are of no consequence.")
    #message("--> A small fraction (<10%) of missing genes, expressed in your context(s), are well-tolerated.")
    message("--> The model is robust to small fractions (<10%) of missing ", 
            "genes that are expressed in your context(s).")
    message("--> Increased numbers of missing expressed genes in your input ", 
            "might result in model performance decline.")
    stopifnot("\n!!!Too many missing features (>10000).  Make certain you provide valid human
Symbol, Ensembl or Entrez gene identifiers as rownames in your input matrix `M`.\n" = (nrow(genes) - max(ID_OLaps) <= 10000))
    idx.commonF <- na.omit(match(rownames(M), toupper(genes[, featureType])))
    idx.commonR <- which(rownames(M) %in% toupper(genes[, featureType]))
    
    ########## Preprocess input (NAs ->0, LibNormalize, LogTransform):
    pc <- 4
    if (processInput) {
        message("Preparing input...")
        #### Manage NA values:
        NA.idx <- which(is.na(M))
        if (length(NA.idx) >0){
            message("!!!Matrix `M` contains ", length(NA.idx),
                    " NA values.", " Those will be set to 0")
            M[NA.idx] <- 0
        }
        ##### Lib normalize and log transform input count matrix:
        M <- sweep(M[rownames(M) %in% toupper(genes[, featureType]), ], 2, 
                   colSums(M), FUN = "/") * 1e+06
        M <- log2(M+pc)
    }
    
    ## Initialize context and delta matrices and populate with the input data:
    if (is.null(MD)) {
        C <- matrix(log2(pc), nrow = length(treatm), ncol = nrow(genes), 
                    dimnames = list(colnames(M)[treatm], rownames(genes)), 
                    byrow = TRUE)
        D <- C * 0
        C[, idx.commonF] <- t(M[idx.commonR, cntr, drop = FALSE] )
        D[, idx.commonF] <- t(M[idx.commonR, treatm, drop = FALSE] ) - 
            t(M[idx.commonR, cntr, drop = FALSE])
    }
    
    if (!is.null(MD)) {
        C <- matrix(log2(pc), nrow = ncol(M), ncol = nrow(genes), 
                    dimnames = list(colnames(M), genes[, 1]), byrow = TRUE)
        D <- C * 0
        C[, idx.commonF] <-  t(M[idx.commonR, , drop = FALSE])
        D[, idx.commonF] <-  t(MD[idx.commonR, , drop = FALSE])
    }
    
    message("Encoding context...")
    LATC <- predict(encoder, list(gene_input = C))
    message("Encoding contrasts...")
    LATD <- predict(encoderD, list(delta_input = D, CONTEXT = LATC))
    message("Decoding contrasts...")
    DEC <- t(predict(generatorD, cbind(LATD, LATC)))
    
    message("Preparing output...")
    INPUT_DLT <- t(D)
    RES <- INPUT_DLT - DEC
    dimnames(DEC) <- dimnames(RES)
    
    VAR_DEC <- cbind(M2Mcor(INPUT_DLT, DEC)**2, M2Mcor(INPUT_DLT,RES)**2)
    VAR_DEC <- cbind(VAR_DEC, rowSums(VAR_DEC) - 1)
    colnames(VAR_DEC) <- c("DECODED","RESIDUAL", "COMMON")
    
    decomposed.contrasts <- list(INPUT_CONTRASTS = INPUT_DLT,
                                 DECODED_CONTRASTS = DEC,
                                 RESIDUAL_CONTRASTS = RES, 
                                 CONTEXT = t(C))
    RESULT <- SummarizedExperiment(assays = decomposed.contrasts,
                                   colData = list(ACCOUNTED_VARIANCE = VAR_DEC))
    rownames(RESULT) <- geneIDs
    message("Done!")
    return(RESULT)
}
