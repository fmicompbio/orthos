#' Preprocess count matrix
#'
#' The given count matrix \code{M} will be preprocessed by:
#' library normalization (calculation of counts per million, using the
#' column sums as library sizes), log2-transformation after addition of
#' \code{pseudocount}, and setting of \code{NA} values to 0.
#' The subset of \code{M} with row names in \code{ids} is returned.
#'
#' @param M Numeric matrix with feature counts (features are in rows,
#'     samples in columns).
#' @param ids Character vector. Only rows of \code{M} with row names in
#'     \code{ids} will be returned.
#' @param verbose Logical scalar. Report on progress if \code{TRUE}.
#' @param pseudocount Numeric scalar, added to the scaled (per million) counts
#'     before log2-transformation.
#'
#' @return Scaled, log2-transformed and subset version of \code{M}.
#'
#' @author Panagiotis Papasaikas, Michael Stadler
#'
#' @keywords internal
#' @noRd
.preprocessInput <- function(M, ids, verbose, pseudocount = 4) {
    if (verbose) {
        message("Preparing input...")
    }

    ## Manage NA values
    indicesNA <- which(is.na(M))

    if (length(indicesNA) > 0) {
        if (verbose) {
            message("Matrix `M` contains ", length(indicesNA),
                    " NA values. Those will be set to 0.")
        }
        M[indicesNA] <- 0
    }
    ## Lib normalize and log2 transform input count matrix
    M <- sweep(M[rownames(M) %in% ids, , drop = FALSE], 2,
               1e+06 / colSums(M), FUN = "*")
    M <- log2(M + pseudocount)
    return(M)
}


#' Load gene annotation table
#'
#' @param organism Character scalar, one of \code{"Human"} or \code{"Mouse"}.
#' @param mustWork Logical scalar. If \code{FALSE} and the gene information
#'     data is not available, return an \code{S4Vectors::DFrame} object with
#'     zero rows. If \code{TRUE} (the default) and the gene information data is
#'     not available, \code{.readGeneInformation} throws an error.
#'  @param verbose Logical scalar indicating whether to print messages along
#'     the way.
#'
#' @return \code{S4Vectors::DFrame} table with gene information.
#'
#' @author Panagiotis Papasaikas, Michael Stadler
#'
#' @importFrom SummarizedExperiment rowData
#' @importFrom S4Vectors DataFrame
#' @importFrom orthosData GetorthosContrastDB
#' 
#' @keywords internal
#' @noRd
.readGeneInformation <- function(organism, mustWork = TRUE, verbose = TRUE) {
    .assertScalar(x = organism, type = "character",
                  validValues = c("Human", "Mouse"))
    .assertScalar(x = mustWork, type = "logical")
    .assertScalar(x = verbose, type = "logical")
    
    geneInfoDir <- orthosData::GetorthosContrastDB(organism = organism,
                                                   mode = "DEMO", 
                                                   verbose = verbose)
    geneInfoFile <- paste0(tolower(organism), "_v212_NDF_c100_DEMOse.rds")
    geneInfoPath <- file.path(geneInfoDir, geneInfoFile)
    
    if (file.exists(geneInfoPath)) {
        genes <- SummarizedExperiment::rowData(
            readRDS(file.path(geneInfoDir, geneInfoFile))
        )
    } else {
        if (mustWork) {
            stop("gene information for '", organism, "' is not available")
        }
        genes <- S4Vectors::DataFrame(seqnames = factor(), start = integer(),
                                      end = integer(), width = integer(),
                                      strand = factor(),
                                      ENSEMBL_GENE_ID = character(),
                                      gene_name = character(),
                                      gene_biotype = character(),
                                      seq_coord_system = character(),
                                      description = character(),
                                      ENSEMBL_GENE_ID_VERSION = character(),
                                      ENSEMBL_CANONICAL_TRANSCRIPT = character(),
                                      GENE_SYMBOL = character(),
                                      ENTREZ_GENE_ID = character(),
                                      ARCHS4_ID = character())
    }
    return(genes)
}


#' Detect type of feature identifier
#'
#' @param featureType Character scalar, one of \code{"AUTO"},
#'     \code{"ENSEMBL_GENE_ID"}, \code{"GENE_SYMBOL"}, \code{"ENTREZ_GENE_ID"}
#'     or \code{"ARCHS4_ID"}.
#' @param genes Gene annotation data frame.
#' @param M Expression matrix (feature identifiers are in row names).
#' @param maxMissing Numeric scalar. This function will throw an error if more
#'     than \code{maxMissing} features in \code{genes} are not available in
#'     \code{M}.
#' @param verbose Logical scalar. If \code{TRUE}, report on progress.
#'
#' @return \code{featureType} as a character scalar.
#'
#' @author Panagiotis Papasaikas, Michael Stadler
#'
#' @keywords internal
#' @noRd
.detectFeatureIdType <- function(featureType, genes, M,
                                 maxMissing = 10000, verbose = TRUE) {
    if (featureType == "AUTO") {
        if (verbose) {
            message("Detecting feature ids-type...")
        }
        IDtypes <- c("ENSEMBL_GENE_ID", "GENE_SYMBOL", "ENTREZ_GENE_ID",
                     "ARCHS4_ID")
        IDoverlaps <- vapply(IDtypes, function(x) {
            sum(rownames(M) %in% toupper(genes[, x]))
        }, FUN.VALUE = numeric(1))
        featureType <- IDtypes[which.max(IDoverlaps)]
        if (verbose) {
            message("Feature ids-type detected: ", featureType)
        }
    } else {
        if (verbose) {
            message("Feature ids-type specified by user: ", featureType)
        }
        IDoverlaps <- sum(rownames(M) %in% toupper(genes[, featureType]))
    }

    if (verbose) {
        message(
            max(IDoverlaps), "/", nrow(M),
            " provided input features mapped against a total of ",
            nrow(genes), " model features.\n", nrow(genes) - max(IDoverlaps),
            " missing features will be set to 0.\n",
            "--> Missing features corresponding to non/lowly expressed genes ",
            "in your context(s) are of no consequence.\n",
            "--> The model is robust to small fractions (<10%) of missing ",
            "genes that are expressed in your context(s).\n",
            "--> Increased numbers of missing expressed genes in your input ",
            "might result in model performance decline.")
    }

    if (nrow(genes) - max(IDoverlaps) > maxMissing) {
        stop("\nToo many missing features (>", maxMissing, "). ",
             "Make certain you provided valid\nSymbol, Ensembl or ",
             "Entrez gene identifiers for the specified organism as ",
             "rownames in your input matrix `M`.\n")
    }

    return(featureType)
}


#' Decompose input contrasts to decoded and residual fractions
#'
#' Decompose input contrasts (gene expression deltas) to decoded (generic)
#' and residual (unique) components according to a contrast encoder-decoder
#' pre-trained on a large corpus of public RNAseq experiments.
#'
#' @author Panagiotis Papasaikas
#' @export
#'
#' @param M Matrix of raw gene counts.
#' @param MD Matrix of gene deltas (optional). If \code{MD} is specified,
#'     \code{M} is assumed to be a raw gene count matrix specifying context for
#'     contrasts specified in \code{MD}. \code{MD} is then a matrix of gene
#'     deltas with the same dimensions as \code{M}. If \code{MD} is specified,
#'     \code{treatm} and \code{cntr} have to be \code{NULL}.
#' @param treatm,cntr Vectors indicating column indices in \code{M}
#'     corresponding to treatments and controls. If \code{treatm} and
#'     \code{cntr} are specified, \code{MD} has to be \code{NULL}.
#' @param processInput If set to \code{TRUE} (default) the count matrix will
#'     be preprocessed (library normalized, log2-transformed after addition of
#'     a pseudocount, NA values will be set to 0).
#' @param organism Selects the autoencoder model trained on data from this
#'     species. One of \code{"Human"} or \code{"Mouse"}.
#' @param featureType Set to \code{"AUTO"} for automatic feature id-type
#'     detection. Alternatively specify the type of supplied id features.
#'     Current supported types are \code{"ENSEMBL_GENE_ID"},
#'     \code{"GENE_SYMBOL"}, \code{"ENTREZ_GENE_ID"} and \code{"ARCHS4_ID"}.
#' @param pseudocount Numerical scalar, added to raw counts in \code{M} when
#'     \code{preprocessInput = TRUE}.
#' @param verbose Logical scalar indicating whether to print messages along
#'     the way.
#'
#' @return A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'     with the decomposed contrasts in the assays and the decomposed variance
#'     as the \code{\link[SummarizedExperiment]{colData}}.
#'
#' @examples
#' \donttest{
#' 
#' MKL1_human <- readRDS(system.file("extdata", "GSE215150_MKL1_Human.rds",
#' package = "orthos"))
#' 
#' # Specifying M, treatm and cntr:
#' dec_MKL1_human <- decomposeVar(M = MKL1_human, treatm = c(2, 3), cntr = c(1, 1), 
#'                               organism = "Human", verbose = FALSE)
#'                               
#'                               
#' # Alternatively by specifying M and MD:
#' pseudocount <- 4 
#' M  <- sweep(MKL1_human, 2,
#'             colSums(MKL1_human), FUN = "/") * 1e+06
#' M  <- log2(M + pseudocount)
#' DeltaM <- M[,c("MKL1","caMKL1")]-M[,"Ctrl"] # Matrix of contrasts
#' ContextM <- M[,c("Ctrl","Ctrl")] # Matrix with context for the specified contrasts
#' colnames(ContextM) <- colnames(DeltaM) # M and MD need identical dimnames                       
#' RES <- decomposeVar(M = ContextM, MD = DeltaM, processInput = FALSE)
#' 
#' }
#'
#' @importFrom stats na.omit
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom basilisk basiliskRun basiliskStart basiliskStop
#'
decomposeVar <- function(M,
                         MD = NULL,
                         treatm = NULL, cntr = NULL,
                         processInput = TRUE,
                         organism = c("Human", "Mouse"),
                         featureType = c("AUTO", "ENSEMBL_GENE_ID",
                                         "GENE_SYMBOL", "ENTREZ_GENE_ID",
                                         "ARCHS4_ID"),
                         pseudocount = 4,
                         verbose = TRUE) {

    ## -------------------------------------------------------------------------
    ## Check input
    ## -------------------------------------------------------------------------
    if (verbose) {
        message("Checking input...")
    }
    .assertVector(x = M, type = "matrix", rngIncl = c(0, Inf))
    .assertVector(x = MD, type = "matrix", allowNULL = TRUE)
    .assertVector(x = treatm, type = "numeric", rngIncl = c(1, ncol(M)),
                  allowNULL = TRUE)
    .assertVector(x = cntr, type = "numeric", len = length(treatm),
                  rngIncl = c(1, ncol(M)), allowNULL = TRUE)
    .assertScalar(x = processInput, type = "logical")
    organism <- match.arg(organism)
    featureType <- match.arg(featureType)
    .assertScalar(x = pseudocount, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = verbose, type = "logical")
    stopifnot("`specify either `MD` OR both `treatm` and `cntr`" =
                  !is.null(MD) | (!is.null(treatm) | !is.null(cntr)))
    stopifnot("`rownames(M)` must be set to gene identifiers" =
                  !is.null(rownames(M)))
    stopifnot(
        "`M` and `MD` matrices have to have the same dimensions and dimnames" =
            is.null(MD) | (identical(dimnames(M), dimnames(MD)) &
                               identical(dim(M), dim(MD))))

    ## -------------------------------------------------------------------------
    ## Read gene information
    ## -------------------------------------------------------------------------
    genes <- .readGeneInformation(organism, verbose = verbose)
    ngenes <- nrow(genes)

    ## -------------------------------------------------------------------------
    ## Detect feature ID type
    ## -------------------------------------------------------------------------
    geneIDs <- rownames(M)
    rownames(M) <- toupper(gsub("\\.\\d+", "", rownames(M)))

    featureType <- .detectFeatureIdType(featureType = featureType,
                                        genes = genes,
                                        M = M,
                                        maxMissing = 10000,
                                        verbose = verbose)
    # Indices of input features in the model feature vector
    idx.commonF1 <- na.omit(match(rownames(M), toupper(genes[, featureType])))
    # Indices of model features in the input feature vector (i.e the rownames
    # of M)
    idx.commonR1 <- which(rownames(M) %in% toupper(genes[, featureType]))
    
    ## -------------------------------------------------------------------------
    ## Preprocess input (NAs -> 0, LibNormalize, LogTransform):
    ## -------------------------------------------------------------------------
    if (processInput) {
        M <- .preprocessInput(M = M,
                              ids = toupper(genes[, featureType]),
                              verbose = verbose,
                              pseudocount = pseudocount)
    }

    # Indices of input features in the model feature vector
    idx.commonF <- na.omit(match(rownames(M), toupper(genes[, featureType])))
    # Indices of model features in the input feature vector (i.e the rownames
    # of M)
    idx.commonR <- which(rownames(M) %in% toupper(genes[, featureType]))

    ## -------------------------------------------------------------------------
    ## Initialize context and delta matrices and populate with the input data:
    ## -------------------------------------------------------------------------
    if (is.null(MD)) {
        C <- matrix(log2(pseudocount), nrow = length(treatm),
                    ncol = nrow(genes),
                    dimnames = list(colnames(M)[treatm], rownames(genes)),
                    byrow = TRUE)
        D <- C * 0
        C[, idx.commonF] <- t(M[idx.commonR, cntr, drop = FALSE])
        D[, idx.commonF] <- t(M[idx.commonR, treatm, drop = FALSE]) -
            t(M[idx.commonR, cntr, drop = FALSE])
    } else {
        ## Manage NA values
        indicesNA <- which(is.na(MD))
        if (length(indicesNA) > 0) {
            if (verbose) {
                message("Matrix `MD` contains ", length(indicesNA),
                        " NA values. Those will be set to 0")
            }
            MD[indicesNA] <- 0
        }

        C <- matrix(log2(pseudocount), nrow = ncol(M), ncol = nrow(genes),
                    dimnames = list(colnames(M), genes[, 1]), byrow = TRUE)
        D <- C * 0
        C[, idx.commonF] <-  t(M[idx.commonR, , drop = FALSE])
        D[, idx.commonF] <-  t(MD[idx.commonR, , drop = FALSE])
    }

    if (verbose) {
        message("Encoding context...")
    }
    cl <- basiliskStart(orthosenv,
                        testload = "tensorflow")
    LATC <- basilisk::basiliskRun(proc = cl,
                                  fun = .predictEncoder,
                                  organism = organism,
                                  gene_input = C)
    if (verbose) {
        message("Encoding and decoding contrasts...")
    }
    res <- basilisk::basiliskRun(proc = cl,
                                 fun = .predictEncoderD,
                                 organism = organism,
                                 delta_input = D, context = LATC)
    basilisk::basiliskStop(cl)
    LATD <- res$LATD
    DEC <- res$DEC

    ## -------------------------------------------------------------------------
    ## Prepare output
    ## -------------------------------------------------------------------------
    if (verbose) {
        message("Preparing output...")
    }
    INPUT_DLT <- t(D)
    RES <- INPUT_DLT - DEC
    dimnames(DEC) <- dimnames(RES)

    VAR_DEC <- cbind(.M2Mcor(INPUT_DLT, DEC)**2, .M2Mcor(INPUT_DLT, RES)**2)
    VAR_DEC <- cbind(VAR_DEC, rowSums(VAR_DEC) - 1)
    colnames(VAR_DEC) <- c("DECODED", "RESIDUAL", "COMMON")

    decomposed.contrasts <- list(INPUT_CONTRASTS = INPUT_DLT,
                                 DECODED_CONTRASTS = DEC,
                                 RESIDUAL_CONTRASTS = RES,
                                 CONTEXT = t(C))
    RESULT <- SummarizedExperiment(assays = decomposed.contrasts,
                                   colData = list(ACCOUNTED_VARIANCE = VAR_DEC))
    rownames(RESULT) <- rownames(genes)
    rowData(RESULT)$User_provided_IDs <- rep(NA, nrow(RESULT))
    rowData(RESULT)$User_provided_IDs[idx.commonF1] <- geneIDs[idx.commonR1]

    if (verbose) {
        message("Done!")
    }

    return(RESULT)
}

#' @keywords internal
#' @noRd
#'
#' @importFrom keras load_model_hdf5
#' @importFrom stats predict
#' @importFrom ExperimentHub ExperimentHub
#' @importFrom AnnotationHub query
#'
.predictEncoder <- function(gene_input, organism) {
    ## Load context encoder model from ExperimentHub:
    query_keys <- c("orthosData", "ContextEncoder_", organism, "ARCHS4")
    hub <- ExperimentHub::ExperimentHub()
    encoder <- AnnotationHub::query(hub, query_keys)[[1]]
    predict(encoder, list(gene_input = gene_input))
}

#' @keywords internal
#' @noRd
#'
#' @importFrom keras load_model_hdf5
#' @importFrom stats predict
#' @importFrom ExperimentHub ExperimentHub
#' @importFrom AnnotationHub query
#'
.predictEncoderD <- function(delta_input, context, organism) {
    ## Load contrast encoder and generator models from ExperimentHub:
    query_keysE <- c("orthosData", "DeltaEncoder_", organism, "ARCHS4")
    query_keysD <- c("orthosData", "DeltaDecoder_", organism, "ARCHS4")
    
    hub <- ExperimentHub::ExperimentHub()
    encoderD <- AnnotationHub::query(hub, query_keysE)[[1]]
    generatorD <- AnnotationHub::query(hub, query_keysD)[[1]]
    ## Encode and decode deltas:
    LATD <- predict(encoderD, list(delta_input = delta_input,
                                   CONTEXT = context))
    DEC <- t(predict(generatorD, cbind(LATD, context)))
    
    list(LATD = LATD, DEC = DEC)
}


