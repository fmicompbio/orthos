#' Decompose input contrasts to decoded and residual fractions
#' 
#' Decompose input contrasts to decode and residual fractions according to a 
#' trained contrast encoder-decoder.
#' 
#' @author Panagiotis Papasaikas
#' @export
#' 
#' @param M Matrix of raw gene counts.
#' @param MD a vector of gene deltas. If MD is specified, M is assumed to be a 
#'     raw gene count matrix specifying context for contrasts specified in MD.
#'     MD is then a matrix of gene deltas with the same dimensions as M. 
#'     If MD is specified, \code{treatm} and \code{cntr} have to be \code{NULL}.
#' @param treatm,cntr Vectors indicating column indices in M corresponding 
#'     to treatments and controls. If treatm and cntr are specified MD has to be 
#'     \code{NULL}. 
#' @param processInput If set to \code{TRUE} (default) the count matrix will 
#'     be preprocessed (library normalized, log-transformed, NA values will 
#'     be set to 0).
#' @param featureType Set to AUTO for automatic feature id-type detection. 
#'     Alternatively specify the type of supplied id features.
#' @param verbose Logical scalar indicating whether to print messages along 
#'     the way.
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
                                         "ARCHS4_ID"), 
                         verbose = TRUE) {
    
    featureType <- match.arg(featureType)
    
    ## -------------------------------------------------------------------------
    ## Read gene information
    ## -------------------------------------------------------------------------
    genes <- readRDS(system.file("extdata", "ARCHS4_deJUNKER_v212_genes.rds", 
                                 package = "deJUNKER"))
    ngenes <- nrow(genes)
    
    ## -------------------------------------------------------------------------
    ## Check input
    ## -------------------------------------------------------------------------
    if (verbose) {
        message("Checking input...")
    }
    stopifnot("`M` must be a (raw) count matrix with features (genes) in the rows and conditions in the columns" = 
                  is.matrix(M) & all(M >= 0, na.rm = TRUE))
    stopifnot("`specify either both `treatm` and `cntr` OR  `MD`" = 
                  !is.null(MD) | (!is.null(treatm) & !is.null(cntr)))
    
    if (is.null(MD)) {
        stopifnot("`treatm` and `cntr` vectors must have the same length" = 
                      identical(length(treatm), length(cntr)))
        stopifnot("`treatm` and `cntr` vectors must be numeric vectors indicating `M` column indices " = 
                      all(c(treatm, cntr) %in% 1:ncol(M)) & !is.null(treatm))
    }
    
    if (!is.null(MD)) {
        stopifnot("`treatm` and `cntr` arguments must be NULL when `MD` is specified" = 
                      is.null(treatm) & is.null(cntr))
        stopifnot("`M` and `MD` matrices have to have the same dimensions and dimnames" = 
                      identical(dimnames(M), dimnames(MD)))
    }
    
    ## -------------------------------------------------------------------------
    ## Detect feature ID type
    ## -------------------------------------------------------------------------
    geneIDs <- rownames(M)
    rownames(M) <- toupper(gsub("\\.\\d.+", "", rownames(M)))
    if (featureType == "AUTO") {
        if (verbose) {
            message("Detecting feature ids-type...")
        }
        IDtypes <- c("ENSEMBL_GENE_ID", "GENE_SYMBOL", "ENTREZ_GENE_ID", "ARCHS4_ID")
        ID_OLaps <- sapply(IDtypes, function(x) {  
            sum(rownames(M) %in% toupper(genes[, x])) 
        })
        featureType <- IDtypes[which.max(ID_OLaps)]
        if (verbose) {
            message("Feature ids-type detected: ", featureType)
        }
    } else {
        if (verbose) {
            message("Feature ids-type specified by user: ", featureType)
        }
        ID_OLaps <- sum(rownames(M) %in% toupper(genes[, featureType]))
    }
    
    if (verbose) {
        message(max(ID_OLaps), "/", nrow(M), 
                " provided input features mapped against a total of ", 
                nrow(genes), " model features.")
        message(nrow(genes) - max(ID_OLaps), " missing features will be set to 0.")
        message("--> Missing features corresponding to non/lowly expressed genes ", 
                "in your context(s) are of no consequence.")
        message("--> The model is robust to small fractions (<10%) of missing ", 
                "genes that are expressed in your context(s).")
        message("--> Increased numbers of missing expressed genes in your input ", 
                "might result in model performance decline.")
    }
    stopifnot("\n!!!Too many missing features (>10000).  Make certain you provided valid 
Symbol, Ensembl or Entrez gene identifiers  for the specified organism as rownames in your input matrix `M`.\n" = 
                  (nrow(genes) - max(ID_OLaps) <= 10000))
    idx.commonF <- na.omit(match(rownames(M), toupper(genes[, featureType])))
    idx.commonR <- which(rownames(M) %in% toupper(genes[, featureType]))
    
    ## -------------------------------------------------------------------------
    ## Preprocess input (NAs ->0, LibNormalize, LogTransform):
    ## -------------------------------------------------------------------------
    pc <- 4
    if (processInput) {
        if (verbose) {
            message("Preparing input...")
        }
        
        #### Manage NA values:
        NA.idx <- which(is.na(M))
        if (length(NA.idx) > 0) {
            if (verbose) {
                message("!!!Matrix `M` contains ", length(NA.idx),
                        " NA values.", " Those will be set to 0")
            }
            M[NA.idx] <- 0
        }
        ##### Lib normalize and log transform input count matrix:
        M <- sweep(M[rownames(M) %in% toupper(genes[, featureType]), ], 2, 
                   colSums(M), FUN = "/") * 1e+06
        M <- log2(M + pc)
    }
    
    ## -------------------------------------------------------------------------
    ## Initialize context and delta matrices and populate with the input data:
    ## -------------------------------------------------------------------------
    if (is.null(MD)) {
        C <- matrix(log2(pc), nrow = length(treatm), ncol = nrow(genes), 
                    dimnames = list(colnames(M)[treatm], rownames(genes)), 
                    byrow = TRUE)
        D <- C * 0
        C[, idx.commonF] <- t(M[idx.commonR, cntr, drop = FALSE])
        D[, idx.commonF] <- t(M[idx.commonR, treatm, drop = FALSE]) - 
            t(M[idx.commonR, cntr, drop = FALSE])
    }
    
    if (!is.null(MD)) {
        C <- matrix(log2(pc), nrow = ncol(M), ncol = nrow(genes), 
                    dimnames = list(colnames(M), genes[, 1]), byrow = TRUE)
        D <- C * 0
        C[, idx.commonF] <-  t(M[idx.commonR, , drop = FALSE])
        D[, idx.commonF] <-  t(MD[idx.commonR, , drop = FALSE])
    }
    
    if (verbose) {
        message("Encoding context...")
    }
    LATC <- basiliskRun(env = dejunkerenv, fun = .predict_encoder, 
                        gene_input = C )
    if (verbose) {
        message("Encoding and decoding contrasts...")
    }
    res <- basiliskRun(env = dejunkerenv, fun = .predict_encoderd, 
                       delta_input = D, context = LATC, ngenes = ngenes)
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
    
    VAR_DEC <- cbind(M2Mcor(INPUT_DLT, DEC)**2, M2Mcor(INPUT_DLT, RES)**2)
    VAR_DEC <- cbind(VAR_DEC, rowSums(VAR_DEC) - 1)
    colnames(VAR_DEC) <- c("DECODED", "RESIDUAL", "COMMON")
    
    decomposed.contrasts <- list(INPUT_CONTRASTS = INPUT_DLT,
                                 DECODED_CONTRASTS = DEC,
                                 RESIDUAL_CONTRASTS = RES, 
                                 CONTEXT = t(C))
    RESULT <- SummarizedExperiment(assays = decomposed.contrasts,
                                   colData = list(ACCOUNTED_VARIANCE = VAR_DEC))
    rownames(RESULT) <- geneIDs
    if (verbose) {
        message("Done!")
    }
    
    return(RESULT)
}

#' @keywords internal
#' @noRd
#' 
#' @importFrom keras backend layer_input layer_dense layer_dropout
#'     layer_concatenate layer_lambda keras_model_sequential keras_model
#'     compile custom_metric loss_mean_squared_error
#'     
.predict_encoder <- function(gene_input) {
    encoder <- keras::load_model_hdf5("/tungstenfs/groups/gbioinfo/papapana/DEEP_LEARNING/Autoencoders/ARCHS4/Trained_models/Model_encoder_deJUNKER_lcpm_ARCHS_v212_human.hdf5",compile=FALSE)
    predict(encoder, list(gene_input = gene_input))
}

#' @keywords internal
#' @noRd
#' 
#' @importFrom keras backend layer_input layer_concatenate layer_dense
#'     layer_dropout keras_model layer_lambda custom_metric k_mean k_sum 
#'     k_square k_sqrt loss_mean_squared_error
#' @importFrom tensorflow tf
#' 
.predict_encoderd <- function(delta_input, context, ngenes) {
    K <- keras::backend()
    if (tensorflow::tf$executing_eagerly())
        tensorflow::tf$compat$v1$disable_eager_execution()
    
    # Parameters --------------------------------------------------------------
    latentD <- 512L # Latent dimension for Delta encoder
    latentC <- 64L # Latent dimension for Context encoder
    drop_rate <- 0.2 # Dropout rate
    gene_dim <- ngenes # Number of features (genes) in the dataset
    epsilon_std <- 0.5 # Standard deviation of the prior latent distribution (vanilla = 1)
    var_prior <- epsilon_std**2
    log_var_prior <- log(var_prior)
    kl_weight <- 0.20 # Weight for the Delta Kulllback-Leibler divergence loss (vanilla = 1)
    
    ### Define encoder input layers:
    xD <- keras::layer_input(shape = c(gene_dim), name = "delta_input")
    cont_input <- keras::layer_input(shape = c(latentC), name = 'CONTEXT')
    all_inputs <- keras::layer_concatenate(list(xD, cont_input))
    
    #### Delta encoder definition:
    hD <- keras::layer_dense(all_inputs, 4 * latentD, activation = "elu")
    hD <- keras::layer_dropout(hD, rate = drop_rate)
    hD <- keras::layer_dense(hD, 2 * latentD, activation = "elu")
    hD <- keras::layer_dropout(hD, rate = drop_rate)
    z_meanD <- keras::layer_dense(hD, latentD)
    z_log_varD <- keras::layer_dense(hD, latentD)
    
    # Define delta encoder
    encoderD <- keras::keras_model(inputs = c(xD, cont_input), outputs = z_meanD)
    
    #### Sampling from the Delta latent space:
    samplingD <- function(arg) {
        z_meanD <- arg[, seq_len(latentD)]
        z_log_varD <- arg[, (latentD + 1):(2 * latentD)]
        epsilonD <- K$random_normal(
            shape = c(K$shape(z_meanD)[[1]]),
            mean = 0.,
            stddev = epsilon_std
        )
        z_meanD + K$exp(z_log_varD/2) * epsilonD
    }
    
    # Lambda layer for variational sampling:
    zD <- keras::layer_concatenate(list(z_meanD, z_log_varD)) %>%
        keras::layer_lambda(samplingD)
    
    # Merge delta latent space with context latent space:
    z_concat <- keras::layer_concatenate(list(zD, cont_input))
    
    # Define layers for the Delta decoder (no batch-norm. Seems to increase overfitting):
    decoder_h1 <- keras::layer_dense(units = 2*latentD, activation = "elu")
    decoder_h2 <- keras::layer_dropout(rate = drop_rate)
    decoder_h3 <- keras::layer_dense(units = 4*latentD, activation = "elu")
    decoder_h4 <- keras::layer_dropout(rate = drop_rate)
    
    decoder_out <- keras::layer_dense(units = ngenes, activation = "linear")
    
    h_p <- decoder_h1(z_concat)
    h_p <- decoder_h2(h_p)
    h_p <- decoder_h3(h_p)
    h_p <- decoder_h4(h_p)
    outputs <- decoder_out(h_p)
    
    # Define full cvae (Input: lcpm, deltas  Target Outputs: Decoded deltas)
    cvae <- keras::keras_model(inputs = c(xD, cont_input), outputs)
    
    # Reuse decoder layers to define the generator (concatenated latent space to gene output) separately
    d_in <- layer_input(shape = latentC + latentD)
    d_h <- decoder_h1(d_in)
    d_h <- decoder_h2(d_h)
    d_h <- decoder_h3(d_h)
    d_h <- decoder_h4(d_h)
    d_out <- decoder_out(d_h)
    generatorD <- keras::keras_model(d_in, d_out)
    
    cvae_loss <- function(x, x_decoded_mean) {
        reconstruction_loss <- keras::loss_mean_squared_error(x, x_decoded_mean)
        kl_loss <- -kl_weight * 0.5 * 
            K$mean(1 + z_log_varD - log_var_prior - K$square(z_meanD) / var_prior - 
                       K$exp(z_log_varD) / var_prior, axis = -1L)  # Delta KL loss
        reconstruction_loss + kl_loss
    }
    
    #### Custom correlation function to keep track of.
    #### All mathematical operation NEED to be performed using the Keras backend (e.g K$mean):
    cor_metric <- function(y_true, y_pred) {
        y_true_dev <- y_true - keras::k_mean(y_true)
        y_pred_dev <- y_pred - keras::k_mean(y_pred)
        r_num <- keras::k_sum(y_true_dev * y_pred_dev)
        r_den <- keras::k_sqrt(k_sum(keras::k_square(y_true_dev)) *
                                   keras::k_sum(keras::k_square(y_pred_dev)))
        r_num / r_den
    }
    
    cvae %>% keras::compile(
        loss = cvae_loss,
        optimizer = "adam",
        metrics = keras::custom_metric("cor", cor_metric)
    )
    
    ### Load model weights:
    # cvae %>% keras::load_model_weights_hdf5("/Users/charlottesoneson/Documents/Rpackages/deJUNKER/data/cvae_deJUNKER_ARCHS4_v212_ftune_512_64_v1bb.hdf5")
    cvae %>%
        keras::load_model_weights_hdf5("/tungstenfs/groups/gbioinfo/papapana/DEEP_LEARNING/Autoencoders/ARCHS4/Trained_models/cvae_deJUNKER_ARCHS4_v212_ftune_512_64_v1bb.hdf5")
    
    LATD <- predict(encoderD, list(delta_input = delta_input, CONTEXT = context))
    DEC <- t(predict(generatorD, cbind(LATD, context)))
    
    list(LATD = LATD, DEC = DEC)
}