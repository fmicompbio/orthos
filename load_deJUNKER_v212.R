reticulate::use_virtualenv("/tungstenfs/groups/gbioinfo/sharedSoft/virtualenvs/r-4.1-bioc-3.13-reticulate-keras-2.4.0-tensorflow-2.4.0-gpu/")
Sys.setenv("CUDA_VISIBLE_DEVICES" = "3" ) # Define visible  GPU devices
ngpus=length(strsplit(Sys.getenv("CUDA_VISIBLE_DEVICES"),",")[[1]])
reticulate::py_config()
library(keras)
K <- backend()
library(tensorflow)
library(dplyr)
library(Matrix)
library(HDF5Array)
library(SummarizedExperiment)
library(coop)
library(ggplot2)
library(ggrastr)
library(ggrepel)
library(ggpubr)
library(cowplot)
library(digest)



### Load model genes:
genes <- readRDS("/tungstenfs/groups/gbioinfo/papapana/DEEP_LEARNING/Autoencoders/ARCHS4/Rdata/ARCHS4_deJUNKER_v212_genes.rds")
ngenes <- nrow(genes)


### Correlation between corresponding columns of two matrices of dimensions n x k. Returns vector of length k
# M2Mcor <- function(M1,M2){
#     sapply(1:ncol(M1),function(i) cor(M1[,i,drop=FALSE],M2[,i,drop=FALSE]  ))
# }


### Decomposes input contrasts to decoded and residual fractions according to a trained contrast encoder-decoder.
### M is a matrix of raw gene counts.
### treatm and cntr are vectors indicating column indices in M corresponding to treatments and controls. If treatm and cntr are specified MD has to be NULL.
### Alternatively (if MD is specified), M is assumed to be a raw gene count matrix specifying context for contrasts specified in MD
### MD is then a matrix of gene deltas with the same dimensions as M. If MD is specified treatm and cntr have to be NULL
### process_input: If set to TRUE (default) the count matrix will be preprocessed (library normalized, log-transformed, NA values will be set to 0).
### feature_type: Set to AUTO for automatic feature id-type detection. Alternatively specify the type of supplied id features.
### Returns a summarized experiment with the decomposed contrasts in the assays slots and
### the decomposed variance as colData
### Ex1: RES <- decompose.var( M=ContextM , MD=DeltaM, process_input = FALSE)
### Ex2: RES <- decompose.var( M=CountM , treatm=c(3,4,5,6), cntr=c(1,1,2,2), process_input = FALSE)
decompose.var <- function(M,MD=NULL,treatm=NULL,cntr=NULL, process_input=TRUE, organism="Human",
                          feature_type=c("AUTO","ENSEMBL_GENE_ID", "GENE_SYMBOL", "ENTREZ_GENE_ID", "ARCHS4_ID")  )    {

    feature_type <- match.arg(feature_type)

    message("Checking input...")
    stopifnot( "`M` must be a (raw) count matrix with features (genes) in the rows and conditions in the columns"  = is.matrix(M) & all(M >= 0, na.rm = TRUE ) )
    stopifnot( "`specify either both `treatm` and `cntr` OR  `MD`"  = !is.null(MD) |  ( !is.null(treatm) & !is.null(cntr) )  )

    if (is.null(MD)) {
        stopifnot( "`treatm` and `cntr` vectors must have the same length"  = identical( length(treatm), length(cntr) ) )
        stopifnot( "`treatm` and `cntr` vectors must be numeric vectors indicating `M` column indices  " = all(c(treatm,cntr) %in% 1:ncol(M))    & !is.null(treatm) )
    }

    if (!is.null(MD)) {
        stopifnot( "`treatm` and `cntr` arguments must be NULL when `MD` is specified"  =  is.null(treatm) & is.null(cntr)   )
        stopifnot( "`M` and `MD` matrices have to have the same dimensions and dimnames"  = identical( dimnames(M), dimnames(MD) ) )
    }

    #### Feature ID type detection:
    geneIDs <- rownames(M)
    rownames(M) <- toupper(gsub("\\.\\d.+","",rownames(M) ))
    if (feature_type=="AUTO") {
        message("Detecting feature ids-type...")
        IDtypes <- c("ENSEMBL_GENE_ID", "GENE_SYMBOL", "ENTREZ_GENE_ID", "ARCHS4_ID")
        ID_OLaps <-  sapply(IDtypes, function(x) {  sum(rownames(M) %in% toupper(genes[,x])) } )
        feature_type <- IDtypes[which.max(ID_OLaps)]
        message("Feature ids-type detected: ", feature_type  )
    }
    else{
        message("Feature ids-type specified by user: ", feature_type  )
        ID_OLaps <-   sum(rownames(M) %in% toupper(genes[,feature_type]) )
    }

    message( max(ID_OLaps),"/", nrow(M) , " provided input features mapped against a total of ", nrow(genes), " model features. " )
    message( nrow(genes)-max(ID_OLaps)," missing features will be set to 0."  )
    message("--> Missing features corresponding to non/lowly expressed genes in your context(s) are of no consequence.")
    #message("--> A small fraction (<10%) of missing genes, expressed in your context(s), are well-tolerated.")
    message("--> The model is robust to small fractions (<10%) of missing genes that are expressed in your context(s).")
    message("--> Increased numbers of missing expressed genes in your input might result in model performance decline.")
    stopifnot("\n!!!Too many missing features (>10000).  Make certain you provide valid human
Symbol, Ensembl or Entrez gene identifiers as rownames in your input matrix `M`.\n" = (nrow(genes)-max(ID_OLaps) <= 10000) )
    idx.commonF <- na.omit(match( rownames(M), toupper(genes[,feature_type] ) ))
    idx.commonR <- which(rownames(M) %in% toupper(genes[,feature_type]))


    ########## Preprocess input (NAs ->0, LibNormalize, LogTransform):
    pc <- 4
    if (process_input){
        message("Preparing input...")
        #### Manage NA values:
        NA.idx <- which(is.na(M))
        if (length(NA.idx) >0){
            message("!!!Matrix `M` contains ",length(NA.idx), " NA values.", " Those will be set to 0"  )
            M[NA.idx] <- 0
        }
        ##### Lib normalize and log transform input count matrix:
        M <- sweep(M[rownames(M) %in% toupper(genes[,feature_type]), ], 2, colSums(M), FUN = "/") * 1e+06
        M <- log2(M+pc)
    }


    ## Initialize context and delta matrices and populate with the input data:
    if (is.null(MD)) {
        C <- matrix( log2(pc),nrow = length(treatm),ncol=nrow(genes), dimnames = list( colnames(M)[treatm], rownames(genes) ) , byrow = TRUE   )
        D <- C*0
        C[ , idx.commonF ] <-  t(M[idx.commonR, cntr, drop=FALSE] )
        D[ , idx.commonF ] <-  t(M[idx.commonR, treatm, drop=FALSE] ) - t(M[idx.commonR, cntr, drop=FALSE] )
    }

    if (!is.null(MD)) {
        C <- matrix( log2(pc),nrow =ncol(M),ncol=nrow(genes), dimnames = list( colnames(M), genes[,1]) , byrow = TRUE   )
        D <- C*0
        C[ , idx.commonF ] <-  t(M[idx.commonR, ,drop=FALSE] )
        D[ , idx.commonF ] <-  t(MD[idx.commonR, ,drop=FALSE] )
    }



    message("Encoding context...")
    LATC <- predict(encoder, list(gene_input=C))
    message("Encoding contrasts...")
    LATD <- predict(encoderD, list(delta_input=D , CONTEXT=LATC))
    message("Decoding contrasts...")
    DEC <- t( predict( generatorD,  cbind(LATD,LATC) )  )

    message("Preparing output...")
    INPUT_DLT <- t(D)
    RES <- INPUT_DLT - DEC
    dimnames(DEC) <- dimnames(RES)

    VAR_DEC <- cbind(M2Mcor(INPUT_DLT,DEC)**2,M2Mcor(INPUT_DLT,RES)**2 )
    VAR_DEC <- cbind(VAR_DEC, rowSums(VAR_DEC)-1 )
    colnames(VAR_DEC) <- c("DECODED","RESIDUAL","COMMON")

    decomposed.contrasts <- (list(INPUT_CONTRASTS=INPUT_DLT,DECODED_CONTRASTS=DEC,RESIDUAL_CONTRASTS=RES, CONTEXT=t(C) ) )
    RESULT <- SummarizedExperiment(assays = decomposed.contrasts,
                                   colData = list(ACCOUNTED_VARIANCE=VAR_DEC)
    )
    rownames(RESULT) <- geneIDs
    message("Done!")
    return(RESULT)
}







### Query the contrast database with a set of CONTRASTS
### Input is a summarized experiment with assays containing contrasts named
### INPUT_CONTRASTS, DECODED_CONTRASTS and RESIDUAL_CONTRASTS (at least one should be present)
### Typically the result of `decompose.var`.
### use determines if all.genes or genes expressed in both query and target context will be used.
### note that "expressed.in.both", though more accurate, is much slower.
### Expr.thr is the quantile in the provided context that determines the expression value above which a gene is considered to be expressed.
### This same value is then used for thresholding the contrast database. Only applies when use="expressed.in.both"
### preserve.in.GlobalEnv specifies whether the contrast database (stored as a SExp in `target.contrasts`)
### should be preserved in the GlobalEnv either for future queries or for accessing its metadata
### detail.topn specifies the number of top hits for which metadata will be returned in the TopHits slot of the results.
### Returns  a list of PearsonRhos, Zscores against the datbase as well as detailed Metadata for the detail.topn hits.
query.with.contrasts <- function( CONTRASTS=NULL, use=c("expressed.in.both", "all.genes"  ),
                                  expr.thr=0.25, organism="Human", preserve.in.GlobalEnv=TRUE,
                                  plot.contrast=c("RESIDUAL", "INPUT","DECODED","NONE"),
                                  detail.topn=10 )    {


    stopifnot( "`CONTRASTS` should be a valid SummarizedExperiment"= is(CONTRASTS, "SummarizedExperiment") )
    valid.contrasts <- c("INPUT_CONTRASTS", "DECODED_CONTRASTS", "RESIDUAL_CONTRASTS")
    present.contrasts <- intersect( valid.contrasts, names(assays(CONTRASTS)))
    stopifnot( "The assays slot in the provided SummarizedExperiment does not contain valid contrast names "=length(present.contrasts)>0 )
    message(paste ("provided contrast: ", present.contrasts, collapse="\n"))


    plot.contrast <-  match.arg(plot.contrast)


    smpl.col <- seq(1,75000,100) # A sample of row/col indexes to quickly check data integrity
    smpl.row <- seq(1,20000,100)
    if(!exists ("target.contrasts") ){
        message("Loading contrast database...")
        target.contrasts <- readRDS("/tungstenfs/groups/gbioinfo/papapana/DEEP_LEARNING/Autoencoders/ARCHS4/Rdata/DECOMPOSED_CONTRASTS_human_v212_uncompressed.rds")
        ### Probably better to move to HDF5 based implementation as this will be both a isgnificant speed-up in terms of loading
        ### and much more lean on memory requirements. However this needs a reworked implementation using DelayedMatrices.
        #   target.contrasts <- loadHDF5SummarizedExperiment(dir="/tungstenfs/groups/gbioinfo/papapana/DEEP_LEARNING/Autoencoders/ARCHS4/Rdata/DECOMPOSED_CONTRASTS_HDF5/",
        #                                               prefix="human_v212")
        if (preserve.in.GlobalEnv){
            assign("target.contrasts", target.contrasts, envir = .GlobalEnv)
        }
    }

    DBhash <- digest::digest(target.contrasts[smpl.row,smpl.col], algo="xxhash64")

    stopifnot("The contrast DB contained in the `target.contrasts` object has not been correctly loaded.
Please remove `target.contrasts` and try again."= DBhash=="f9abc421d4e93c08")

    stopifnot( "Incompatible rownames in the provided summarized experiment.
Rownames should be the same as in the contrast database.
You can make sure by generating your SE generated using `decompose.var`"
               = identical( rownames(CONTRASTS), rownames(target.contrasts) ) )
    CONTRASTS <- assays(CONTRASTS)[present.contrasts]


    if(use=="expressed.in.both")
    {
        #Set a global expression threshold according to a quantile in the query data context
        message("Thresholding genes...")
        thr <- quantile(CONTRASTS[["CONTEXT"]],expr.thr)
        CONTRASTS[["CONTEXT"]] [CONTRASTS[["CONTEXT"]] <= thr] <- NA
        assays(target.contrasts)[["CONTEXT"]][  assays(target.contrasts)[["CONTEXT"]] <= thr  ] <- NA

        pearson.rhos <- sapply( present.contrasts, function(x){
            message(paste0("Querying contrast database with ",x,"...")  )
            cor(CONTRASTS[[x]], assays(target.contrasts)[[x]], use = "pairwise.complete.obs" )  }, simplify = FALSE,USE.NAMES = TRUE
        )
    }


    if(use=="all.genes"){
        pearson.rhos <- sapply( present.contrasts, function(x){
            message(paste0("Querying contrast database with ",x,"...")  )
            cor(CONTRASTS[[x]], assays(target.contrasts)[[x]] )  }, simplify = FALSE,USE.NAMES = TRUE
        )
    }


    ### Z-score and pval calculation for each decomposed component query:
    message("Compiling query statistics...")
    zscores <- sapply(present.contrasts, function(x){
        t(scale(t(query.res$pearson.rhos[[x]])))
    },simplify = FALSE,USE.NAMES = TRUE
    )


    TopHits <- sapply(present.contrasts, function(contr){
        apply(zscores[[contr]],1,function(x) { N <- names(sort(x,decreasing=TRUE)[1:detail.topn])
        colData(target.contrasts)[N,c(3,12,22,24,29,31,33)]
        }    )
    },simplify = FALSE,USE.NAMES = TRUE
    )


    if(plot.contrast!="NONE" && paste0(plot.contrast,"_CONTRASTS") %in% present.contrasts){
        message("Generating plots...")
        plot.data <- zscores[[paste0(plot.contrast,"_CONTRASTS")]]
        PLOTLIST <- sapply( 1:nrow(plot.data),function(i){
            plot.query_results(plot.data[i,])
        },simplify=FALSE )
        suppressWarnings({
            print(
                cowplot::plot_grid(plotlist=PLOTLIST,labels=rownames(plot.data))
            )
        })
    }


    message("Done!")
    return(list(pearson.rhos=pearson.rhos, zscores=zscores,TopHits=TopHits )   )
}











##### Helper plotting function
plot.query_results <- function(scores, topn = 10) {
    min.score <- sort(scores, decreasing = TRUE)[topn]
    DF <-
        data.frame(
            idx = 1:length(scores),
            score = scores,
            ACC = names(scores)
        )
    blank.theme <-
        theme(
            axis.line = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.position = "none",
            panel.background = element_blank(),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.background = element_blank()
        )


    dens.plot <-
        ggpubr::ggdensity(DF, "score", fill = "#33638DFF") + clean_theme() +
        geom_point(
            data = DF %>% filter(score >= min.score),
            aes(x = score, y = 0),
            size = 1.5,
            color = "red"
        ) +
        ggpubr::rotate() +  theme(plot.margin = unit(c(1, 0, 1, 0), "cm")) +
        blank.theme

    manh.plot <-
        ggplot(data = DF, aes(x = idx, y = score)) + geom_bin2d(bins = 200) +
        scale_fill_continuous(type = "viridis") + theme_bw() +
        geom_point(
            data = subset(DF, score >= min.score),
            color = "red",
            size = 1.5
        ) +
        geom_text_repel(
            data = subset(DF, score >= min.score),
            aes(x = idx, y = score , label = ACC),
            size = 3
        ) +
        theme(legend.position = 'none') +  theme(plot.margin = unit(c(1, 0, 1, 0), "cm")) +
        ggpubr::rremove("x.ticks") + theme(panel.border = element_blank())

    P <-
        cowplot::plot_grid(
            manh.plot,
            NULL,
            dens.plot,
            ncol = 3,
            align = "hv",
            rel_widths = c(4, -0.2, 1)
        )
    return(P)
}







######## cVAE (conditional Variational Autoencoder for contrasts, conditioned on transcriptomic context)
######## Model definition:
if (tensorflow::tf$executing_eagerly())
    tensorflow::tf$compat$v1$disable_eager_execution()
if(!exists ("cvae")){

    # Parameters --------------------------------------------------------------
    latentD <- 512L # Latent dimension for Delta encoder
    latentC <- 64L  # Latent dimension for Context encoder
    drop_rate <- 0.2 # Dropout rate
    gene_dim <- ngenes  #Number of features (genes) in the dataset
    epsilon_std <- 0.5  ##Standard deviation of the prior latent distribution (vanilla =1)
    var_prior <- epsilon_std**2
    log_var_prior <- log(var_prior)
    kl_weight <- 0.20   #Weight for the Delta Kulllback-Leibler divergence loss (vanilla =1 )

    ###  Define encoder input layers:
    xD <- layer_input(shape = c(gene_dim),name="delta_input") #
    cont_input <- layer_input(shape=c(latentC),name='CONTEXT')
    all_inputs <- layer_concatenate(list(xD, cont_input  ))

    #### Delta encoder definition:
    hD <- layer_dense(all_inputs, 4 * latentD, activation = "elu")
    hD <- layer_dropout(hD, rate = drop_rate)
    hD <- layer_dense(hD, 2 * latentD, activation = "elu")
    hD <- layer_dropout(hD, rate = drop_rate)
    z_meanD <- layer_dense(hD, latentD)
    z_log_varD <- layer_dense(hD, latentD)

    # Define delta encoder
    encoderD <- keras_model(inputs=c(xD,cont_input), outputs=z_meanD)

    #### Sampling from the Delta latent space:
    samplingD <- function(arg){
        z_meanD <- arg[, 1:(latentD)]
        z_log_varD <- arg[, (latentD + 1):(2 * latentD)]
        epsilonD <- K$random_normal(
            shape = c(K$shape(z_meanD)[[1]]),
            mean=0.,
            stddev=epsilon_std
        )
        z_meanD + K$exp(z_log_varD/2)*epsilonD
    }

    # Lambda layer for variational sampling:
    zD <- layer_concatenate(list(z_meanD, z_log_varD)) %>%
        layer_lambda(samplingD)

    # Merge delta latent space with context latent space:
    z_concat<- layer_concatenate(list(zD,cont_input))

    # Define layers for the Delta decoder (no batch-norm. Seems to increase overfitting):
    decoder_h1 <-  layer_dense(units=2*latentD ,activation="elu")
    decoder_h2 <-  layer_dropout(rate = drop_rate)
    decoder_h3 <-  layer_dense( units=4*latentD ,activation="elu")
    decoder_h4 <-  layer_dropout(rate = drop_rate)

    decoder_out <- layer_dense(units = ngenes, activation = "linear")

    h_p <- decoder_h1(z_concat)
    h_p <- decoder_h2( h_p )
    h_p <- decoder_h3( h_p )
    h_p <- decoder_h4( h_p )
    outputs <- decoder_out(h_p)

    # Define full cvae (Input: lcpm, deltas  Target Outputs: Decoded deltas)
    cvae <- keras_model(inputs=c(xD,cont_input), outputs)

    # Reuse decoder layers to define the generator (concatenated latent space to gene output) separately
    d_in <- layer_input(shape = latentC+latentD)
    d_h <- decoder_h1(d_in)
    d_h <- decoder_h2(d_h)
    d_h <- decoder_h3(d_h)
    d_h <- decoder_h4(d_h)
    d_out <- decoder_out(d_h)
    generatorD <- keras_model(d_in, d_out)

    cvae_loss <- function(x, x_decoded_mean){
        reconstruction_loss  <-  loss_mean_squared_error(x, x_decoded_mean)
        kl_loss <- -kl_weight*0.5*K$mean(1 + z_log_varD-log_var_prior - K$square(z_meanD)/var_prior - K$exp(z_log_varD)/var_prior, axis = -1L)  # Delta KL loss
        reconstruction_loss + kl_loss
    }

    #### Custom correlation function to keep track of.
    #### All mathematical operation NEED to be performed using the Keras backend (e.g K$mean):
    cor_metric <- function(y_true, y_pred) {
        y_true_dev <- y_true - k_mean(y_true)
        y_pred_dev <- y_pred - k_mean(y_pred)
        r_num <- k_sum(y_true_dev * y_pred_dev)
        r_den <- k_sqrt(k_sum(k_square(y_true_dev)) *
                            k_sum(keras::k_square(y_pred_dev)))
        r_num / r_den
    }

    cvae %>% compile(
        loss = cvae_loss,
        optimizer = "adam",
        metrics = custom_metric("cor",cor_metric)
    )

    ### Load model weights:
    cvae %>% load_model_weights_hdf5("/tungstenfs/groups/gbioinfo/papapana/DEEP_LEARNING/Autoencoders/ARCHS4/Trained_models/cvae_deJUNKER_ARCHS4_v212_ftune_512_64_v1bb.hdf5")

}





######## Variational autoencoder model for context encoding:
######## Model definition:

if(!exists ("vae")){
    # Parameters --------------------------------------------------------------
    neck <- 64L # Latent (bottleneck) dimension 256//64
    drop_rate <- 0.1 #
    gene_dim <- ngenes  #Number of features (genes) in your dataset
    latent_dim <- neck
    epsilon_std <- 0.5  ##Standard deviation of the prior latent distribution (vanilla =1)
    var_prior <- epsilon_std**2
    log_var_prior <- log(var_prior)
    kl_weight <- 0.2   #Weight for the Kulllback-Leibler divergence loss (vanilla =1 )

    # Encoder definition  (using the functional API) :
    x <- layer_input(shape = c(gene_dim),name="gene_input") #
    h <- layer_dense(x, 8 * neck, activation = "elu")
    h <- layer_dropout(h, rate = drop_rate)
    h <- layer_dense(h, 4 * neck, activation = "elu")
    h <- layer_dropout(h, rate = drop_rate)

    z_mean <- layer_dense(h, latent_dim)
    z_log_var <- layer_dense(h, latent_dim)

    #### Sampling from the latent space:
    sampling <- function(arg){
        z_mean <- arg[, 1:(latent_dim)]
        z_log_var <- arg[, (latent_dim + 1):(2 * latent_dim)]
        epsilon <- K$random_normal(
            shape = c(K$shape(z_mean)[[1]]),
            mean=0.,
            stddev=epsilon_std
        )
        z_mean + K$exp(z_log_var/2)*epsilon
    }

    # Lambda layer for variational sampling:
    z <- layer_concatenate(list(z_mean, z_log_var)) %>%
        layer_lambda(sampling)


    # we instantiate the decoder separately so as to reuse it later
    decoder_h <- keras_model_sequential()
    decoder_h %>%
        layer_dense(units=4 * neck,activation="elu") %>%
        layer_dropout(rate = drop_rate) %>%
        layer_dense(8 * neck, activation = "elu") %>%
        layer_dropout(rate = drop_rate)

    decoder_mean <- layer_dense(units = gene_dim, activation = "relu")
    h_decoded <- decoder_h(z)
    x_decoded_mean <- decoder_mean(h_decoded)

    # end-to-end autoencoder (again notice the use of the functional API):
    vae <- keras_model(x, x_decoded_mean)

    # encoder, from inputs to latent space, also using the functional API:
    encoder <- keras_model(x, z_mean)

    # generator, from latent space to reconstructed inputs
    decoder_input <- layer_input(shape = latent_dim)
    h_decoded_2 <- decoder_h(decoder_input)
    x_decoded_mean_2 <- decoder_mean(h_decoded_2)
    generator <- keras_model(decoder_input, x_decoded_mean_2)

    vae_loss <- function(x, x_decoded_mean){
        reconstruction_loss  <-  loss_mean_squared_error(x, x_decoded_mean)
        kl_loss <- -kl_weight*0.5*K$mean(1 + z_log_var-log_var_prior - K$square(z_mean)/var_prior - K$exp(z_log_var)/var_prior, axis = -1L)  # More general formula
        reconstruction_loss + kl_loss
    }

    vae %>% compile(
        loss = vae_loss,
        optimizer = "adam",
        metrics = custom_metric("cor",cor_metric)
    )

    vae %>% load_model_weights_hdf5("Trained_models/vae_deJUNKER_lcpm_ARCHS_v212_64.hdf5")

}


