#### deJUNKER
Variance decomposition using conditional variational auto-encoders


## Function decompose.var:
decompose.var <- function(M,MD=NULL,treatm=NULL,cntr=NULL, process_input=TRUE, organism="Human",
                      feature_type=c("AUTO","ENSEMBL_GENE_ID", "GENE_SYMBOL", "ENTREZ_GENE_ID", "ARCHS4_ID")  )

# Decomposes input contrasts to decoded and residual fractions according to a trained contrast encoder-decoder.
# M is a matrix of raw gene counts. 
# treatm and cntr are vectors indicating column indices in M corresponding to treatments and controls. If treatm and cntr are specified MD has to be NULL.
# Alternatively (if MD is specified), M is assumed to be a raw gene count matrix specifying context for contrasts specified in MD
# MD is then a matrix of gene deltas with the same dimensions as M. If MD is specified treatm and cntr have to be NULL
# process_input: If set to TRUE (default) the count matrix will be preprocessed (library normalized, log-transformed, NA values will be set to 0).
# feature_type: Set to AUTO for automatic feature id-type detection. Alternatively specify the type of supplied id features.
# Returns a summarized experiment with the decomposed contrasts in the assays slots and
# the decomposed variance as colData
# Ex1: RES <- decompose.var( M=ContextM , MD=DeltaM, process_input = FALSE)
# Ex2: RES <- decompose.var( M=CountM , treatm=c(3,4,5,6), cntr=c(1,1,2,2), process_input = FALSE) 

## Function query.with.contrasts:
query.with.contrasts <- function( CONTRASTS=NULL, use=c("expressed.in.both", "all.genes"  ),
                                  expr.thr=0.25, organism="Human", preserve.in.GlobalEnv=TRUE,
                                  plot.contrast=c("RESIDUAL", "INPUT","DECODED","NONE"),
                                  detail.topn=10 )    
                                  
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
