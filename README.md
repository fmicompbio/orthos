# orthos
<img src="orthos.png" alt="orthos" width="200"/>

The measurement of effects on transcript abundance upon cellular perturbations is the most commonly used tool for probing molecular mechanisms of action due to its simplicity and low cost. However the interpretation of these experiments is hindered by the complex and nuanced effects of such perturbations on cellular processes. `orthos` is a generative modelling based approach that can disentangle specific/direct from generic/indirect effects of gene expression perturbations. Utilizing a large corpus of gene expression contrasts compiled from the [ARCHS4](https://maayanlab.cloud/archs4/)  database of uniformly processed RNAseq experiments and generative variational modelling we trained a model that accurately captures generic effects present in contrasts derived from RNAseq experiments while accounting for context. The model allows us to tease out specific perturbation components in the form of a contrast residual. The residual component obtained from this variance decomposition offers a more unequivocal signature of perturbations and is more closely related to the direct effects of a treatment when compared to the raw signal. `orthos` also provides a means to query the contrast database with arbitrary contrasts and identify related perturbations with similar specific effects providing among others the possibility to map treatments to mechanisms of action.


## Function `decompose.var`
`decompose.var` is a function that decomposes input contrasts to decoded and residual fractions according to a trained contrast encoder-decoder.

Input:
- `M`: a matrix of raw gene counts. If MD is specified, M is assumed to be a raw gene count matrix specifying context for the contrasts in `MD`.
If `MD` is not specified treatment and control conditions are specified using `treatm` and `cntr` respectively.
- `MD`: a matrix of gene deltas with the same dimensions as `M`
- `treatm`: vector indicating column indices in `M` corresponding to treatments. Only required if `MD` is not present.
- `cntr`: vector indicating column indices in `M` corresponding to controls. Only required if `MD` is not present.
- `process_input`: if `TRUE` (default), the count matrix will be preprocessed (library normalized, log-transformed, NA values will be set to 0). When 
set to `FALSE` it is assumed that `M` is a gene countr matrix that has undergone these steps of preprocessing.
- `organism`: organism of interest (default: Human)
-  `feature_type:` Set to AUTO for automatic feature id-type detection. Alternatively specify the type of supplied id features. Possible values are:
`ENSEMBL_GENE_ID`, `GENE_SYMBOL`, `ENTREZ_GENE_ID` and `ARCHS4_ID`.

Output: a summarized experiment object with the decomposed contrasts in the assays slots and the decomposed variance as colData.

Example 1: 

`RES <- decompose.var(M = ContextM, MD = DeltaM, process_input = FALSE)`

Example 2: 

`RES <- decompose.var(M = CountM, treatm = c(3,4,5,6), cntr = c(1,1,2,2), process_input = FALSE)`


## Function `query.with.contrasts`

`query.with.contrasts` is a function that queries the contrast database with a set of contrasts.

Input: 
- `CONTRASTS`: a summarized experiment with assays containing contrasts named INPUT_CONTRASTS, DECODED_CONTRASTS and RESIDUAL_CONTRASTS (at least one should be present). Typically the result of `decompose.var`.
- `use`: determines if all genes (`all.genes`) or genes expressed in both query and target context (`expressed.in.both`) will be used.
Note that the default `expressed.in.both` is more accurate but significantly slower.
- `expr.thr`: quantile in the provided context that determines the expression value above which a gene is considered to be expressed. 
 This same value is then used for thresholding the contrast database (only applies when `use` is set to `expressed.in.both`). 
- `organism`: organism of interest (default: Human)
- `preserve.in.GlobalEnv`: whether the contrast database (stored as a SExp in `target.contrasts`) should be preserved in the GlobalEnv either for future queries or for accessing its metadata
- `plot.contrast`: which of the contrast query results should be plotted. Should be one of `RESIDUAL`, `INPUT`,`DECODED` or `NONE`
- `detail.topn`: number of top hits for which metadata will be returned in the TopHits slot of the results

Output: a list of PearsonRhos, Zscores against the database as well as detailed Metadata for the top `detail.topn` hits.

Example: 

`RESULT <- query.with.contrasts(CONTRASTS=RES)`

## Installation

`orthos` can be installed from GitHub using:

```
devtools::install_github("fmicompbio/orthos")
```









