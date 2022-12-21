#' Calculate correlation between corresponding columns of two matrices
#'
#' @param M1,M2 Numeric matrices of dimensions n x k.
#'
#' @return Vector of length k with correlations between corresponding columns.
#'
#' @author Panagiotis Papasaikas, Charlotte Soneson
#'
#' @importFrom stats cor
M2Mcor <- function(M1, M2) {
    vapply(seq_len(ncol(M1)), function(i) {
        stats::cor(M1[, i, drop = FALSE], M2[, i, drop = FALSE])
    }, NA_real_)
}



#' Calculate correlation between a numeric matrix and a (large) delayed arrayed using grid access
#'
#' @param query Numeric matrix n x k.
#' @param hdf5 HDF5Matrix/DelayedMatrix n x l, where l is typically >> k
#' @param chunk_size column dimension for the grid used to read blocks from the HDF5 Matrix. Should be smaller than/equal to the chunk size used for storing the hdf5 Matrix
#' @param use character string specifying the handling of missing data. See R's standard correlation function cor.
#' @param workers Number of workers used for parallelization
#'
#' @return Correlation matrix n x k
#'
#' @author Panagiotis Papasaikas
#'
#' @importFrom stats cor
#' @importFrom HDF5Array HDF5Array
#' @importFrom DelayedArray colAutoGrid read_block
#' @importFrom BiocParallel bplapply MulticoreParam
#' 
grid_cor <- function(query, hdf5, use="pairwise.complete.obs", chunk_size=1000,
                     workers=16){
    full_dim <- dim(hdf5)
    full_grid <- DelayedArray::colAutoGrid(hdf5, ncol=min(chunk_size, ncol(hdf5))) #grid contains entire columns
    nblock <- length(full_grid) 
    res <- BiocParallel::bplapply(seq_len(nblock), function(b){
        ref_block <- DelayedArray::read_block(hdf5, full_grid[[b]])
        cor_res <- stats::cor(query, ref_block, use=use)
        return(cor_res)}, BPPARAM = BiocParallel::MulticoreParam(workers = workers))
    return( do.call(cbind, res) )
}