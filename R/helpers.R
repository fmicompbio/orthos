#' Calculate correlation between corresponding columns of two matrices
#'
#' @param M1,M2 Numeric matrices of dimensions n x k.
#'
#' @return Vector of length k with correlations between corresponding columns.
#'
#' @author Panagiotis Papasaikas, Charlotte Soneson
#'
#' @importFrom stats cor
.M2Mcor <- function(M1, M2) {
    vapply(seq_len(ncol(M1)), function(i) {
        stats::cor(M1[, i, drop = FALSE], M2[, i, drop = FALSE])
    }, NA_real_)
}



#' Calculate correlation between a numeric matrix and a (large) HDF5MAtrix/DelayedMatrix using grid access. 
#' Both Matrices can contain NAs. If thr is specified, substitution of values<=thr with NAs will be performed 
#' on the HDF5Matrix
#'
#' @param query Numeric matrix n x k.
#' @param hdf5 HDF5Matrix/DelayedMatrix n x l, where l is typically >> k
#' @param chunk_size column dimension for the grid used to read blocks from the HDF5 Matrix. Should be smaller than/equal to the ncol chunkdim used to write the data on disk.
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
.grid_cor_wNAs <- function(query, hdf5, chunk_size=500,
                     workers=16, thr=FALSE){
    full_dim <- dim(hdf5)
    full_grid <- DelayedArray::colAutoGrid(hdf5, ncol=min(chunk_size, ncol(hdf5))) #grid contains entire columns
    nblock <- length(full_grid) 
    res <- BiocParallel::bplapply(seq_len(nblock), function(b){
        ref_block <- DelayedArray::read_block(hdf5, full_grid[[b]])
        if(thr){
        ref_block[ref_block <= thr] <- NA
        }
        cor_res <- stats::cor(query, ref_block, use="pairwise.complete.obs")
        return(cor_res)}, BPPARAM = BiocParallel::MulticoreParam(workers = workers))
    return( do.call(cbind, res) )
}





.grid_cor_woNAs <- function(query, hdf5, chunk_size=500,
                          workers=16){
    full_dim <- dim(hdf5)
    full_grid <- DelayedArray::colAutoGrid(hdf5, ncol=min(chunk_size, ncol(hdf5))) #grid contains entire columns
    nblock <- length(full_grid) 
    res <- BiocParallel::bplapply(seq_len(nblock), function(b){
        ref_block <- DelayedArray::read_block(hdf5, full_grid[[b]])
        cor_res <- stats::cor(query, ref_block, use="everything" )
        return(cor_res)}, BPPARAM = BiocParallel::MulticoreParam(workers = workers))
    return( do.call(cbind, res) )
}