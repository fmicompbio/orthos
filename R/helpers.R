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
