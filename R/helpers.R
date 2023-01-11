#' Calculate correlation between corresponding columns of two matrices
#'
#' @param M1,M2 Numeric matrices of dimensions n x k.
#'
#' @return Vector of length k with correlations between corresponding columns.
#'
#' @author Panagiotis Papasaikas, Charlotte Soneson
#'
#' @importFrom stats cor
#' 
#' @keywords internal
.M2Mcor <- function(M1, M2) {
    .assertVector(x = M1, type = "matrix")
    .assertVector(x = M2, type = "matrix")
    vapply(seq_len(ncol(M1)), function(i) {
        stats::cor(M1[, i, drop = FALSE], M2[, i, drop = FALSE])
    }, NA_real_)
}



#' Calculate correlation between a numeric matrix and a (large)
#' HDF5MAtrix/DelayedMatrix using grid access. Both Matrices can contain NAs.
#' If \code{thr} is specified, substitution of values less than \code{thr} with
#' \code{NA}s will be performed on the HDF5Matrix
#'
#' @param query Numeric matrix n x k.
#' @param hdf5 HDF5Matrix/DelayedMatrix n x l, where l is typically >> k
#' @param chunk_size column dimension for the grid used to read blocks from the
#'   HDF5 Matrix. Should be larger than/equal to the ncol chunkdim used to write
#'   the data on disk.
#' @param workers Number of workers used for parallelization
#' @param thr If specified, a low bound on expression. Values lower than that
#'   are substituted by \code{NA}s on the HDF5 Matrix
#'
#' @return Correlation matrix k x l
#'
#' @author Panagiotis Papasaikas
#'
#' @importFrom stats cor
#' @importFrom HDF5Array HDF5Array
#' @importFrom DelayedArray colAutoGrid read_block
#' @importFrom BiocParallel bplapply MulticoreParam
#' 
#' @keywords internal
.grid_cor_wNAs <- function(query, hdf5, chunk_size = 1000,
                           workers = 16, thr = NULL){
    .assertVector(x = query, type = "matrix")
    .assertVector(x = hdf5, type = "HDF5Matrix")
    .assertScalar(x = chunk_size, type = "numeric") # add limit using `rngIncl`?
    .assertScalar(x = workers, type = "numeric", rngExcl = c(0, Inf))
    .assertScalar(x = thr, type = "numeric", allowNULL = TRUE)
    full_dim <- dim(hdf5)
    full_grid <- DelayedArray::colAutoGrid(hdf5, ncol = min(chunk_size, ncol(hdf5))) #grid contains entire columns
    nblock <- length(full_grid) 
    res <- BiocParallel::bplapply(seq_len(nblock), function(b){
        ref_block <- DelayedArray::read_block(hdf5, full_grid[[b]])
        if (!is.null(thr)) {
            ref_block[ref_block <= thr] <- NA
        }
        cor_res <- stats::cor(query, ref_block, use = "pairwise.complete.obs")
        return(cor_res)}, BPPARAM = BiocParallel::MulticoreParam(workers = workers))
    return( do.call(cbind, res) )
}


#' Calculate correlation between a numeric matrix and a (large)
#' HDF5MAtrix/DelayedMatrix using grid access. It is assumed that both Matrices
#' do not contain \code{NA}s. 
#'
#' @param query Numeric matrix n x k.
#' @param hdf5 HDF5Matrix/DelayedMatrix n x l, where l is typically >> k
#' @param chunk_size column dimension for the grid used to read blocks from the
#'   HDF5 Matrix. Should be larger than/equal to the ncol chunkdim used to write
#'   the data on disk.
#' @param workers Number of workers used for parallelization
#'
#' @return Correlation matrix k x l
#'
#' @author Panagiotis Papasaikas
#'
#' @importFrom stats cor
#' @importFrom HDF5Array HDF5Array
#' @importFrom DelayedArray colAutoGrid read_block
#' @importFrom BiocParallel bplapply MulticoreParam
#' 
#' @keywords internal
.grid_cor_woNAs <- function(query, hdf5, chunk_size = 1000,
                            workers = 16) {
    .assertVector(x = query, type = "matrix")
    .assertVector(x = hdf5, type = "HDF5Matrix")
    .assertScalar(x = chunk_size, type = "numeric") # add limit using `rngIncl`?
    .assertScalar(x = workers, type = "numeric", rngExcl = c(0, Inf))
    full_dim <- dim(hdf5)
    full_grid <- DelayedArray::colAutoGrid(hdf5, ncol = min(chunk_size, ncol(hdf5))) #grid contains entire columns
    nblock <- length(full_grid) 
    res <- BiocParallel::bplapply(seq_len(nblock), function(b){
        ref_block <- DelayedArray::read_block(hdf5, full_grid[[b]])
        cor_res <- stats::cor(query, ref_block, use = "everything" )
        return(cor_res)}, BPPARAM = BiocParallel::MulticoreParam(workers = workers))
    return( do.call(cbind, res) )
}

#' Utility function to check validity of scalar variable values.
#' 
#' This function provides a convenient way e.g. to check that provided
#' arguments to functions satisfy required criteria. 
#'  
#' @param x The variable to be checked.
#' @param type The desired type of \code{x}.
#' @param rngIncl The allowed range of the (numeric) variable \code{x}, 
#'     including the endpoints.
#' @param rngExcl The allowed range of the (numeric) variable \code{x}, 
#'     excluding the endpoints.
#' @param validValues A vector with the allowed values of \code{x}.
#' @param allowNULL Logical, whether or not \code{NULL} is an acceptable 
#'     value for \code{x}.
#' 
#' @author Michael Stadler, Charlotte Soneson
#' @noRd
#' @keywords internal
#' @importFrom methods is
.assertScalar <- function(x,
                          type = NULL,
                          rngIncl = NULL,
                          rngExcl = NULL,
                          validValues = NULL,
                          allowNULL = FALSE) {
    
    .assertVector(x = x, type = type, rngIncl = rngIncl,
                  rngExcl = rngExcl, validValues = validValues,
                  len = 1, rngLen = NULL, allowNULL = allowNULL)
    
}

#' Utility function to check validity of vector variable values.
#' 
#' This function provides a convenient way e.g. to check that provided
#' arguments to functions satisfy required criteria. 
#'  
#' @param x The variable to be checked
#' @param type The desired type of \code{x}.
#' @param rngIncl The allowed range of the (numeric) variable \code{x}, 
#'     including the endpoints.
#' @param rngExcl The allowed range of the (numeric) variable \code{x}, 
#'     excluding the endpoints.
#' @param validValues A vector with the allowed values of \code{x}.
#' @param len The required length of \code{x}.
#' @param rngLen The allowed range for the length of \code{x}.
#' @param allowNULL Logical, whether or not \code{NULL} is an acceptable 
#'     value for \code{x}.
#' 
#' @author Michael Stadler, Charlotte Soneson
#' @noRd
#' @keywords internal
#' @importFrom methods is
.assertVector <- function(x,
                          type = NULL,
                          rngIncl = NULL,
                          rngExcl = NULL,
                          validValues = NULL,
                          len = NULL, 
                          rngLen = NULL,
                          allowNULL = FALSE) {
    sc <- sys.calls()
    mycall <- sc[[length(sc)]]
    if (length(sc) >= 2 &&
        identical(as.character(sc[[length(sc) - 1]])[1], ".assertScalar")) {
        mycall <- sc[[length(sc) - 1]]
    }
    args <- lapply(mycall, as.character)[-1]
    xname <- if ("x" %in% names(args)) args$x else "argument"
    
    ## Check arguments
    stopifnot(is.null(type) || (length(type) == 1L && is.character(type)))
    stopifnot(is.null(rngIncl) || (length(rngIncl) == 2L && is.numeric(rngIncl)))
    stopifnot(is.null(rngExcl) || (length(rngExcl) == 2L && is.numeric(rngExcl)))
    stopifnot(is.null(len) || (length(len) == 1L && is.numeric(len)))
    stopifnot(is.null(rngLen) || (length(rngLen) == 2L && is.numeric(rngLen)))
    stopifnot(is.logical(allowNULL) && length(allowNULL) == 1L)
    if (!is.null(rngIncl) && !is.null(rngExcl)) {
        stop("'rngIncl' and 'rngExcl' can not both be specified")
    }
    
    ## If there are too many valid values, print only the first 15
    if (length(validValues) > 15) {
        vvPrint <- paste(c(validValues[seq_len(15)], 
                           "...(truncated)"),
                         collapse = ", ")
    } else {
        vvPrint <- paste(validValues, collapse = ", ")
    }
    
    if (is.null(x)) {
        if (allowNULL) {
            return(invisible(TRUE))
        } else {
            stop("'", xname, "' must not be NULL", call. = FALSE)
        }
    }
    
    if (is.null(type) && (!is.null(rngIncl) || !is.null(rngExcl))) {
        type <- "numeric"
    }
    
    if (!is.null(type) && !methods::is(x, type)) {
        stop("'", xname, "' must be of class '", type, "'", call. = FALSE)
    }
    
    if (!is.null(rngIncl)) {
        if (!is.null(validValues)) {
            if (any((x < rngIncl[1] | x > rngIncl[2]) & !(x %in% validValues))) {
                stop("'", xname, "' must be within [", rngIncl[1], ",", 
                     rngIncl[2], "] (inclusive), or one of: ", vvPrint,
                     call. = FALSE)
            }
        } else {
            if (any(x < rngIncl[1] | x > rngIncl[2])) {
                stop("'", xname, "' must be within [", rngIncl[1], ",", 
                     rngIncl[2], "] (inclusive)", call. = FALSE)
            }
        }
    } else if (!is.null(rngExcl)) {
        if (!is.null(validValues)) {
            if (any((x <= rngExcl[1] | x >= rngExcl[2]) & !(x %in% validValues))) {
                stop("'", xname, "' must be within (", rngExcl[1], ",", 
                     rngExcl[2], ") (exclusive), or one of: ", vvPrint,
                     call. = FALSE)
            }
        } else {
            if (any(x <= rngExcl[1] | x >= rngExcl[2])) {
                stop("'", xname, "' must be within (", rngExcl[1], ",", 
                     rngExcl[2], ") (exclusive)", call. = FALSE)
            }
        }
    } else {
        if (!is.null(validValues) && !all(x %in% validValues)) {
            stop("All values in '", xname, "' must be one of: ", vvPrint,
                 call. = FALSE)
        }
    }
    
    
    if (!is.null(len) && length(x) != len) {
        stop("'", xname, "' must have length ", len, call. = FALSE)
    }
    
    if (!is.null(rngLen) && (length(x) < rngLen[1] || length(x) > rngLen[2])) {
        stop("length of '", xname, "' must be within [", rngLen[1], ",", 
             rngLen[2], "] (inclusive)", call. = FALSE)
    }
    
    return(invisible(TRUE))
}

#' Utility function that makes sure that packages are available
#' 
#' The function tries loading the namespaces of the packages given in
#' \code{pkgs}, and throws an exception with an informative error message if
#' that is not the case.
#' 
#' @param pkgs Character vector with package names. Can be either just a
#'   package name or a string of the form \code{"githubuser/packagename"} for
#'   packages hosted on GitHub.
#' @param suggestInstallation Logical scalar. If \code{TRUE}, include an
#'   expression to install the missing package(s) as part of the generated
#'   error message.
#' 
#' @author Michael Stadler, Charlotte Soneson
#' 
#' @noRd
#' @keywords internal
.assertPackagesAvailable <- function(pkgs, suggestInstallation = TRUE) {
    stopifnot(exprs = {
        is.character(pkgs)
        is.logical(suggestInstallation)
        length(suggestInstallation) == 1L
    })
    
    avail <- unlist(lapply(sub("^[^/]+/", "", pkgs),
                           function(pkg) {
                               requireNamespace(pkg, quietly = TRUE)
                           }))
    
    if (any(!avail)) {
        caller <- deparse(sys.calls()[[sys.nframe() - 1]])
        callerfunc <- sub("\\(.+$", "", caller)
        haveBioc <- requireNamespace("BiocManager", quietly = TRUE)
        msg <- paste0("The package", ifelse(sum(!avail) > 1, "s '", " '"),
                      paste(sub("^[^/]+/", "", pkgs[!avail]), collapse = "', '"),
                      "' ",
                      ifelse(sum(!avail) > 1, "are", "is"), " required for ",
                      callerfunc, "(), but not installed.\n")
        if (suggestInstallation) {
            msg <- paste0(msg,
                          "Install ", ifelse(sum(!avail) > 1, "them", "it"), " using:\n",
                          ifelse(haveBioc, "", "install.packages(\"BiocManager\")\n"),
                          "BiocManager::install(c(\"",
                          paste(pkgs[!avail], collapse = "\", \""), "\"))")
        }
        stop(msg, call. = FALSE)
    }
    
    invisible(TRUE)
}
