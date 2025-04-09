# Set (and unset) the RETICULATE_ENABLE_PYTHON_FINALIZER environment 
# variable to avoid errors related to pybind11 and GIL starting with 
# reticulate 1.42.0

.onLoad <- function(libname, pkgname) {
    ORTHOS_RETICULATE_ENABLE_PYTHON_FINALIZER <- Sys.getenv("RETICULATE_ENABLE_PYTHON_FINALIZER")
    Sys.setenv("RETICULATE_ENABLE_PYTHON_FINALIZER" = "yes")
}

.onUnload <- function(libname, pkgname) {
    if (exists(ORTHOS_RETICULATE_ENABLE_PYTHON_FINALIZER)) {
        Sys.setenv("RETICULATE_ENABLE_PYTHON_FINALIZER" = ORTHOS_RETICULATE_ENABLE_PYTHON_FINALIZER)
    }
}
