# Set (and unset) the RETICULATE_ENABLE_PYTHON_FINALIZER environment 
# variable to avoid errors related to pybind11 and GIL starting with 
# reticulate 1.42.0

# keras uses two environment variables controlling the support for legacy 
# models: TF_USE_LEGACY_KERAS and PY_REQUIRE_LEGACY_KERAS.
# If we run into issues loading the models, consider setting one (or both) of 
# these to "1".

.orthos_pkg_env <- new.env(parent = emptyenv())

.onLoad <- function(libname, pkgname) {
    .orthos_pkg_env$ORTHOS_RETICULATE_ENABLE_PYTHON_FINALIZER <- 
        Sys.getenv("RETICULATE_ENABLE_PYTHON_FINALIZER", unset = NA)
    Sys.setenv("RETICULATE_ENABLE_PYTHON_FINALIZER" = "yes")
}

.onUnload <- function(libpath) {
    if (!is.na(.orthos_pkg_env$ORTHOS_RETICULATE_ENABLE_PYTHON_FINALIZER)) {
        Sys.setenv("RETICULATE_ENABLE_PYTHON_FINALIZER" = 
                       .orthos_pkg_env$ORTHOS_RETICULATE_ENABLE_PYTHON_FINALIZER)
    } else {
        Sys.unsetenv("RETICULATE_ENABLE_PYTHON_FINALIZER")
    }
}
