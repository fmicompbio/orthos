#' @importFrom basilisk BasiliskEnvironment
dejunkerenv <- basilisk::BasiliskEnvironment(
    envname = "dejunker", pkgname = "deJUNKER",
    packages = c("tensorflow-hub", "scipy==1.7.3",
                 "requests", "Pillow", "h5py", "pydot",
                 "pandas==1.3.5", "tensorflow==2.11.0",
                 "keras==2.11.0"),
    channels = c("bioconda", "conda-forge"),
    pip = c("")
)
