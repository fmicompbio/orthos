#'
#' @author Panagiotis Papasaikas, Charlotte Soneson
#'
#'
#' @importFrom basilisk BasiliskEnvironment
dejunkerenv <- basilisk::BasiliskEnvironment(
    envname = "dejunker", pkgname = "deJUNKER",
    packages = c("tensorflow-hub==0.12.0", "scipy==1.7.3",
                 "requests==2.28.1", "Pillow==9.2.0", "h5py==3.7.0",
                 "pandas==1.3.5", "tensorflow==2.10.0",
                 "keras==2.10.0", "pydot==1.4.2"),
    channels = c("bioconda", "conda-forge")
)
