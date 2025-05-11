## List of packages exported from working environment
## - pinned versions were obtained from
## https://github.com/csoneson/orthos_package_list

if (basilisk::isLinuxAarch64()) {
    .orthos_dependencies <- c(
        "python==3.10.14",
        "numpy==1.26.4",
        "keras==2.15.0",
        "tensorflow==2.15.0",
        "pandas==2.2.2"
    )
    .pip_dependencies <- character(0)
} else if (basilisk::isLinux()) {
    .orthos_dependencies <- c(
        "python==3.11.11",
        "numpy==1.26.4",
        "keras==2.15.0",
        "tensorflow==2.15.0",
        "pandas==2.2.3"
    )
    .pip_dependencies <- character(0)
} else if (basilisk::isMacOSXArm()) {
    .orthos_dependencies <- c(
        "python==3.11.11",
        "numpy==1.26.4",
        "keras==2.15.0",
        "tensorflow==2.15.0",
        "pandas==2.2.3"
    )
    .pip_dependencies <- character(0)
} else if (basilisk::isMacOSX()) {
    .orthos_dependencies <- c(
        "python==3.11.11",
        "numpy==1.26.4",
        "keras==2.15.0",
        "tensorflow==2.15.0",
        "pandas==2.2.3"
    )
    .pip_dependencies <- character(0)
} else if (basilisk::isWindows()) {
    .orthos_dependencies <- c(
        "python==3.8.19",
        "numpy==1.22.4",
        "keras==2.10.0",
        "tensorflow==2.10.0",
        "pandas==1.3.5"
    )
    .pip_dependencies <- character(0)
}

#' @author Charlotte Soneson
#'
#' @importFrom basilisk BasiliskEnvironment
orthosenv <- basilisk::BasiliskEnvironment(
    envname = "orthos", pkgname = "orthos",
    packages = .orthos_dependencies,
    pip = .pip_dependencies,
    channels = "conda-forge"
)
