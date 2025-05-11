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
        "absl-py==2.2.2", 
        "astunparse==1.6.3", 
        "cachetools==5.5.2", 
        "certifi==2025.4.26",
        "charset-normalizer==3.4.2", 
        "flatbuffers==25.2.10", 
        "gast==0.6.0",
        "google-auth==2.40.1",
        "google-auth-oauthlib==1.2.2", 
        "google-pasta==0.2.0",
        "grpcio==1.71.0", 
        "h5py==3.13.0", 
        "idna==3.10", 
        "keras==2.15.0",
        "libclang==18.1.1", 
        "Markdown==3.8",
        "MarkupSafe==3.0.2", 
        "ml-dtypes==0.2.0", 
        "numpy==1.26.4",
        "oauthlib==3.2.2",
        "opt_einsum==3.4.0", 
        "packaging==25.0", 
        "pandas==2.2.3",
        "pip==25.1.1", 
        "protobuf==4.25.7", 
        "pyasn1==0.6.1",
        "pyasn1_modules==0.4.2",
        "python==3.11.11",
        "python-dateutil==2.9.0.post0", 
        "pytz==2025.2", 
        "requests==2.32.3", 
        "requests-oauthlib==2.0.0", 
        "rsa==4.9.1",
        "setuptools==80.4.0",
        "six==1.17.0", 
        "tensorboard==2.15.2", 
        "tensorboard-data-server==0.7.2",
        "tensorflow==2.15.0", 
        "tensorflow-estimator==2.15.0", 
        "tensorflow-io-gcs-filesystem==0.37.1", 
        "termcolor==3.1.0", 
        "typing_extensions==4.13.2", 
        "tzdata==2025.2", 
        "urllib3==2.4.0",
        "Werkzeug==3.1.3", 
        "wheel==0.45.1", 
        "wrapt==1.14.1"
    )
    .pip_dependencies <- character(0)
} else if (basilisk::isMacOSXArm()) {
    .orthos_dependencies <- c(
        "python==3.11.11",
        "absl-py==2.2.2",
        "astunparse==1.6.3", 
        "cachetools==5.5.2", 
        "certifi==2025.4.26", 
        "charset-normalizer==3.4.2",
        "flatbuffers==25.2.10", 
        "gast==0.6.0", 
        "google-auth==2.40.1", 
        "google-auth-oauthlib==1.2.2", 
        "google-pasta==0.2.0", 
        "grpcio==1.71.0", 
        "h5py==3.13.0", 
        "idna==3.10", 
        "keras==2.15.0", 
        "libclang==18.1.1",
        "Markdown==3.8", 
        "MarkupSafe==3.0.2", 
        "ml-dtypes==0.2.0",
        "numpy==1.26.4", 
        "oauthlib==3.2.2",
        "opt_einsum==3.4.0", 
        "packaging==25.0",
        "pandas==2.2.3", 
        "pip==25.1.1", 
        "protobuf==4.25.7", 
        "pyasn1==0.6.1", 
        "pyasn1_modules==0.4.2", 
        "python-dateutil==2.9.0.post0", 
        "pytz==2025.2", 
        "requests==2.32.3", 
        "requests-oauthlib==2.0.0", 
        "rsa==4.9.1",
        "setuptools==80.4.0", 
        "six==1.17.0", 
        "tensorboard==2.15.2", 
        "tensorboard-data-server==0.7.2",
        "tensorflow==2.15.0", 
        "tensorflow-estimator==2.15.0", 
        "tensorflow-io-gcs-filesystem==0.37.1", 
        "tensorflow-macos==2.15.0", 
        "termcolor==3.1.0", 
        "typing_extensions==4.13.2", 
        "tzdata==2025.2", 
        "urllib3==2.4.0",
        "Werkzeug==3.1.3", 
        "wheel==0.45.1", 
        "wrapt==1.14.1"
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
        "python==3.8.10",
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
