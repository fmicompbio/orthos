#' Test conda environment
#'
#' @return A list indicating whether keras is available, and the version
#'     of TensorFlow.
#'
#' @author Charlotte Soneson, Panagiotis Papasaikas
#'
#' @examples
#' testOrthosEnv()
#'
#' @export
#'
#' @importFrom reticulate import
#' @importFrom basilisk basiliskStart basiliskRun basiliskStop
testOrthosEnv <- function() {
    cl <- basiliskStart(orthosenv,
                        testload = "tensorflow")
    keras_tf_version <- basiliskRun(cl, function() {
        
        # Set up env variable for explicit report of libload and debugging
        Sys.setenv(
            LD_DEBUG = "libs",         # prints every dlopen() to stderr
            TF_CPP_MIN_LOG_LEVEL = "0" # show all TF log messages
        )
        
        ## capture stderr into a temp-file 
        tf_trace <- tempfile("ld_debug_", fileext = ".log")
        trace_con <- file(tf_trace, open = "wt")
        sink(trace_con, type = "message")      # redirect C-level stderr
        on.exit({ sink(type = "message"); close(trace_con) }, add = TRUE)
        
        
        
        # Collect diagnostics
        tf_version <- tensorflow::tf$version$VERSION
        build_info <- tensorflow::tf$sysconfig$get_build_info()
        devices <- tensorflow::tf$config$list_physical_devices()
        keras_ok <- keras::is_keras_available("2.10.0")
        ld_path <- Sys.getenv("LD_LIBRARY_PATH")
        
        
        ## stop capturing, read & filter trace
        sink(type = "message")   # restore stderr
        close(trace_con)
        
        trace_txt <- readLines(tf_trace, warn = FALSE)
        missing   <- grep("=> .*ENOENT|undefined symbol", trace_txt,
                          value = TRUE)
        
        
        # Print diagnostic info
        cat("=== TensorFlow/Keras Diagnostic ===\n")
        cat("TF Version: ", tf_version, "\n")
        cat("LD_LIBRARY_PATH:\n", ld_path, "\n")
        cat("Build Info:\n")
        print(build_info)
        cat("Physical Devices:\n")
        print(devices)
        cat("Keras Available: ", keras_ok, "\n")
        cat("---- Missing/failed libraries (ENOENT / undefined) ----\n")
        if (length(missing)) {
            cat(paste0(missing, collapse = "\n"), "\n")
        } else {
            cat("None – linker found every requested .so \n")
        }
        cat("===================================\n")
        
        list(
            keras_available = keras_ok,
            build_info = build_info, 
            tf_version = tf_version,
            cuda_version = build_info$cuda_version,
            cudnn_version = build_info$cudnn_version,
            devices = unlist(lapply(devices, function(x) as.character(x$name))),
            ld_library_path = ld_path
        )
    })
    basiliskStop(cl)
    keras_tf_version
        

}
