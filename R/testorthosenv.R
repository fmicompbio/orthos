#' Test conda environment
#'
#' @return A list indicating whether keras is available, and the version
#'     of TensorFlow.
#'
#' @author Charlotte Soneson
#'
#' @examples
#' testOrthosEnv()
#'
#' @export
#'
#' @importFrom reticulate import
#' @importFrom basilisk basiliskStart basiliskRun basiliskStop
testOrthosEnv <- function() {
    cl <- basiliskStart(orthosenv)
    keras_tf_version <- basiliskRun(cl, function() {
        list(keras_available = keras::is_keras_available("2.10.0"),
             tf_version = tensorflow::tf$version$VERSION)
    })
    basiliskStop(cl)
    keras_tf_version
}
