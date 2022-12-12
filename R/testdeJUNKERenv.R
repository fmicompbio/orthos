#' Test conda environment
#'
#' @return Nothing is returned
#' @author Charlotte Soneson
#'
#' @examples
#' testdeJUNKERenv()
#'
#' @export
#'
#' @importFrom reticulate import
#' @importFrom basilisk basiliskStart basiliskRun basiliskStop
testdeJUNKERenv <- function() {
    cl <- basiliskStart(dejunkerenv)
    basiliskRun(cl, function() {
        K <- keras::backend()
        tensorflow::tf$version$VERSION
    })
    basiliskStop(cl)
}
