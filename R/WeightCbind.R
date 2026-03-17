#' Weighting the tables in a list and concatenating them by columns
#'
#' @param List A list of K data matrices with matching rows
#' @param weight K weights, one for each table
#'
#' @returns A big matrix with all weighted tables concatenated by columns
#' @export
#'
#' @examples
WeightCbind <- function(List, weight){
  weight <- as.vector(weight)
  if (length(List) != length(weight)){
    stop("The number of tables does not match with the length of weight.")
  }
  no.Row <- nrow(List[[1]])
  weightedList <- Map("*", List, weight)
  CbindList <- do.call(cbind, weightedList)
  return(CbindList)
}
