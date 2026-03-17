#' Compute RV for two big matrices (with many rows)
#'
#' @param Amat a matrix with matching number of columns with Bmat
#' @param Bmat a matrix with matching number of columns with Amat
#'
#' @returns an RV coefficient between Amat and Bmat with matching columns
#' @export
#'
#' @examples
RV4BigMatrix <- function(Amat,Bmat){
  AtA <- crossprod(Amat)
  BtB <- crossprod(Bmat)
  AtB <- t(Amat) %*% Bmat
  rv.coeff <- sum(diag(crossprod(AtB)))/sqrt(sum(diag(crossprod(AtA))) * sum(diag(crossprod(BtB))))
  return(rv.coeff)
}
