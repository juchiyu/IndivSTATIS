#' @title MFAnormMat
#'
#' @param Y A matrix to be MFA-normalized
#' @param BigShareDim FALSE (default). whether its a high dimensional matrix that requires regularized svd
#'
#' @return An MFA-normalized Y
#' @export
#'
#' @examples
#' ## a random 5 x 6 matrix
#' X <- matrix(rnorm(30), nrow = 5)
#' MFAnormMat(X)
MFAnormMat <- function (Y, BigShareDim = FALSE) {
  if (BigShareDim){
    sv = rsvd::rsvd(Y, k = 1)
  }else{sv = svd(Y)}
    sv1 = sv$d[1]
    Ynormed = Y/sv1
    return(Ynormed)
}
