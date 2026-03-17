#' Compute RV from a list with big matrices
#'
#' @param data.list a list of matrices
#'
#' @returns an RV matrix that stores the pair-wise RV coefficients of tables in the data.list
#' @export
#'
#' @examples
GetCmat_for_BigShareDim <- function(data.list){
  N.dim <- length(data.list)
  Cmat <- matrix(nrow = N.dim, ncol = N.dim)

  for (mat.walk1 in 1:(N.dim-1)){
    for (mat.walk2 in c((1 + mat.walk1):N.dim)){
      Cmat[mat.walk1, mat.walk2] = RV4BigMatrix(data.list[[mat.walk1]], data.list[[mat.walk2]])
    }
  }
  diag(Cmat) = 1
  Cmat[lower.tri(Cmat)] <- Cmat[upper.tri(Cmat)]
  return(Cmat)
}
