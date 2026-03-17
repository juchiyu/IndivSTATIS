#' A new version metric mulidimensional scaling
#'
#' @param DistanceMatrix a distance matrix
#' @param masses NULL (default). Weights for each observations of the distance matrix.
#' @param is.distance If TRUE (default), the matrix will be processed as a distance matrix; if not, it will be processed as a covariance matrix
#' @param double.center If is.distance = FALSE, this is an option of deciding whether to double-center the matrix.
#'
#' @return
#' @export
#'
#' @examples
mmds.new <- function (DistanceMatrix, masses = NULL, is.distance = TRUE, double.center = TRUE)
{
  D <- DistanceMatrix
  nI <- nrow(D)
  if (is.null(masses)) {
    masses <- rep(1/nI, nI)
  }
  m <- masses
  LeM = t(kronecker(m, t(rep(1, nI))))
  Xhi <- diag(1, nI) - LeM

  if (is.distance){
     S <- -0.5 * sqrt(LeM) * Xhi %*% D %*% t(Xhi) * sqrt(t(LeM))
  }else{
    if (double.center){
      S <- sqrt(LeM) * Xhi %*% D %*% t(Xhi) * sqrt(t(LeM))
    }else{
      S <- D
    }
  }
  eig <- eigen(S, symmetric = TRUE)
  Cleaner = which(eig$value > 0)
  U <- eig$vector[, Cleaner]
  L <- eig$values[Cleaner]
  Nom2Dim = paste("dim", 1:length(L))
  names(L) <- Nom2Dim
  Nom2Row = rownames(D)
  LeF <- kronecker(1/sqrt(m), t(rep(1, length(L)))) * t(t(U) *
                                                          sqrt(L))
  rownames(LeF) <- Nom2Row
  colnames(LeF) <- Nom2Dim
  tau <- round(100 * (L/sum(L)), digits = 2)
  names(tau) <- Nom2Dim
  Ctr <- kronecker(m, t(rep(1, length(L)))) * t(t(LeF^2)/L)
  rownames(Ctr) <- Nom2Row
  colnames(Ctr) <- Nom2Dim
  return(list(FactorScores = LeF, eigenvalues = L, Contributions = Ctr,
              percentage = tau, M = LeM))
}
