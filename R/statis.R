# Preamble ----
#_____________________________________________________________________
# statis: function statis
# Private functions used in this file
#  >>> See file names matching the function names
# -- private functions
#   o MFAnormMat
#   o CP2MFAnormedMat
#   o GetCmat
#   o ComputeSplus
#   o rdiag & ldiag
#________________________
# Last update  09 / 06 / 2022 by Ju-Chi
#_____________________________________________________________________
# statis Preamble ----
#' multitable method with  "STATIS" optimization
#' procedure for asymmetric, rectangular matrices
#' @title statis
#'
#' @param LaGrandeTable A list of observations
#' \eqn{\times}{*} variables matrices. These matrices describe the same
#' observations with sets of variables (which do not necessarily have
#' to match across matrices).
#' @param Norm Type of normalization
#' used for each matrix.
#' Current options are \code{NONE} (do nothing),
#' \code{SUMPCA} (normalize by the square root of the total sum of squares)
#' or \code{MFA} (\code{default}) that normalizes each matrix so
#' that its first eigenvalue is equal to one.
#' @param center If \code{TRUE},
#' (\code{default}) the columns of each matrix are centered
#' (i.e., to have a mean of 0).
#' If \code{FALSE}, the columns will \emph{not} be centered.
#' @param scale If \code{TRUE}, the columns of each matrix are normalized
#' to have a standard deviation equal 1.
#' If \code{"SS1"}, (\code{default}) the columns of each matrix are
#' normalized to have a sum of squares of 1,
#' so that they contribute equally.
#' If \code{FALSE}, the columns will \emph{not} be normalized.
#' @param RV if \code{TRUE} (\code{default})
#' we use the \eqn{R_V}{Rv} coefficient to
#' compute the \eqn{\alpha}{weights},
#' if \code{FALSE}
#' we use the matrix scalar product
#' @param nfact2keep (default: \code{3}) Number of factors
#' to keep for the computation of the
#' factor scores of the observations.
#' @param compact if \code{FALSE} (default),
#'  \code{statis} provides detailed output, if
#'  \code{TRUE},  \code{statis} sends back
#' only the \eqn{\alpha}{alpha} weights
#'  (this option is used to make the
#' bootstrap routine compute RV for two big matrices (with many rows).
#' @param BigShareDim whether the shared dimensions between tables are huge.
#' If \code{TRUE}, the function will use \code{GetCmat_for_BigShareDim}
#' to compute the RV matrix. Default: \code{FALSE}
#' @param haveC A predetermined RV matrix, if there exist one. This can be useful
#' for big datasets, when it takes a long time to generate the RV matrix.
#' \code{\link{BootFromCompromise}} more efficient).
#' @return  \code{statis} sends back the results
#' \emph{via} two lists:
#' \code{res.Cmat}
#' and \code{res.Splus}.
#' Note that items with a * are the only ones sent back
#' when using the \code{compact = TRUE} option.
#'
#' @return A list with:
#' \describe{
#' \item{res.Cmat}{Results for the between
#'           distance matrices analysis.}
#' \itemize{
#'   \item \code{res.Cmat$C}
#'   The \eqn{I\times I}{I*I} \bold{C} matrix
#'   of scalar products (or \eqn{R_V}{Rv} between matrices).
#'   \item
#'   \code{res.Cmat$vectors} The eigenvectors of the \bold{C} matrix
#'   \item
#'  \code{res.Cmat$alpha} * The \eqn{\alpha}{alpha} weights
#'  \item
#'  \code{res.Cmat$value} The eigenvalues of the \bold{C} matrix
#'  \item
#'  \code{res.Cmat$G} The factor scores for the \bold{C} matrix
#'  \item
#'   \code{res.Cmat$ctr} The contributions for \code{res.Cmat$G},
#'  \item
#'    \code{res.Cmat$cos2} The squared cosines for \code{res.Cmat$G}
#'  \item
#'  \code{res.Cmat$d2} The squared
#'  Euclidean distance  for \code{res.Cmat$G}.
#'    }
#'
#' \item{res.Splus}{Results for the between observation analysis.}
#' \itemize{
#' \item \code{res.Splus$SCP} an \eqn{I\times I\times K}{I*I*K} array.
#' Contains
#' the (normalized if needed)
#' cross product matrices corresponding to the
#' distance matrices.
#' \item \code{res.Splus$Splus} * The compromise
#' (optimal linear
#' combination of the SCP's').
#'  \item \code{res.Splus$eigValues} *
#'   The eigenvalues of the compromise).
#'  \item \code{res.Splus$eigVectors} *
#'   The eigenvectors of the compromise).
#' \item \code{res.Splus$tau} * The percentage
#' of explained inertia of the eigenValues).
#' \item \code{res.Splus$Fi} The factor scores for the rows.
#' \item \code{res.Splus$Fj} The factor scores for the columns.
#' \item \code{res.Splus$Q} The generalized singular vectors (singular vectors Weighted by alphas) for the columns.
#'  \item
#'   \code{res.Splus$ctr.i} The contributions for \code{res.Cmat$Fi}.
#'   \item
#'   \code{res.Splus$ctr.j} The contributions for \code{res.Cmat$Fj}.
#'  \item
#'    \code{res.Splus$cos2.i} The squared cosines for \code{res.Cmat$Fi}.
#'  \item
#'    \code{res.Splus$cos2.j} The squared cosines for \code{res.Cmat$Fj}.
#'  \item
#'  \code{res.Splust$d2.i} The squared
#'  Euclidean distance  for \code{res.Cmat$Fi}.
#'  \item
#'  \code{res.Splust$d2.j} The squared
#'  Euclidean distance  for \code{res.Cmat$Fj}.
#' \item \code{res.Splus$PartialFi} an
#' \eqn{I \times \code{nf2keep} \times K}{I*nf2keep*K} array.
#' Contains the partial factors for the distance
#' matrices.
#' \item \code{res.Splus$ProjectionMatrix} The
#' projection matrix used to compute factor
#' scores and partial factor scores.
#' \item \code{res.Splus$Splus} The
#' alpha-weighted linear combination of the cross-product matrices that
#' creates the compromise.
#' }
#'}
#'
#' @importFrom ExPosition expo.scale
#' @export
#'
#' @examples
#' ## Wine data from Abdi et al. (2012)
#' ## A list of tables with matching rows and different variables
#' data(wines2012List)
#'
#' ## This replicate results from the paper
#' statis.res <- statis(wines2012List, Norm = "SUMPCA", scale = TRUE, RV = FALSE)
statis <- function(LaGrandeTable,
                   Norm = 'MFA',
                   center = TRUE,
                   scale = "SS1",
                   RV = TRUE,
                   nfact2keep = 3,
                   BigShareDim = TRUE,
                   haveC = NULL,
                   compact = FALSE) {

  ## center and scale the data tables
  LaGrandeTable.preproc <- lapply(LaGrandeTable, function(x) expo.scale(x, center = center, scale = scale))

  # perform MFA normalization ----
  if (Norm == 'MFA') {
    LaGrandeTable.preproc <- lapply(LaGrandeTable.preproc, MFAnormMat, BigShareDim = BigShareDim)
  }
  if (Norm == 'SUMPCA') {
    LaGrandeTable.preproc <- lapply(LaGrandeTable.preproc, SUMPCAnormMat)
  }

  if (is.null(haveC)){
    if (BigShareDim){
      C <- GetCmat_for_BigShareDim(LaGrandeTable.preproc)
    }else{
      ## Put data into an array
      if (is.list(LaGrandeTable.preproc)){
        diff.row <- length(unique(lapply(LaGrandeTable.preproc, nrow)))
        diff.col <- length(unique(lapply(LaGrandeTable.preproc, ncol)))

        if (diff.row == 1 & diff.col > 1 | diff.row == 1 & diff.col == 1 & diff.row < diff.col){
          ## create crossproduct with rows
          data.list <- lapply(LaGrandeTable.preproc, function(x) crossprod(t(x)))
        }else if (diff.row > 1 & diff.col == 1 | diff.row == 1 & diff.col == 1 & diff.row > diff.col){
          data.list <- lapply(LaGrandeTable.preproc, crossprod)
        }else{
          stop("Either the rows or the columns of the data tables have to match.")
        }
        leCube <- array(as.numeric(unlist(data.list)),
                        dim=c(dim(data.list[[1]])[1], dim(data.list[[1]])[2], length(data.list)),
                        dimnames = list(rownames(data.list[[1]]), colnames(data.list[[1]]), names(data.list)))
      }else if(is.array(LaGrandeTable.preproc)){
        leCube <- LaGrandeTable.preproc
      }else{
        stop("The tables should be in a list or an array.")
      }

      # Compute C matrix ----
      C <- DistatisR::GetCmat(leCube, RV = RV) # (to match matlab implementation of distatis)
    }
  }else{
    # If there's a C to read
    C <- haveC
  }
  # eigen of C ----
  eigC <- eigen(C, symmetric = TRUE) # Eigen-decomposition of C
  # alpha weights ----
  alpha <- eigC$vectors[, 1] / sum(eigC$vectors[, 1])

  if (compact == FALSE) {
    # All C stuff ----
    eigC$vectors[, 1] <- abs(eigC$vectors[, 1])
    rownames(eigC$vectors) <- rownames(C)
    nfC <- ncol(eigC$vectors) # number of eigenvectors of C
    Nom2Dim <- paste('dim', 1:nfC)
    colnames(eigC$vectors) <- Nom2Dim
    # make sure that first eigenvector is positive
    eigC$tau <- round(100 * (eigC$values / sum(eigC$values)))
    names(eigC$tau)    <- Nom2Dim
    names(eigC$values) <- Nom2Dim
    # factors scores for RV mat
    eigC$G   <- t(apply(eigC$vectors, 1, '*', t(t(sqrt(abs(eigC$values) )))))
    rownames(eigC$G) <- rownames(C)
    colnames(eigC$G) <- Nom2Dim
    G2        <- eigC$G^2
    # new C stuff ----
    eigC$ctr  <- rdiag(G2, 1/eigC$values)
    eigC$d2G  <- rowSums(G2)
    eigC$cos2 <- ldiag(1/eigC$d2G, G2)

    # Save the output ----
    res.Cmat <- list(
      C = C,
      eigVector = eigC$vector,
      eigValues = eigC$values,
      tau = eigC$tau,
      G = eigC$G,
      ctr  = eigC$ctr,
      cos2 = eigC$cos2,
      d2   = eigC$d2G,
      alpha = alpha,
      compact = compact
    )
    class(res.Cmat) <- c("Cmat", "list")
  }else{
    res.Cmat <- list(alpha = alpha, compact = compact)
    class(res.Cmat) <- c("Cmat", "list")
  }

  if (BigShareDim == TRUE){
    ## Create weighted grand table
    BigX <- WeightCbind(LaGrandeTable.preproc, sqrt(alpha))
    rownames(BigX) <- rownames(LaGrandeTable.preproc[[1]])
    colnames(BigX) <- unlist(lapply(LaGrandeTable.preproc, colnames))
    if (compact == FALSE){
      # BigX stuff ----
      # SVD of BigX ----
      svdBigX <- rsvd::rsvd(BigX, k = nfact2keep)
      # Eigenvalues
      svdBigX$eig <- svdBigX$d^2
      # Percent of Inertia
      svdBigX$tau <- round(100 *
                             svdBigX$eig / sum(BigX^2))
      # Get the factor scores
      Fi <- t(t(svdBigX$u)*svdBigX$d)
      rownames(Fi) <- rownames(BigX)

      ## Add on ctr + cos
      Fi2 <- Fi^2
      d2Fi    <-  rowSums(Fi2)
      ctrFi   <- rdiag(Fi2, 1/ abs(svdBigX$eig) )
      cos2Fi  <- ldiag(1/d2Fi, Fi2)
      # Keep only the interesting factors
      Nom2Factors <- paste('Factor', 1:nfact2keep)
      colnames(Fi) <- Nom2Factors
      colnames(ctrFi)  <- Nom2Factors
      colnames(cos2Fi) <- Nom2Factors

      # Column factor scores ----
      nCol4Tab <- unlist(lapply(LaGrandeTable.preproc, ncol))
      alphaForCol <- rep(sqrt(alpha), nCol4Tab)
      qj <- svdBigX$v/alphaForCol
      colnames(qj) <- Nom2Factors
      rownames(qj) <- colnames(BigX)

      fj.noalpha <- t(t(qj)*svdBigX$d[1:nfact2keep])
      fj <- fj.noalpha*alphaForCol

      Fj <- Fj_noalpha <- Q <- LaGrandeTable.preproc
      for (i in 1:length(Q)){
        if (i == 1){
          nstart = 1
        }else{
          nstart = sum(nCol4Tab[1:(i-1)])+1
        }
        nend = sum(nCol4Tab[1:i])
        Q[[i]] <- qj[nstart:nend,]
        Fj_noalpha[[i]] <- fj.noalpha[nstart:nend,]
        Fj[[i]] <- fj[nstart:nend,]
      }

      Fj2     <-  lapply(Fj, '^', 2)
      d2Fj    <-  lapply(Fj, rowSums)
      ctrFj   <- lapply(Fj2, function(x) rdiag(x, 1/ abs(svdBigX$eig)))
      cos2Fj  <- list()
      for (ln in 1:length(d2Fj)){
        cos2Fj[[ln]] <- ldiag(1/d2Fj[[ln]], Fj2[[ln]])
      }
      names(cos2Fj) <- names(Fj2)

      # Get the partial projections
      # PartialF ----
      PartialFi <-  array(dim = c(dim(BigX)[1], nfact2keep, length(LaGrandeTable)))
      for (i in 1:length(LaGrandeTable)){
        PartialFi[,,i] <- LaGrandeTable.preproc[[i]] %*% Q[[i]]
      }
      rownames(PartialFi) <- rownames(BigX)
      colnames(PartialFi) <- Nom2Factors
      dimnames(PartialFi)[[3]] <- names(LaGrandeTable)

      # Save the output ----
      res.Splus <- list(
        GrandTable = BigX,
        svd = svdBigX,
        eigValues  = svdBigX$eig,
        tau   = svdBigX$tau,
        Fi = Fi,
        Fj = Fj,
        Q = Q,
        ctr.i = ctrFi,
        ctr.j = ctrFj,
        cos2.i = cos2Fi,
        cos2.j = cos2Fj,
        d2.i = d2Fi,
        d2.j = d2Fj,
        PartialFi = PartialFi,
        BigX = BigX,
        compact = compact
      )
      class(res.Splus) <- c("Splus", "list")
      res.statis <- list(res4Cmat = res.Cmat,
                         res4Splus = res.Splus,
                         compact = compact,
                         params = list(
                           Norm = Norm,
                           RV = RV))
      class(res.statis) <- c("StatisR", "list")
    }else{
      # When "compact" is TRUE, send back only the compact information
      res.Splus <- list(Splus = Splus, compact = compact)
      class(res.Splus) <- c("Splus", "list")
      res.statis <- list(res4Cmat = res.Cmat,
                         res4Splus = res.Splus,
                         compact = compact)
      class(res.statis) <- c("StatisR", "list")
    }
  }else{
    # compute compromise ----
    Splus <- DistatisR::ComputeSplus(leCube, alpha)
    if (compact == FALSE) {
      # S+ stuff ----
      # eigen of S+ ----
      # Eigen decomposition of Splus
      eigenSplus <- eigen(Splus, symmetric = TRUE)
      # Percent Of Inertia
      eigenSplus$tau <- round(100 *
                                eigenSplus$values / sum(eigenSplus$values))
      # singular values
      eigenSplus$SingularValues <- sqrt(abs(eigenSplus$values))
      # Get the factor scores (ugly way, but works)
      Fi <- t(apply(eigenSplus$vectors, 1, '*', t(t(
        eigenSplus$SingularValues ))))
      rownames(Fi) <- rownames(Splus)
      # Add on ctr + cos.
      # new S stuff----
      Fi2     <-  Fi^2
      d2Fi    <-  rowSums(Fi2)
      ctrFi   <- rdiag(Fi2, 1/ abs(eigenSplus$values) )
      cos2Fi  <- ldiag(1/d2Fi, Fi2)
      # Keep only the interesting factors
      Nom2Factors <- paste('Factor', 1:nfact2keep)
      Fi <- Fi[, 1:nfact2keep]
      colnames(Fi) <- Nom2Factors
      ctrFi            <- ctrFi[, 1:nfact2keep]
      colnames(ctrFi)  <- Nom2Factors
      cos2Fi           <- cos2Fi[, 1:nfact2keep]
      colnames(cos2Fi) <- Nom2Factors
      # Projection matrix ----
      ProjMat <- t(apply(eigenSplus$vectors, 1, '*',
                         1 / t(t(eigenSplus$SingularValues))))
      Proj <- ProjMat[, 1:nfact2keep]
      colnames(Proj) <- Nom2Factors
      rownames(Proj) <- rownames(Splus)
      # Get the partial projections
      # PartialF ----
      PartialFi <- array(apply(leCube, 3, '%*%', Proj),
                         dim = c(dim(leCube)[1], nfact2keep, dim(leCube)[3]))
      rownames(PartialFi) <- rownames(Splus)
      colnames(PartialFi) <- Nom2Factors
      dimnames(PartialFi)[[3]] <- rownames(C)
      # Column factor scores ----

      Q <- lapply(LaGrandeTable.preproc, function(x){
        qj <- t(x) %*% eigenSplus$vectors[,1:nfact2keep] %*% diag(1/eigenSplus$SingularValues[1:nfact2keep])
        colnames(qj) <- Nom2Factors
        return(qj)
      })

      Fj_noalpha <- lapply(LaGrandeTable.preproc, function(x){ ## missing alpha
        fj <- t(x) %*% eigenSplus$vectors[,1:nfact2keep]
        colnames(fj) <- Nom2Factors
        return(fj)
      })
      Fj <- mapply('*', Fj_noalpha, sqrt(alpha), SIMPLIFY = FALSE)

      Fj2     <-  lapply(Fj, '^', 2)
      d2Fj    <-  lapply(Fj, rowSums)
      ctrFj   <- lapply(Fj2, function(x) rdiag(x, 1/ abs(eigenSplus$values[1:nfact2keep])))
      cos2Fj  <- list()
      for (ln in 1:length(d2Fj)){
        cos2Fj[[ln]] <- ldiag(1/d2Fj[[ln]], Fj2[[ln]])
      }
      names(cos2Fj) <- names(Fj2)
      # pack up the information to send back.
      # May try as some point to keep a structure similar
      # to MExPosition
      # in the meantime go for fast and match original matlab program.
      # pack return list ----
      res.Splus <- list(
        SCP = leCube,
        eigValues  = eigenSplus$values,
        eigVectors = eigenSplus$vectors,
        tau   = eigenSplus$tau,
        Fi = Fi,
        Fj = Fj,
        Q = Q,
        ctr.i = ctrFi,
        ctr.j = ctrFj,
        cos2.i = cos2Fi,
        cos2.j = cos2Fj,
        d2.i = d2Fi,
        d2.j = d2Fj,
        PartialFi = PartialFi,
        ProjectionMatrix = Proj,
        Splus = Splus,
        compact = compact
      )
      class(res.Splus) <- c("Splus", "list")
      res.statis <- list(res4Cmat = res.Cmat,
                         res4Splus = res.Splus,
                         compact = compact,
                         params = list(
                           Norm = Norm,
                           RV = RV))
      class(res.statis) <- c("StatisR", "list")
      # End of if compact == FALSE
    } else {
      # When "compact" is TRUE, send back only the compact information
      res.Splus <- list(Splus = Splus, compact = compact)
      class(res.Splus) <- c("Splus", "list")
      res.statis <- list(res4Cmat = res.Cmat,
                         res4Splus = res.Splus,
                         compact = compact)
      class(res.statis) <- c("StatisR", "list")
    }
  }
  # return ----
  return(res.statis) # et voila ----
}
