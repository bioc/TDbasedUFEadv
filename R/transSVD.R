#' @title Convert SVD to that for the case where 
#' samples are shared between two matrices
#'
#' @param SVD input SVD object generated from computeSVD function
#'
#' @return converted SVD objects
#' @export
#'
#' @examples
#' matrix1 <- matrix(runif(200),20)
#' matrix2 <- matrix(runif(400),20)
#' SVD <- computeSVD(matrix1,matrix2)
#' SVD <- transSVD(SVD)
transSVD <- function(SVD) {
  # Augument check
  stopifnot("`SVD` must be a list." = is.list(SVD))
  #
  u0 <- SVD$SVD$u
  v0 <- SVD$SVD$v
  SVD$SVD$u <- scale(SVD$u)
  SVD$SVD$v <- scale(SVD$v)
  SVD$u <- scale(u0)
  SVD$v <- scale(v0)
  return(SVD)
}