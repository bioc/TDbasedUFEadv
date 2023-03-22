#' Title Perform SVD toward reduced matrix generated from a tensor with 
#' partial summation
#'
#' @param matrix1 The first original matrix that generates a tensor
#' @param matrix2 The second original matrix that generates a tensor
#' @param dim  The number of singular value vectors to be computed
#' @param scale If matrix should be scaled or not
#'
#' @return Singular value vectors attributed to two sets of objects associated 
#' with singular value vectors attributed to features, by multiplying  
#' @export
#'
#' @examples
#' matrix1 <- matrix(runif(200),20)
#' matrix2 <- matrix(runif(400),20)
#' SVD <- computeSVD(matrix1,matrix2)
computeSVD <- function(matrix1, matrix2, dim = 10L, scale = TRUE) {
  # Argument check
  stopifnot("`matrix1` must be a matrix." = is.matrix(matrix1))
  stopifnot("`matrix2` must be a matrix." = is.matrix(matrix2))
  stopifnot("`dim` must be a integer." = is.integer(dim))
  stopifnot("`scale` must be a logical." = is.logical(scale))
  stopifnot(
    "# of rows must be common between `matrix1` and `matrix2`" 
        = nrow(matrix1) == nrow(matrix2)
  )
  #
  matrix1[is.na(matrix1)] <- 0
  matrix2[is.na(matrix2)] <- 0
  Z <- t(matrix1) %*% matrix2
  if (scale) Z <- Z / mean(Z)
  SVD <- svd(Z, dim, dim)
  u <- matrix1 %*% SVD$u
  v <- matrix2 %*% SVD$v
  return(list(SVD = SVD, u = u, v = v))
}