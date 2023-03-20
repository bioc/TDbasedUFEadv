#' @title Generate tensor from two matrices
#'
#' @param matrix1 the first input matrix
#' @param matrix2 the second input matrix
#'
#' @return A tensor generated from the first and second matricies 
#' @export
#'
#' @examples
#' Z <- prepareTensorfromMatrix(matrix(runif(100),10),matrix(runif(100),10))
prepareTensorfromMatrix <- function(matrix1, matrix2) {
  # Argument Check
  stopifnot("`matrix1` must be a matrix." = is.matrix(matrix1))
  stopifnot("`matrix2` must be a matrix." = is.matrix(matrix2))
  #
  col_num1 <- dim(matrix1)[2]
  col_num2 <- dim(matrix2)[2]
  matrix1[is.na(matrix1)] <- 0
  matrix2[is.na(matrix2)] <- 0
  Z <- apply(
    cbind(matrix1, matrix2), 1,
    function(x) {
      outer(
        x[seq_len(col_num1)],
        x[col_num1 + seq_len(col_num2)]
      )
    }
  )
  dim(Z) <- c(col_num1, col_num2, dim(Z)[2])
  Z <- aperm(Z, c(3, 1, 2))
  return(Z)
}