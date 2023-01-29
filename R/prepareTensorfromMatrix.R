#' Title Generate tensor from two matricies
#'
#' @param matrix1 the first input matrix
#' @param matrix2 the second input matrix
#'
#' @return A tensor generated from the first and second matricies 
#' @export
#'
#' @examples
#' Z <- prepareTensorfromMatrix(matrix(runif(100),10),matrix(runif(100),10))
prepareTensorfromMatrix <-function(matrix1,matrix2){
    L1 <- dim(matrix1)[2]
    L2 <- dim(matrix2)[2]
    matrix1[is.na(matrix1)] <-0 
    matrix1[is.na(matrix2)] <-0 
    Z <- apply(cbind(matrix1,matrix2),1,function(x){outer(x[seq_len(L1)],x[L1+seq_len(L2)])})
    dim(Z) <- c(L1,L2,dim(Z)[2])
    Z <- aperm(Z,c(3,1,2))
    return(Z)
}