#' Title Prepare tensor from a list that includes multiple profiles
#'
#' @param Multi  a list that includes multiple profiles
#' @param L the number of projection dimensions
#'
#' @return a tensor as a bundle of singular value vectors obtained by 
#' applying  SVD to individual omics
#' @export
#'
#' @examples
#' require(MOFAdata)
#' data("CLL_data")
#' data("CLL_covariates")
#' Z <- prepareTensorfromList(CLL_data,10)
prepareTensorfromList <- function(Multi,L)
{
    LIST <- lapply(Multi,function(x){x[is.na(x)]<-0;return(x)})
    SVD_lst <- lapply(LIST,function(x){SVD <- svd(x,L);return(t(SVD$v[,seq_len(L)]))})
    Z <- array(unlist(SVD_lst),c(dim(SVD_lst[[1]]),length(Multi)))
    return(Z)
}