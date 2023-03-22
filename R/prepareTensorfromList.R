#' @title Prepare tensor from a list that includes multiple profiles
#'
#' @param Multi  a list that includes multiple profiles
#' @param proj_dim the number of projection dimensions
#'
#' @return a tensor as a bundle of singular value vectors obtained by 
#' applying  SVD to individual omics
#' @export
#'
#' @examples
#' library(MOFAdata)
#' data("CLL_data")
#' data("CLL_covariates")
#' Z <- prepareTensorfromList(CLL_data,10L)
prepareTensorfromList <- function(Multi, proj_dim) {
  # Argument check
  stopifnot("`Multi` must be a list." = is.list(Multi))
  stopifnot("`Proj_dim` must be an integer." = is.integer(proj_dim))
  #
  Multi_list <- lapply(Multi, function(x) {
    x[is.na(x)] <- 0
    return(x)
  })
  SVD_lst <- lapply(
    Multi_list,
    function(x) {
      SVD <- svd(x, proj_dim)
      return(t(SVD$v[, seq_len(proj_dim)]))
    }
  )
  Z <- array(unlist(SVD_lst), c(dim(SVD_lst[[1]]), length(Multi)))
  return(Z)
}