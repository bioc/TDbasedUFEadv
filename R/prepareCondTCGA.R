#' Title Prepare Sample label for TCGA data
#'
#' @param Multi_sample list of sample ids
#' @param Clinical List of clinical data matrix from RTCGA.clinical
#' @param k Column numbers used for conditions
#' @param j Column numbers that include corresponding sample ids 
#' in clinical data
#'
#' @return list of sample labels
#' @export
#'
#' @examples
#' require(RTCGA.clinical)
#' require(RTCGA.rnaseq)
#' Clinical <- list(BLCA.clinical,BRCA.clinical,CESC.clinical,COAD.clinical)
#' Multi_sample <- list(BLCA.rnaseq[seq_len(100),1,drop=FALSE],
#'                     BRCA.rnaseq[seq_len(100),1,drop=FALSE],
#'                     CESC.rnaseq[seq_len(100),1,drop=FALSE],
#'                     COAD.rnaseq[seq_len(100),1,drop=FALSE])
#' k <- c(770,1482,773,791)
#' j <- c(20,20,12,14)
#' cond <- prepareCondTCGA(Multi_sample,Clinical,k,j)
prepareCondTCGA <- function(Multi_sample,Clinical,k,j)
{
    Cond <- rep(list(NA),length(Multi_sample))
    for (i in seq_len(length(Multi_sample)))
    {
        index <- match(tolower(substring(Multi_sample[[i]][,1],1,12)),
                       Clinical[[i]][,j[i]])
        Cond[[i]]<- Clinical[[i]][index,k[i]]
    }
    return(Cond)
}