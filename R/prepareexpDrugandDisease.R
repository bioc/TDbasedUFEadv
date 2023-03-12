#' Title Generating gene expression of drug treated cell lines and a disease
#' cell line
#'
#' @param LIST list that includes indivisual data set from RTCGA.rnaseq 
#'
#' @return list of expDrug and expDisease 
#' @export
#'
#' @examples
#'  \donttest{
#' require(RTCGA.rnaseq)
#' LIST <- list(ACC.rnaseq,
#'             BLCA.rnaseq,
#'             BRCA.rnaseq)
#' dummy <- prepareexpDrugandDisease(LIST)
#' }
prepareexpDrugandDisease <- function(LIST){
    dir <- system.file("extdata", package="TDbasedUFEadv")
    TCGA <- read.csv(file.path(dir,"drug_response.txt"), sep="\t")
    TCGAdb <- unique(TCGA[,1])
    if (length(TCGAdb)>16) TCGAdb <- TCGAdb[-16]
    RTCGA_all <- matrix(NA, nrow=dim(TCGA)[1],ncol=dim(LIST[[1]])[2])
    LISTp <- rep(list(NA),length(LIST))
    toTCGA <- function(x){
        ID <- lapply(strsplit(as.character(x[,1]),"-"),"[",seq_len(3))
        ID <- t(data.frame(ID))
        ID <- paste(ID[,1],ID[,2],ID[,3],sep="-")
        return( x[ID %in% TCGA[,2],])
    }
    LISTp <- lapply(LIST,toTCGA)
    expDrug <- NULL
    for (i in seq_len(length(LISTp)))
    {
        expDrug <- rbind(expDrug,data.frame(LISTp[[i]]))
    }
    expDrug <- convertTCGA(expDrug)
    ID <- lapply(strsplit(as.character(LIST[[3]][,1]),"-"),"[",seq_len(3))
    ID <- t(data.frame(ID))
    ID <- paste(ID[,1],ID[,2],ID[,3],sep="-")
    expDisease <- LIST[[3]][!(ID %in% TCGA[,2]),]
    expDisease <- convertTCGA(expDisease)
    return(list(expDrug=expDrug,expDisease=expDisease))
}