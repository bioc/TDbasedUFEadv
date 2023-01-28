#' Title Prepare condition matrix for expDrug
#'
#' @param expDrug input gene expression profile
#'
#' @return Condition matrix for expDrug
#' @export
#'
#' @examples
#' dummy <- prepareexpDrugandDisease()
#' expDrug <- dummy[[1]]
#' Cond <- prepareCondDrugandDisease(expDrug)
prepareCondDrugandDisease <- function(expDrug){
    dir <- system.file("extdata", package="TDbasedUFEadv")
    TCGA <- read.csv(file.path(dir,"drug_response.txt"), sep="\t")
    Cond <- table(TCGA[,2],TCGA[,3])
    ID <- lapply(strsplit(as.character(colnames(expDrug)),"-"),"[",1:3)
    ID <- t(data.frame(ID))
    ID <- paste(ID[,1],ID[,2],ID[,3],sep="-")
    Cond <- Cond[match(ID,rownames(Cond)),]
    return(Cond)
}