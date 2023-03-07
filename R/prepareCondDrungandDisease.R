#' Title Prepare condition matrix for expDrug
#'
#' @param expDrug input gene expression profile
#'
#' @return Condition matrix for expDrug
#' @export
#'
#' @examples
#' \donttest{
#' require(RTCGA.rnaseq)
#'LIST <- list(ACC.rnaseq,
#'             BLCA.rnaseq,
#'             BRCA.rnaseq,
#'             CESC.rnaseq,
#'             COAD.rnaseq,
#'             ESCA.rnaseq,
#'             GBM.rnaseq,
#'             HNSC.rnaseq,
#'             KICH.rnaseq,
#'             KIRC.rnaseq,
#'             KIRP.rnaseq,
#'             LGG.rnaseq,
#'             LIHC.rnaseq,
#'             LUAD.rnaseq,
#'             LUSC.rnaseq,
#'             OV.rnaseq,
#'             PAAD.rnaseq,
#'             PCPG.rnaseq,
#'             PRAD.rnaseq,
#'             READ.rnaseq,
#'             SARC.rnaseq,
#'             SKCM.rnaseq,
#'             STAD.rnaseq,
#'             TGCT.rnaseq,
#'             THCA.rnaseq,
#'             UCEC.rnaseq,
#'             UCS.rnaseq)
#' dummy <- prepareexpDrugandDisease(LIST)
#' expDrug <- dummy[[1]]
#' Cond <- prepareCondDrugandDisease(expDrug)
#' }
prepareCondDrugandDisease <- function(expDrug){
    dir <- system.file("extdata", package="TDbasedUFEadv")
    TCGA <- read.csv(file.path(dir,"drug_response.txt"), sep="\t")
    Cond <- table(TCGA[,2],TCGA[,3])
    ID <- lapply(strsplit(as.character(colnames(expDrug)),"-"),"[",seq_len(3))
    ID <- t(data.frame(ID))
    ID <- paste(ID[,1],ID[,2],ID[,3],sep="-")
    Cond <- Cond[match(ID,rownames(Cond)),]
    return(Cond)
}