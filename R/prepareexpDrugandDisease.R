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
#' }
prepareexpDrugandDisease <- function(LIST){
    dir <- system.file("extdata", package="TDbasedUFEadv")
    TCGA <- read.csv(file.path(dir,"drug_response.txt"), sep="\t")
    TCGAdb <- unique(TCGA[,1])
    TCGAdb <- TCGAdb[-16]
    RTCGA_all <- matrix(NA, nrow=dim(TCGA)[1],ncol=dim(LIST[[1]])[2])
    LISTp <- rep(list(NA),length(LIST))
    toTCGA <- function(x){
        ID <- lapply(strsplit(as.character(x[,1]),"-"),"[",seq_len(3))
        ID <- t(data.frame(ID))
        ID <- paste(ID[,1],ID[,2],ID[,3],sep="-")
        return( x[ID %in% TCGA[,2],])
    }
    LISTp <- lapply(LIST,toTCGA)
    expDrug <- rbind(data.frame(LISTp[[1]]),data.frame(LISTp[[2]]),
        data.frame(LISTp[[3]]),data.frame(LISTp[[4]]),data.frame(LISTp[[5]]),
        data.frame(LISTp[[6]]),data.frame(LISTp[[7]]),data.frame(LISTp[[8]]),
        data.frame(LISTp[[9]]),data.frame(LISTp[[10]]),data.frame(LISTp[[11]]),
        data.frame(LISTp[[12]]),data.frame(LISTp[[13]]),data.frame(LISTp[[14]]),
        data.frame(LISTp[[15]]),data.frame(LISTp[[16]]),data.frame(LISTp[[17]]),
        data.frame(LISTp[[18]]),data.frame(LISTp[[19]]),data.frame(LISTp[[20]]),
        data.frame(LISTp[[21]]),data.frame(LISTp[[22]]),data.frame(LISTp[[23]]),
        data.frame(LISTp[[24]]),data.frame(LISTp[[25]]),data.frame(LISTp[[26]]),
                     data.frame(LISTp[[27]]))
    expDrug <- convertTCGA(expDrug)
    ID <- lapply(strsplit(as.character(LIST[[3]][,1]),"-"),"[",seq_len(3))
    ID <- t(data.frame(ID))
    ID <- paste(ID[,1],ID[,2],ID[,3],sep="-")
    expDisease <- LIST[[3]][!(ID %in% TCGA[,2]),]
    expDisease <- convertTCGA(expDisease)
    return(list(expDrug=expDrug,expDisease=expDisease))
}