#' @title Generating gene expression of drug treated cell lines and a disease
#' cell line
#'
#' @param Cancer_cell_lines <- list(ACC.rnaseq,BLCA.rnaseq,BRCA.rnaseq)
#'  list that includes individual data set from RTCGA.rnaseq 
#'
#' @return list of expDrug and expDisease 
#' @export
#'
#' @examples
#'  \donttest{
#' require(RTCGA.rnaseq)
#' Cancer_cell_lines <- list(ACC.rnaseq,BLCA.rnaseq,BRCA.rnaseq)
#' Drug_and_Disease <- prepareexpDrugandDisease(Cancer_cell_lines)
#' }
prepareexpDrugandDisease <- function(Cancer_cell_lines) {
  # Argument check
  stopifnot(
    "`Cancer_cell_lines` must be an list." = is.list(Cancer_cell_lines)
  )
  #
  dir <- system.file("extdata", package = "TDbasedUFEadv")
  TCGA <- read.csv(file.path(dir, "drug_response.txt"), sep = "\t")
  TCGAdb <- unique(TCGA$cancers)
  if (length(TCGAdb) > 16) TCGAdb <- TCGAdb[-16]
  RTCGA_all <- matrix(NA,
    nrow = dim(TCGA)[1],
    ncol = dim(Cancer_cell_lines[[1]])[2]
  )
  Cancer_cell_lines_p <- rep(list(NA), length(Cancer_cell_lines))
  toTCGA <- function(x) {
    ID <- lapply(strsplit(as.character(x[, 1]), "-"), "[", seq_len(3))
    ID <- t(data.frame(ID))
    ID <- paste(ID[, 1], ID[, 2], ID[, 3], sep = "-")
    return(x[ID %in% TCGA$patient.arr, ])
  }
  Cancer_cell_lines_p <- lapply(Cancer_cell_lines, toTCGA)
  expDrug <- NULL
  for (i in seq_len(length(Cancer_cell_lines_p)))
  {
    expDrug <- rbind(expDrug, data.frame(Cancer_cell_lines_p[[i]]))
  }
  expDrug <- convertTCGA(expDrug)
  ID <- lapply(
    strsplit(as.character(Cancer_cell_lines[[3]][, 1]), "-"),
    "[", seq_len(3)
  )
  ID <- t(data.frame(ID))
  ID <- paste(ID[, 1], ID[, 2], ID[, 3], sep = "-")
  expDisease <- Cancer_cell_lines[[3]][!(ID %in% TCGA$patient.arr), ]
  expDisease <- convertTCGA(expDisease)
  return(list(expDrug = expDrug, expDisease = expDisease))
}