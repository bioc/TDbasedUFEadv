#' @title Prepare condition matrix for expDrug
#'
#' @param expDrug input gene expression profile
#'
#' @return Condition matrix for expDrug
#' @export
#'
#' @examples
#' input <- readRDS(system.file(
#'   "extdata", "drug_disease_test_inputs.rds",
#'   package = "TDbasedUFEadv", mustWork = TRUE
#' ))
#' Cancer_cell_lines <- input$Cancer_cell_lines[seq_len(3)]
#' Drug_and_Disease <- prepareexpDrugandDisease(Cancer_cell_lines)
#'
#' Cond <- prepareCondDrugandDisease(Drug_and_Disease$expDrug)
prepareCondDrugandDisease <- function(expDrug) {
  # Arugument Check
  stopifnot(
    "`expDrug` must be an ExpressionSet." = is(expDrug, "ExpressionSet")
  )
  #
  dir <- system.file("extdata", package = "TDbasedUFEadv")
  TCGA <- read.csv(file.path(dir, "drug_response.txt"), sep = "\t")
  Cond <- table(TCGA$patient.arr, TCGA$drug.name)
  ID <- lapply(strsplit(as.character(colnames(expDrug)), "-"), "[", seq_len(3))
  ID <- t(data.frame(ID))
  ID <- paste(ID[, 1], ID[, 2], ID[, 3], sep = "-")
  Cond <- Cond[match(ID, rownames(Cond)), ]
  return(Cond)
}