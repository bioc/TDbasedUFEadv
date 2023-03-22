#' @title Class definitions
#'
#' @slot sample character. 
#' @slot feature list. 
#' @slot value array. 
#' @slot featureRange GRanges. 
#' @slot sampleData list. 
setClass("TensorRect",slot=list(sample="character",
                                                feature="list",value="array",
                                                featureRange="GRanges",sampleData="list"))