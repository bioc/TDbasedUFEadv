#' Title
#'
#' @slot sample character. 
#' @slot feature list. 
#' @slot value array. 
#' @slot featureRange GRanges. 
#' @slot sampleData list. 
setClass("SummarizedExperimentTensorRect",slot=list(sample="character",
                                                feature="list",value="array",
                                                featureRange="GRanges",sampleData="list"))