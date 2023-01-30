#' Title Prepare tensor generated from two matricies that share samples
#'
#' @param sample Chracter vecor of sample names
#' @param feature list of features from two matrices
#' @param value array, contents of 
#' @param featureRange Genomic Ranges to be associated with features
#' @param sampleData List of conditional labeling associated with samples
#'
#' @return Tensor generated from two matricies that share samples
#' @export
#'
#' @examples
#' matrix1 <- matrix(runif(10000),200) #row features, column samples
#' matrix2 <- matrix(runif(20000),400) #row features, column samples
#' Z <- prepareTensorfromMatrix(t(matrix1),t(matrix2))
#' Z <- PrepareSummarizedExperimentTensorRect(sample=as.character(seq_len(50)),
#' feature=list(as.character(seq_len(200)),as.character(seq_len(400))),
#' sampleData=list(rep(seq_len(2),each=25)),value=Z)
PrepareSummarizedExperimentTensorRect <- function(sample,feature,value,
    featureRange=GRanges(NULL),sampleData=list(NULL)){
    new("SummarizedExperimentTensorRect",sample=sample,
        feature=feature,value=value,
        featureRange=featureRange,sampleData=sampleData)
}