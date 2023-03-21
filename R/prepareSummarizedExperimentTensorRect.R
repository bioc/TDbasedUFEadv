#' @title Prepare tensor generated from two matrices that share samples
#'
#' @param sample Character vector of sample names
#' @param feature list of features from two matrices
#' @param value array, contents of 
#' @param featureRange Genomic Ranges to be associated with features
#' @param sampleData List of conditional labeling associated with samples
#'
#' @return Tensor generated from two matrices that share samples
#' @export
#'
#' @examples
#' matrix1 <- matrix(runif(1000),200) #row features, column samples
#' matrix2 <- matrix(runif(2000),400) #row features, column samples
#' Z <- prepareTensorfromMatrix(t(matrix1),t(matrix2))
#' Z <- prepareSummarizedExperimentTensorRect(sample=as.character(seq_len(50)),
#' feature=list(as.character(seq_len(200)),as.character(seq_len(400))),
#' sampleData=list(rep(seq_len(2),each=25)),value=Z)
prepareSummarizedExperimentTensorRect <- function(
    sample, feature, value,
    featureRange = GRanges(NULL), sampleData = list(NULL)) {
  # Argument check
  stopifnot("`sample` must be a character." = is.character(sample))
  stopifnot("`feature` must be a list." = is.list(feature))
  stopifnot("`value` must be a array." = is.array(value))
  stopifnot("`featureRange` must be an GRanges." = is(featureRange, "GRanges"))
  stopifnot("`sampleData` must be a list." = is.list(sampleData))
  #
  new("SummarizedExperimentTensorRect",
    sample = sample,
    feature = feature, value = value,
    featureRange = featureRange, sampleData = sampleData
  )
}