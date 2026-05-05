.convertTCGA <- function(dataSet) {
    if (!is.data.frame(dataSet)) {
        stop("'dataSet' must be a data.frame.", call. = FALSE)
    }

    if (ncol(dataSet) < 2L) {
        stop("'dataSet' must have at least two columns.", call. = FALSE)
    }

    rnaseqMatrix <- t(as.matrix(dataSet[, -1, drop = FALSE]))
    colnames(rnaseqMatrix) <- as.character(dataSet[[1]])

    Biobase::ExpressionSet(assayData = rnaseqMatrix)
}
