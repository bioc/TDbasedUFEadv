#' @title Prepare Sample label for TCGA data
#'
#' @param Multi_sample list of sample ids
#' @param Clinical List of clinical data matrix from RTCGA.clinical
#' @param ID_column_of_Multi_sample Column numbers used for conditions
#' @param ID_column_of_Clinical Column numbers that include corresponding
#'  sample ids in clinical data
#'
#' @return list of sample labels
#' @export
#'
#' @examples
#' example <- readRDS(system.file(
#'   "extdata", "enrichment_example.rds",
#'   package = "TDbasedUFEadv", mustWork = TRUE
#' ))
#' cond <- prepareCondTCGA(
#'   example$Multi_sample, example$Clinical,
#'   rep(2, 4), rep(1, 4)
#' )
prepareCondTCGA <- function(Multi_sample, Clinical,
                            ID_column_of_Multi_sample,
                            ID_column_of_Clinical) {
  # Argument check
  stopifnot("`Multi_sample` must be an list." = is.list(Multi_sample))
  stopifnot("`Clinical` must be an list." = is.list(Multi_sample))
  stopifnot(
    "`ID_column_of_Clinical` must be a vector" = is.vector(ID_column_of_Clinical)
  )
  stopifnot(
    "`ID_column_of_Multi_sample` must be a vector" = is.vector(ID_column_of_Multi_sample)
  )
  stopifnot(
    "Length must be common among `Multi_sample` " = var(unlist(lapply(list(
      Multi_sample, Clinical,
      ID_column_of_Multi_sample,
      ID_column_of_Clinical
    ), length))) == 0
  )
  #
  Cond <- rep(list(NA), length(Multi_sample))
  for (i in seq_along(Multi_sample))
  {
    index <- match(
      tolower(substring(Multi_sample[[i]][, 1], 1, 12)),
      Clinical[[i]][, ID_column_of_Clinical[i]]
    )
    Cond[[i]] <- Clinical[[i]][index, ID_column_of_Multi_sample[i]]
  }
  return(Cond)
}
