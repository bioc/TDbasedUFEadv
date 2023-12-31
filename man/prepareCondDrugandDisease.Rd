% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepareCondDrugandDisease.R
\name{prepareCondDrugandDisease}
\alias{prepareCondDrugandDisease}
\title{Prepare condition matrix for expDrug}
\usage{
prepareCondDrugandDisease(expDrug)
}
\arguments{
\item{expDrug}{input gene expression profile}
}
\value{
Condition matrix for expDrug
}
\description{
Prepare condition matrix for expDrug
}
\examples{
\donttest{
library(RTCGA.rnaseq)
Cancer_cell_lines <- list(ACC.rnaseq,BLCA.rnaseq,BRCA.rnaseq)
Drug_and_Disease <- prepareexpDrugandDisease(Cancer_cell_lines)
Cond <- prepareCondDrugandDisease(Drug_and_Disease$expDrug)
}
}
