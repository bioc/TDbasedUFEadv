% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepareTensorfromList.R
\name{prepareTensorfromList}
\alias{prepareTensorfromList}
\title{Prepare tensor from a list that includes multiple profiles}
\usage{
prepareTensorfromList(Multi, proj_dim)
}
\arguments{
\item{Multi}{a list that includes multiple profiles}

\item{proj_dim}{the number of projection dimensions}
}
\value{
a tensor as a bundle of singular value vectors obtained by 
applying  SVD to individual omics
}
\description{
Prepare tensor from a list that includes multiple profiles
}
\examples{
library(MOFAdata)
data("CLL_data")
data("CLL_covariates")
Z <- prepareTensorfromList(CLL_data,10L)
}
