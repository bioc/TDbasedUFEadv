
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TDbasedUFEadv

<!-- badges: start -->
<!-- badges: end -->

The goal of TDbasedUFEadv is to …

## Installation



``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("TDbasedUFEadv")
```

You can install the latest release of TDbasedUFEadv from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("tagtag/TDbasedUFEadv@0.99.0-1")
```

You can install the development version of TDbasedUFEadv from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("tagtag/TDbasedUFEadv",build_vignettes = TRUE)
```

or if it does not work, try

``` r
# install.packages("devtools")
devtools::install_github("tagtag/TDbasedUFEadv")
```

## Introduction

It is an advanced version of
[TDbasedUFE](https://bioconductor.org/packages/devel/bioc/html/TDbasedUFE.html).
Thus people who would like to use TDbasedUFEadv, please first install
and check TDbasedUFE.

## Vignettes

How to use it.

vignette(“QuickStart”)

vignette(“QuickStart2”)

vignette(“Enrichment”)

For more theoretical background

vignette(“TDbasedUFEadv”)
