
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DivComAnalyses

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/fconstancias/DivComAnalyses.svg?branch=master)](https://travis-ci.com/fconstancias/DivComAnalyses)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/fconstancias/DivComAnalyses?branch=master&svg=true)](https://ci.appveyor.com/project/fconstancias/DivComAnalyses)
[![Codecov test
coverage](https://codecov.io/gh/fconstancias/DivComAnalyses/branch/master/graph/badge.svg)](https://codecov.io/gh/fconstancias/DivComAnalyses?branch=master)
<!-- badges: end -->

The goal of DivComAnalyses is to facilitate community analysis (e.g.,
16S community metabarcoding, metaphlan taxonomic profiles from
metagenomes). In data filtering (i.e., potential contaminants identified
using decontam R package), help in normalizing data (i.e., 16S qPCR
normalization), generate phylogenetic and taxonomic alpha-diversity
metrics, visualize community composition and structure using ordination,
facilitate statistical analysis based on DEseq2, ALDex, vegan R
packages.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("fconstancias/DivComAnalyses")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(DivComAnalyses)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub\!
