# fastGeneMI: an R package for computing the mutual information between the expression profiles of genes in a microarray #

## Requirements

* Rcpp 
* RcppArmadillo
* infotheo
* PRROC
* minet
* testthat

## Installation

Requires `devtools`. Install from R using

```devtools::install_github("jonathanishhorowicz/fastGeneMI")```.

## Features

`fastGeneMI` is a package that implements a suite of mutual information estimators in parallelised C++ with an R interface.

These were implemented for inferring gene regulatory networks from DNA microarrays but are suitable for any continuous data.

The following estimators are implemented:

* Maximum Likelihood
* Chao-Shen
* Miller-Madow
* Shrinkage
* B-spline

Also included are the inference algorithms CLR, MRNET and ARACNe, as well as functions to evaluate their predictions using precision-recall curves. Included data is from the DREAM5 Network Inference Challenge. See the Supplementary Material of the accompanying paper at https://www.biorxiv.org/content/early/2017/07/26/132647 for descriptions of each estimator.

Full documentation is available in fastGeneMI-manual.pdf
