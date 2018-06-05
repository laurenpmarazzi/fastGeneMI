# fastGeneMI: an R package for computing the mutual information between the expression profiles of genes in a microarray #

## Installation

Requires `devtools`. Insatll from the R command line using

```devtools::install_bitbucket('Jonathan-Ish-Horowicz/fastGeneMI')```.

## Features

`fastGeneMI` is implements a suite of mutual information estimators in parallised C++ with an R interface.

These were implemented for inferring gene regulatory networks from DNA microarrays, but are suitable for any continuous data.

The following estimators are implemented:

* Maximum Likelihood
* Chao-Shen
* Miller-Madow
* Shrinkage
* B-spline

Also included are the inference algorithms CLR, MRNET and ARACNe (see https://www.biorxiv.org/content/early/2017/07/26/132647 for references), as well as functions to evaluate their predictions using precision-recall curves.

Full documentation is available in fastGeneMI-manual.pdf