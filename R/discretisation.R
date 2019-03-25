# MIT License

# Copyright (c) 2018 Jonathan Ish-Horowicz

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

#' Discretise continuous expression values
#'
#' @keywords internal

get.disc.data <- function(expr.data, discretisation, n.bins, n.cores) {
  # Default number of bins is (int)n_samples^(1/3) if not using 
  # Bayesian Blocks
  if(discretisation=="bb")
    disc.data <- disc_dataset_bb_cpp(expr.data, n.cores)
  
  # Equal width or equal frequency bins
  else {
    if(as.integer(n.bins)==n.bins & n.bins > 0)
      disc.data <- infotheo::discretize(expr.data, disc = discretisation,
                                        nbins = n.bins)
    else
      stop(paste(discretisation,
                 "discretisation requires a positive integer number of bins"))
  }
  
  return(disc.data)
}

#' Bayesian Blocks discretisation of expression data
#'
#' Returns the bin edges found using the Bayesian Blocks algorithm for the expression of a single gene
#'
#' @param expr.prof A vector of expression values of a single gene.
#'
#' @return A vector of bin edges.
#'
#' @examples
#' # Compute the bin edges for the first variable of the iris
#' # dataset according to Bayesian Blocks
#' data("iris")
#' bin.edges <- get.bayesian.blocks.bins(iris[,1])
#' 
#' @references Scargle, J.D. et al, 2013. Studies in astronomical time series analysis. VI. Bayesian block representations. The Astrophysical Journal, 764(2), p.167.
#' 
#' @export

get.bayesian.blocks.bins <- function(expr.prof) {
  return(get_bb_bin_edges_cpp(expr.prof))
}