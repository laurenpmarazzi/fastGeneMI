# Wrappers for the functions in RcppExports.R

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

#' Maximum Likelihood mutual information estimate
#'
#' Returns a symmetric matrix of pairwise Maximum Likelihood mutual information estimates from a samples x genes matrix of expression values. The \eqn{ij^{th}} element is the mutual information estimate between the expression of genes \eqn{i} and \eqn{j}.
#'
#' @param expr.data A samples x genes matrix of expression values.
#' @param discretisation The types of bins into which the expression data will be discretised. Can choose equal width bins (\code{"equalwidth"}), equal frequency bins (\code{"equalfreq"}) or Bayesian Blocks bins (\code{"bb"})
#' @param n.bins The number of bins into which the expression data will be discretised. The default is \eqn{N^{1/3}} for \eqn{N} samples. Ignored if using Bayesian Blocks discretisation.
#' @param n.cores The number of cores to use for the computation. The default is 1.
#'
#' @return A symmetric matrix of pairwise mutual information between the expression of pairs of genes. These are the mutual information estimtates between pairs of columns of \code{expr.data}.
#'
#' @export

get.mim.ML <- function(expr.data, discretisation=c("equalwidth", "equalfreq", "bb"), n.bins=as.integer(nrow(expr.data)^(1/3)), n.cores=1L) {
  discretisation <- match.arg(discretisation)
  if(!as.integer(n.cores)==n.cores)
    stop("n.cores must be an integer")
  disc.data <- get.disc.data(expr.data, discretisation, n.bins, n.cores)
  mim <- mim_ML_cpp(as.matrix(disc.data), n.cores)
  rownames(mim) <- colnames(expr.data)
  colnames(mim) <- colnames(expr.data)
  return(mim)
}

#' Miller-Madow mutual information estimate
#'
#' Returns a symmetric matrix of pairwise Miller-Madow mutual information estimates from a samples x genes matrix of expression values. The \eqn{ij^{th}} element is the mutual information estimate between the expression of genes \eqn{i} and \eqn{j}.
#'
#' @param expr.data A samples x genes matrix of expression values.
#' @param discretisation The types of bins into which the expression data will be discretised. Can choose equal width bins (\code{"equalwidth"}), equal frequency bins (\code{"equalfreq"}) or Bayesian Blocks bins (\code{"bb"})
#' @param n.bins The number of bins into which the expression data will be discretised. The default is \eqn{N^{1/3}} for \eqn{N} samples. Ignored if using Bayesian Blocks discretisation.
#' @param n.cores The number of cores to use for the computation. The default is 1.
#'
#' @return A symmetric matrix of pairwise mutual information between the expression of pairs of genes. These are the mutual information estimtates between pairs of columns of \code{expr.data}.
#'
#' @export

get.mim.MM <- function(expr.data, discretisation=c("equalwidth", "equalfreq", "bb"), n.bins=as.integer(nrow(expr.data)^(1/3)), n.cores=1L) {
  discretisation <- match.arg(discretisation)
  if(!as.integer(n.cores)==n.cores)
    stop("n.cores must be an integer")
  disc.data <- get.disc.data(expr.data, discretisation, n.bins, n.cores)
  mim <- mim_MM_cpp(as.matrix(disc.data), n.cores)
  rownames(mim) <- colnames(expr.data)
  colnames(mim) <- colnames(expr.data)
  return(mim)
}

#' Chao-Shen Likelihood mutual information estimate
#'
#' Returns a symmetric matrix of pairwise Chao-Shen mutual information estimates from a samples x genes matrix of expression values. The \eqn{ij^{th}} element is the mutual information estimate between the expression of genes \eqn{i} and \eqn{j} in nats.
#'
#' @param expr.data A samples x genes matrix of expression values.
#' @param discretisation The types of bins into which the expression data will be discretised. Can choose equal width bins (\code{"equalwidth"}), equal frequency bins (\code{"equalfreq"}) or Bayesian Blocks bins (\code{"bb"})
#' @param n.bins The number of bins into which the expression data will be discretised. The default is \eqn{N^{1/3}} for \eqn{N} samples. Ignored if using Bayesian Blocks discretisation.
#' @param n.cores The number of cores to use for the computation. The default is 1.
#'
#' @return A symmetric matrix of pairwise mutual information between the expression of pairs of genes. These are the mutual information estimtates between pairs of columns of \code{expr.data}.
#'
#' @export

get.mim.CS <- function(expr.data, discretisation=c("equalwidth", "equalfreq", "bb"), n.bins=as.integer(nrow(expr.data)^(1/3)), n.cores=1L) {
  discretisation <- match.arg(discretisation)
  if(!as.integer(n.cores)==n.cores & !n.cores>0)
    stop("n.cores must be a positive integer")
  disc.data <- get.disc.data(expr.data, discretisation, n.bins, n.cores)
  mim <- mim_CS_cpp(as.matrix(disc.data), n.cores)
  rownames(mim) <- colnames(expr.data)
  colnames(mim) <- colnames(expr.data)
  return(mim)
}

#' Shrinkage mutual information estimate
#'
#' Returns a symmetric matrix of pairwise Shrinkage mutual information estimates from a samples x genes matrix of expression values. The \eqn{ij^{th}} element is the mutual information estimate between the expression of genes \eqn{i} and \eqn{j} in nats.
#'
#' @param expr.data A samples x genes matrix of expression values.
#' @param discretisation The types of bins into which the expression data will be discretised. Can choose equal width bins (\code{"equalwidth"}), equal frequency bins (\code{"equalfreq"}) or Bayesian Blocks bins (\code{"bb"})
#' @param n.bins The number of bins into which the expression data will be discretised. The default is \eqn{N^{1/3}} for \eqn{N} samples. Ignored if using Bayesian Blocks discretisation.
#' @param n.cores The number of cores to use for the computation. The default is 1.
#'
#' @return A symmetric matrix of pairwise mutual information between the expression of pairs of genes. These are the mutual information estimtates between pairs of columns of \code{expr.data}.
#'
#' @export

get.mim.shrink <- function(expr.data, discretisation=c("equalwidth", "equalfreq", "bb"), n.bins=as.integer(nrow(expr.data)^(1/3)), n.cores=1L) {
  discretisation <- match.arg(discretisation)
  if(!as.integer(n.cores)==n.cores)
    stop("n.cores must be an integer")
  disc.data <- get.disc.data(expr.data, discretisation, n.bins, n.cores)
  mim <- mim_shrink_cpp(as.matrix(disc.data), n.cores)
  rownames(mim) <- colnames(expr.data)
  colnames(mim) <- colnames(expr.data)
  return(mim)
}

#' B-spline mutual information estimate
#'
#' Returns a symmetric matrix of pairwise B-spline mutual information estimates from a samples x genes matrix of expression values. The \eqn{ij^{th}} element is the mutual information estimate between the expression of genes \eqn{i} and \eqn{j} in nats.
#'
#' @param expr.data A samples x genes matrix of expression values.
#' @param order The order of the B-spline used to smooth the histograms. This must be greater than 1 + the number of bins.
#' @param n.bins The number of bins into which the expression data will be discretised. The default is \eqn{N^{1/3}} for \eqn{N} samples.
#' @param n.cores The number of cores to use for the computation. The default is 1.
#'
#' @return A symmetric matrix of pairwise mutual information between the expression of pairs of genes. These are the mutual information estimtates between pairs of columns of \code{expr.data}.
#'
#' @export

get.mim.bspline <- function(expr.data, order, n.bins=as.integer(nrow(expr.data)^(1/3)), n.cores=1L) {
  # Check arguments are valid
  if(!as.integer(order)==order)
    stop("Spline order must be an integer")

  if(!as.integer(n.bins)==n.bins)
    stop("number of bins must be an integer")
  
  if(!as.integer(n.cores)==n.cores)
    stop("n.cores must be an integer")
  
  if(all(n.bins+1<=order))
    stop(paste("spline order must be less than number of bins + 1"))
  
  mim <- mim_bspline_cpp(as.matrix(expr.data), order, n.bins, n.cores)
  rownames(mim) <- colnames(expr.data)
  colnames(mim) <- colnames(expr.data)
  return(mim)
}

# Discretise expression data
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

