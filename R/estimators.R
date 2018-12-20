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

# Wrappers for the estimator functions in RcppExports.R

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
#' @examples
#' # Compute the mutual information matrix for the variables in the 
#' # iris dataset using default discretisation parameters
#' # (N^1/3 equal width bins for N samples). Computation is performed
#' # sequentially be default (i.e. using a single core)
#' data("iris")
#' expr.data <- iris[,1:4]
#' mim <- get.mim.ML(expr.data, "equalwidth")
#' # Now use N^1/2 bins and use 2 cores
#' mim <- get.mim.ML(expr.data, "equalwidth", as.integer(nrow(expr.data)^0.5),
#'                    n.cores=2)
#' 
#' @references Scargle, J.D. et al, 2013. Studies in astronomical time series analysis. VI. Bayesian block representations. The Astrophysical Journal, 764(2), p.167.
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
#' @examples
#' # Compute the mutual information matrix for the variables in the 
#' # iris dataset using default discretisation parameters
#' # (N^1/3 equal width bins for N samples). Computation is performed
#' # sequentially be default (i.e. using a single core)
#' data("iris")
#' expr.data <- iris[,1:4]
#' mim <- get.mim.MM(expr.data, "equalwidth")
#' # Now use N^1/2 bins and use 2 cores
#' mim <- get.mim.MM(expr.data, "equalwidth", as.integer(nrow(expr.data)^0.5),
#'                    n.cores=2)
#' 
#' @references Miller, G., 1955. Note on the bias of information estimates. Information Theory in Psychology: Problems and Methods.
#' @references Scargle, J.D. et al, 2013. Studies in astronomical time series analysis. VI. Bayesian block representations. The Astrophysical Journal, 764(2), p.167.
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
#' @examples
#' # Compute the mutual information matrix for the variables in the 
#' # iris dataset using default discretisation parameters
#' # (N^1/3 equal width bins for N samples). Computation is performed
#' # sequentially be default (i.e. using a single core)
#' data("iris")
#' expr.data <- iris[,1:4]
#' mim <- get.mim.CS(expr.data, "equalwidth")
#' # Now use N^1/2 bins and use 2 cores
#' mim <- get.mim.CS(expr.data, "equalwidth", as.integer(nrow(expr.data)^0.5),
#'                    n.cores=2)
#' # Infer the regulatory network using CLR. No evaluation is performed
#' # as no gold standard is provided
#' net <- infer.net(mim, "clr")
#' 
#' @references Chao, A. and Shen, T.J., 2003. Nonparametric estimation of Shannon's index of diversity when there are unseen species in sample. Environmental and Ecological Statistics, 10(4), pp.429-443.
#' @references Scargle, J.D. et al, 2013. Studies in astronomical time series analysis. VI. Bayesian block representations. The Astrophysical Journal, 764(2), p.167.
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
#' @examples
#' # Compute the mutual information matrix for the variables in the 
#' # iris dataset using default discretisation parameters
#' # (N^1/3 equal width bins for N samples). Computation is performed
#' # sequentially be default (i.e. using a single core)
#' data("iris")
#' expr.data <- iris[,1:4]
#' mim <- get.mim.shrink(expr.data, "equalwidth")
#' # Now use N^1/2 bins and use 2 cores
#' mim <- get.mim.shrink(expr.data, "equalwidth", as.integer(nrow(expr.data)^0.5),
#'                    n.cores=2)
#' # Infer the regulatory network using CLR. No evaluation is performed
#' # as no gold standard is provided
#' net <- infer.net(mim, "clr")
#' 
#' @references Hausser, J. and Strimmer, K., 2009. Entropy inference and the James-Stein estimator, with application to nonlinear gene association networks. Journal of Machine Learning Research, 10(Jul), pp.1469-1484.
#' @references Scargle, J.D. et al, 2013. Studies in astronomical time series analysis. VI. Bayesian block representations. The Astrophysical Journal, 764(2), p.167. 
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
#' @references Daub, C.O. et al, 2004. Estimating mutual information using B-spline functions-an improved similarity measure for analysing gene expression data. BMC Bioinformatics, 5(1), p.118.
#'
#' @return A symmetric matrix of pairwise mutual information between the expression of pairs of genes. These are the mutual information estimtates between pairs of columns of \code{expr.data}.
#'
#' @examples
#' # Compute the mutual information matrix for the variables in the 
#' # iris dataset using default discretisation parameters
#' # (N^1/3 equal width bins for N samples) and spline order (2).
#' # Computation is performed sequentially be default (i.e. using a
#' # single core)
#' data("iris")
#' expr.data <- iris[,1:4]
#' mim <- get.mim.bspline(expr.data, 2)
#' # Now use N^1/2 bins and use 2 cores
#' mim <- get.mim.bspline(expr.data, 2, as.integer(nrow(expr.data)^0.5),
#'                        n.cores=2)
#' # Infer the regulatory network using CLR. No evaluation is performed
#' # as no gold standard is provided
#' net <- infer.net(mim, "clr")
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

