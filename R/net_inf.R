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

#' Infer a regulatory network from a matrix of mutual information values between pairs of genes using one of CLR, MRNET or ARACNe
#'
#' From a symmetric matrix of mutual information values between the expression of pairs of genes, infer a regulatory network. If a gold standard network is provided the inferred network is evaluated and the resulting precision, recall and the area under precision-recall is computed and returned along with the inferred network.
#'
#' @param mim A symmetric matrix of mutual information values. The \eqn{ij^{th}} element is the mutual information between the expression of genes \eqn{i} and \eqn{j}
#' @param inf.algo The inference algorithm used to infer the network from the mutual information estimates. Must be one of \code{"clr"}, \code{"mrnet"} or \code{"aracne"}. Default is \code{"CLR"}.
#' @param gs.net The symmetric adjacency matrix of the gold standard regulatory network. Should contain a 0 for no edge, 1 for an edge and \code{NA} if unknown. Note that unknown edges are not included in the evaluation.
#' @param n.reg The number of genes that are designated as potential regulators. If it is an integer then the first \code{n.reg} genes are marked as regulators. If a vector then those gene indexes contained in the vector are regulators. Only interactions involving regulators are used to compute the precision and recall. Default is \code{NULL}, in which case all the genes are potential regulators (note that this is typically not the biological reality).
#' @param plot Logical controlling whether or not to return the precision-recall curve.
#'
#' @return \item{network}{The inferred regulatory network. A matrix whose \eqn{ij^{th}} element represents the confidence of an edge between genes \eqn{i} and \eqn{j}.} \item{pr}{A two-column matrix of the precision and recall values of the precision-recall curve resulting from evaluating the inferred network against a gold standard. The first column contains the recall and the second column contains the precision. Only returned if a gold standard network is provided.} \item{auprc}{The area under the precision-recall curve. Only returned if a gold standard network is provided.} \item{plot}{The precision-recall curve.}
#'
#' @references Faith, J.J. et al, 2007. Large-scale mapping and validation of Escherichia coli transcriptional regulation from a compendium of expression profiles. PLoS Biology, 5(1), p.e8.
#' @references Meyer, P.E. et al, 2007. Information-theoretic inference of large transcriptional regulatory networks. EURASIP Journal on Bioinformatics and Systems Biology, 2007, pp.8-8.
#' @references Margolin, A.A. et al, 2006, March. ARACNE: an algorithm for the reconstruction of gene regulatory networks in a mammalian cellular context. BMC Bioinformatics (Vol. 7, No. 1, p. S7).
#' 
#'
#' @export

infer.net <- function(mim, inf.algo=c("clr", "mrnet", "aracne"), gs.net=NULL, n.reg=NULL, plot=FALSE)
{
  # Infer network
  inf.algo <- match.arg(inf.algo)
  inf.algo.funcs <- list("clr" = minet::clr, "mrnet" = minet::mrnet, "aracne" = minet::aracne)
  inf.algo <- inf.algo.funcs[[inf.algo]]
  network <-inf.algo(mim)
  predictions <- network
  
  # Return just the inferred network if no gold standard is provided
  if(is.null(gs.net))
    return(network)
  else # Or evaluate against gold standard if one is provided
  {
    # Checks argument dimensions
    if(!all(dim(mim)==dim(gs.net)))
      stop("mutual information matrix and gold standard must have the same dimensions")
    
    if(!isSymmetric(gs.net))
      stop("gold standard must be a symmetric matrix - only undirected networks are supported")
    
    # Check gold standard is binary
    if( (length(gs.net[is.na(gs.net)]) + length(gs.net[gs.net==0]) + length(gs.net[gs.net==1]) != length(gs.net)) )
      stop("gold standard can only contain 0, 1 or NA")
    
    n.genes <- nrow(mim)
    predictions <- network
    
    # Only evaluate against regulator genes
    if(!is.null(n.reg)){
      if(length(n.reg)==1) {
        if(n.reg > n.genes)
          stop("number of regulators must be less than number of genes")
        
        predictions <- predictions[,1:n.reg]
        gs.net <- gs.net[,1:n.reg]
      }
      else {
        if(length(n.reg) > n.genes)
          stop("number of regulators must be less than the number of genes")
        
        predictions <- predictions[,n.reg]
        gs.net <- gs.net[,n.reg]
      }
    }
    else
      warning("All genes are designated as potential regulators - this is not recommended.")
    
    # Get predictions and labels, removing NAs from both predictions and gold standard
    predictions <- predictions[lower.tri(predictions)]
    labels <- gs.net[lower.tri(gs.net)]
    not.na.idxs <- which(!is.na(labels), arr.ind=TRUE)
    predictions <- predictions[not.na.idxs]
    labels <- labels[not.na.idxs]
    
    # Compute precision, recall and AUPR
    pr.curve <- PRROC::pr.curve(
      scores.class0 = predictions,
      weights.class0 = labels,
      curve = TRUE
    )
    
    return(list(network = network,
                pr = pr.curve$curve,
                auprc = pr.curve$auc.davis.goadrich,
                plot = plot(pr.curve)))
  }
}

#' Compute the precision and recall of a network inferred from pairwise mutual information estimates by comparing to a gold standard
#'
#' From a symmetric matrix of mutual information values between the expression of pairs of genes, infer a network then compute the precision, recall and the area under precision-recall curve by comparing the inferred network to the gold standard. If the regulatory relationship between two genes is unknown then use \code{NA} as the appropriate element of \code{gs.net}. These are ignnored when computing the precision and recall. This is a thin wrapper around \code{\link{infer.net}} that does not return the inferred network and is included for backward compatibility with a previous version.
#' 
#' @param mim A symmetric matrix of mutual information values. The \eqn{ij^{th}} element is the mutual information between the expression of genes \eqn{i} and \eqn{j}
#' @param gs.net The symmetric adjacency matrix of the gold standard regulatory network. Should contain a 0 for no edge, 1 for an edge and \code{NA} if unknown. \code{NA} values are ignored when computing the precsicion and recall.
#' @param inf.algo The inference algorithm used to infer the network from the mutual information estimates. Must be one of \code{"clr"}, \code{"mrnet"} or \code{"aracne"}
#' @param n.reg The number of genes that are designated as potential regulators. If it is an integer then the first \code{n.reg} genes are marked as regulators. If a vector then those gene indexes contained in the vector are regulators. Only interactions involving regulators are used to evaluate the precision and recall. Default is \code{NULL}, in which case all the genes are potential regulators (note that this is typically not the biological reality).
#'
#' @return \item{pr}{A two-column matrix of the precision and recall values of the precision-recall curve. The first column contains the recall and the second column contains the precision.} \item{auprc}{The area under the precision-recall curve.}
#'
#' @references Faith, J.J. et al, 2007. Large-scale mapping and validation of Escherichia coli transcriptional regulation from a compendium of expression profiles. PLoS Biology, 5(1), p.e8.
#' @references Meyer, P.E. et al, 2007. Information-theoretic inference of large transcriptional regulatory networks. EURASIP Journal on Bioinformatics and Systems Biology, 2007, pp.8-8.
#' @references Margolin, A.A. et al, 2006, March. ARACNE: an algorithm for the reconstruction of gene regulatory networks in a mammalian cellular context. BMC Bioinformatics (Vol. 7, No. 1, p. S7).
#'
#' @export

get.pr <- function(mim, gs.net, inf.algo=c("clr", "mrnet", "aracne"), n.reg=NULL) {
  warning("This function is deprecated - use infer.net instead")
  out <- infer.net(mim=mim, inf.algo=inf.algo, gs.net=gs.net, n.reg=n.reg)
  return(list(pr = out$pr, auprc = out$auprc))
}