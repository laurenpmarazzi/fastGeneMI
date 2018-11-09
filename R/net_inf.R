# Return the AUPRC for a matrix of pairwise mutual information values. N.reg is the number
# of genes that are known to be regulators

#' Compute the precision and recall from a matrix of mutual information estimates
#'
#' From a symmetric matrix of mutual information values between the expression of pairs of genes, infer a network then compute the precision, recall and the area under precision-recall curve by comparing the inferred network to the gold standard
#'
#' @param mim A symmetric matrix of mutual information values. The \eqn{ij^{th}} element is the mutual information between the expression of genes \eqn{i} and \eqn{j}
#' @param gs.net The symmetric adjacency matrix of the gold standard regulatory network. Should contain a 0 for no edge, 1 for an edge and \code{NA} if unknown.
#' @param inf.algo The inference algorithm used to infer the network from the mutual information estimates. Must be one of \code{"clr"}, \code{"mrnet"} or \code{"aracne"}
#' @param n.reg The number of genes that are designated as potential regulators. If it is an integer then the first \code{n.reg} genes are marked as regulators. If a vector then those gene indexes contained in the vector are regulators. Only interactions involving regulators are used to evaluate the precision and recall. Default is \code{NULL}, in which case all the genes are potential regulators.
#'
#' @return \item{pr}{A two-column matrix of the precision and recall values of the precision-recall curve. The first column contains the recall and the second column contains the precision.} \item{auprc}{The area under the precision-recall curve.}
#'
#' @export

get.pr <- function(mim, gs.net, inf.algo=c("clr", "mrnet", "aracne"), n.reg=NULL) {
  
  inf.algo <- match.arg(inf.algo)
  
  # Put in checks of dimensions
  if(!all(dim(mim)==dim(gs.net)))
    stop("mutual information matrix and gold standard must have the same dimensions")
  
  if(!isSymmetric(gs.net))
    stop("gold standard must be a symmetric matrix")
  
  n.genes <- nrow(mim)
  
  # Infer network
  inf.algo.funcs <- list("clr" = minet::clr, "mrnet" = minet::mrnet, "aracne" = minet::aracne)
  inf.algo <- inf.algo.funcs[[inf.algo]]
  adj.mat <- inf.algo(mim)
  
  # Only evaluate against regulator genes
  if(!is.null(n.reg)){
    if(length(n.reg)==1) {
      if(n.reg > n.genes)
        stop("number of regulators must be less than number of genes")
      
      adj.mat <- adj.mat[,1:n.reg]
      gs.net <- gs.net[,1:n.reg]
    }
    else {
      if(length(n.reg) > n.genes)
        stop("number of regulators must be less than the number of genes")
        
      adj.mat <- adj.mat[,n.reg]
      gs.net <- gs.net[,n.reg]
    }
  }
  
  # Get predictions and labels, removing NAs from both predictions and gold standard
  preds <- adj.mat[lower.tri(adj.mat)]
  labels <- gs.net[lower.tri(gs.net)]
  not.na.idxs <- which(!is.na(labels), arr.ind=TRUE)
  preds <- preds[not.na.idxs]
  labels <- labels[not.na.idxs]
  
  # Compute precision, recall and AUPR
  pr.curve <- PRROC::pr.curve(
    scores.class0 = preds,
    weights.class0 = labels,
    curve = TRUE
  )
  
  res <- list("pr" = pr.curve$curve,
              "auprc" = pr.curve$auc.davis.goadrich)
  return(res)
}