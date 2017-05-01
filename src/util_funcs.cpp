// ----------------------------------------------------------------------------------
//  Utility functions
// ----------------------------------------------------------------------------------

#include "fastGeneMI.h"

using namespace Rcpp;

// ----------------------------------------------------------------------------------
//  For converting between R and C++/Armadillo objects
// ----------------------------------------------------------------------------------

// Convert a NumericMatrix (R object) to an Armadillo matrix (C++ object)
arma::mat R2armaMat_num(NumericMatrix rMat)
{
  arma::mat armaMat = arma::mat(rMat.begin(), rMat.nrow(), rMat.ncol(), false);
  return armaMat;
}

// Convert a NumericVector (R object) to an Armadillo vector (C++ object)
arma::vec R2armaVec_num(NumericVector rVec)
{
  arma::vec armaVec = arma::vec(rVec.begin(), rVec.length(), false);
  return armaVec;
}

// ----------------------------------------------------------------------------------
//  Maximum Likelihood Entropy functions
// ----------------------------------------------------------------------------------


// Maximum Likelihood marginal entropy from marginal probability
double get_marginal_ml_entropy(const arma::vec& p_marg)
{
  return -arma::sum(p_marg % arma::log(p_marg + 1e-16));
}

// Maximum likelihood joint entropy from joint probability
double get_joint_ml_entropy(const arma::mat& p_joint)
{
  return -arma::accu(p_joint % arma::log(p_joint + 1e-16));
}

// ----------------------------------------------------------------------------------
//  Build a matrix whose entries are the index of flattened 1d array of
//  the upper-triangle of a n_genes x n_genes matrix
// ----------------------------------------------------------------------------------

arma::Mat<int> get_idx_lookup_mat(const int n_genes)
{
  arma::Mat<int> idx_lookup = arma::Mat<int>(n_genes, n_genes, arma::fill::zeros);
  int ij(0);
  for(int i(0); i<n_genes; ++i)
  {
    for(int j(i); j<n_genes; ++j)
    {
      idx_lookup(i,j) = ij;
      ++ij;
    }
  }
  return idx_lookup;
}


