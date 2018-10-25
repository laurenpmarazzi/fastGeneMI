// ----------------------------------------------------------------------------------
//  Utility functions
// ----------------------------------------------------------------------------------

#include "fastGeneMI.h"

using namespace Rcpp;

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

int get_n_gene_pairs(const int n_genes)
{
    return (int)(((double)n_genes/2.0) * (double)(n_genes+1));
}

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

std::vector<std::pair <int,int> > get_ij_list(const int n_genes)
{
  const int n_pairs = get_n_gene_pairs(n_genes);
  std::vector<std::pair <int,int> > ij_list(n_pairs);
  int ij(0);

  for(int i=0; i<n_genes; ++i)
  {
    for(int j=i; j<n_genes; ++j)
    {
      ij_list[ij] = std::make_pair(i,j);
      ++ij;
    }
  }

  return ij_list;
}  



