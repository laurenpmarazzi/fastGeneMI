// ----------------------------------------------------------------------------------
//  These functions get empirical probability distributions from discretised data
// ----------------------------------------------------------------------------------

#include "fastGeneMI.h"

using namespace Rcpp;

// Return empirical marginal distributions of a discretised dataset as an
// std::vector of arma::vec
arma::vec get_emp_marg_dist(const arma::Mat<int>& disc_data_col)
{
  const int n_samples(disc_data_col.n_rows);
  
  const int n_bins = disc_data_col.max() + 1;
  arma::vec p_marginal(n_bins);
    
  // Loop through each bin
  for(int k(0); k<n_bins; ++k)
  {
    arma::umat s = arma::sum(disc_data_col==k);
    arma::mat s2 = arma::conv_to<arma::mat>::from(s);
    p_marginal(k) = s2(0)/(double)n_samples;
  }
  return p_marginal;
}

// Return the empirical joint distributions of a discrete dataset
arma::mat get_emp_joint_dist(const arma::Mat<int>& disc_data_col_i, const arma::Mat<int>& disc_data_col_j)
{
  const int n_samples(disc_data_col_i.n_rows);
  
  // This requires every bin to have at least one entry
  const int n_bins_i = disc_data_col_i.max() + 1; 
  const int n_bins_j = disc_data_col_j.max() + 1;
        
  arma::mat counts2d = arma::mat(n_bins_i, n_bins_j, arma::fill::zeros);
  for(int k(0); k<n_samples; ++k)
  {
    counts2d(disc_data_col_i(k), disc_data_col_j(k)) += 1.0;
  }

  return counts2d/(double)n_samples;
}
