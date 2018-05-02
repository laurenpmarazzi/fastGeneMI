// MIT License

// Copyright (c) 2018 Jonathan Ish-Horowicz

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

// ----------------------------------------------------------------------------------
//  Implementations of the Mutual Information estimators using Maximum likelihood,
//  Miller-Madow, Chao-Shen, Jack-knififed (unvalidated) and Shrinkage entropy 
//  estimation (problem in joint entropy computation). For references see 
//  Ish-Horowicz and Reid, (2018)
// ----------------------------------------------------------------------------------

#include "fastGeneMI.h"

using namespace Rcpp;

// ----------------------------------------------------------------------------------
//  Maximum Likelihood, Miller-Madow, Chao-Shen and Shrinkage Mutual Information
//  Estimators
// ----------------------------------------------------------------------------------

// Maximum likelihood mutual information

// [[Rcpp::export]]
arma::mat mim_ML_cpp(NumericMatrix disc_expr_data, int n_cores)
{
  arma::Mat<int> data = arma::conv_to<arma::Mat<int> >::from(R2armaMat_num(disc_expr_data));
  const int n_genes(data.n_cols), n_samples(data.n_rows);
  const int n_pairs = get_n_gene_pairs(n_genes);

  //Rcout << "There are " << n_genes << " genes and "
        //<< n_samples << " samples" << std::endl;
  
  // Change from R indexing to C++ indexing
  data -= 1;
  
  // Compute marginal entropies
  //Rcout << "Computing marginal entropies...";
  std::vector<double> h_marginals(n_genes);
  for(int j; j<n_genes; ++j)
  {
    arma::vec p_marginal = get_emp_marg_dist(data.col(j));
    h_marginals[j] = get_marginal_ml_entropy(p_marginal);
  }
  //Rcout << "done" << std::endl;
    
  // Compute joint entropies in parallel
  //Rcout << "Computing joint entropies...";
  std::vector<double> h_joints(n_pairs); // Change to armadillo vector for returning to R
  
  // Number of cores to use
  // const int n_cores = sysconf(_SC_NPROCESSORS_ONLN);
  omp_set_num_threads(n_cores);
  
  //const int n_threads = omp_get_num_threads();
  //Rcout << "There are " << n_threads << " threads\n";

  const std::vector<std::pair <int,int> > ij_pairs = get_ij_list(n_genes);

  #pragma omp parallel for shared(data, h_joints)
  for(int ij=0; ij<n_pairs; ++ij)
  {
    std::pair <int,int> ij_pair = ij_pairs[ij];
    int i(ij_pair.first), j(ij_pair.second);
    arma::mat p_joint = get_emp_joint_dist(data.col(i), data.col(j));
    h_joints[ij] = get_joint_ml_entropy(p_joint);
  }
  //Rcout << "done" << std::endl;
    
  // Compute mutual information
  //Rcout << "Computing mutual information...";
  arma::mat mim = arma::mat(n_genes, n_genes, arma::fill::zeros);
  int ij = 0;
  for(int i(0); i<n_genes; ++i)
  {
    for(int j(i); j<n_genes; ++j)
    {
      mim(i,j) = h_marginals[i] + h_marginals[j] - h_joints[ij];
      mim(j,i) = mim(i,j);
      ++ij;
    }
  }
  //Rcout << "done" << std::endl;
  
  return mim;
}


// Mutual information using maximum likelihood entropy estimate and
// Miller-Madow bias correction

// [[Rcpp::export]]
arma::mat mim_MM_cpp(NumericMatrix disc_expr_data, int n_cores)
{
  arma::Mat<int> data = arma::conv_to<arma::Mat<int> >::from(R2armaMat_num(disc_expr_data));
  const int n_genes(data.n_cols), n_samples(data.n_rows);
  const int n_pairs = get_n_gene_pairs(n_genes);

  //Rcout << "There are " << n_genes << " genes and "
        //<< n_samples << " samples" << std::endl;
  
  // Change from R indexing to C++ indexing
  data -= 1;
  
  // Compute marginal entropies
  //Rcout << "Computing marginal entropies...";
  std::vector<double> h_marginals(n_genes);
  for(int j; j<n_genes; ++j)
  {
    // Compute the Miller-Madow correction to the entropy
    arma::vec p_marginal = get_emp_marg_dist(data.col(j));
    int nonzero_bins = arma::size(arma::find(p_marginal))(0);
    double mm_corr = (double)(nonzero_bins - 1)/(2.0 * (double)n_samples);
    
    h_marginals[j] = get_marginal_ml_entropy(p_marginal) + mm_corr;
  }
  //Rcout << "done" << std::endl;
  
  // Compute joint entropies in parallel
  //Rcout << "Computing joint entropies...";
  std::vector<double> h_joints(n_pairs);
  
  // Number of cores to use
  //const int n_cores = sysconf(_SC_NPROCESSORS_ONLN);
  omp_set_num_threads(n_cores);

  const std::vector<std::pair <int,int> > ij_pairs = get_ij_list(n_genes);

  #pragma omp parallel for shared(data, h_joints)
  for(int ij=0; ij<n_pairs; ++ij)
  {
    std::pair <int,int> ij_pair = ij_pairs[ij];
    int i(ij_pair.first), j(ij_pair.second);
    arma::mat p_joint = get_emp_joint_dist(data.col(i), data.col(j));
    int nonzero_bins = arma::size(arma::find(p_joint))(0);
    double mm_corr = (double)(nonzero_bins - 1)/(2.0 * (double)n_samples);
    h_joints[ij] = get_joint_ml_entropy(p_joint) + mm_corr;
  }
  //Rcout << "done" << std::endl;
  
  // Compute mutual information
  //Rcout << "Computing mutual information...";
  arma::mat mim = arma::mat(n_genes, n_genes, arma::fill::zeros);
  int ij(0);
  for(int i(0); i<n_genes; ++i)
  {
    for(int j(i); j<n_genes; ++j)
    {
      mim(i,j) = h_marginals[i] + h_marginals[j] - h_joints[ij];
      mim(j,i) = mim(i,j);
      ++ij;
    }
  }
  //Rcout << "done" << std::endl;

  // Zero negative values
  mim.elem( find(mim < 0.0) ).fill(0.0);
  
  return mim;
}


// Chao-Shen Estimator
// [[Rcpp::export]]
arma::mat mim_CS_cpp(NumericMatrix disc_expr_data, int n_cores)
{
  arma::Mat<int> data = arma::conv_to<arma::Mat<int> >::from(R2armaMat_num(disc_expr_data));
  const int n_genes(data.n_cols), n_samples(data.n_rows);
  const int n_pairs = get_n_gene_pairs(n_genes);

  //Rcout << "There are " << n_genes << " genes and "
        //<< n_samples << " samples" << std::endl;
    
  // Change from R indexing to C++ indexing
  data -= 1;
    
  // Compute marginal entropies
  //Rcout << "Computing marginal entropies...";
  std::vector<double> h_marginals(n_genes);
  for(int j; j<n_genes; ++j)
  {
    // The number of bins with a single count
    arma::vec p_marginal = get_emp_marg_dist(data.col(j));
    int sing_count_bins = arma::size(arma::find(p_marginal==1.0/(double)n_samples))(0);
    double samp_cov = 1.0 - (double)sing_count_bins/(double)n_samples;
    arma::vec cs_corr = 1.0/(1.0 - arma::pow(1.0-samp_cov*p_marginal,n_samples));
      
    // Set inf for empty bins to zero
    cs_corr.elem(arma::find(cs_corr==arma::datum::inf)).zeros();
      
    h_marginals[j] = -arma::sum(samp_cov*p_marginal % arma::log(samp_cov*p_marginal + 1e-16) % cs_corr);
  }
  //Rcout << "done" << std::endl;
    
  // Compute joint entropies in parallel
  //Rcout << "Computing joint entropies...";
  std::vector<double> h_joints(n_pairs);
  
  // Number of cores to use
  //const int n_cores = sysconf(_SC_NPROCESSORS_ONLN);
  omp_set_num_threads(n_cores);
  
  const std::vector<std::pair <int,int> > ij_pairs = get_ij_list(n_genes);

  #pragma omp parallel for shared(data, h_joints)
  for(int ij=0; ij<n_pairs; ++ij)
  {
    std::pair <int,int> ij_pair = ij_pairs[ij];
    int i(ij_pair.first), j(ij_pair.second);
    arma::mat p_joint = get_emp_joint_dist(data.col(i), data.col(j));
    int sing_count_bins = arma::size(arma::find(p_joint==1.0/(double)n_samples))(0);
    double samp_cov = 1.0 - (double)sing_count_bins/(double)n_samples;
    arma::mat cs_corr= 1.0/(1.0 - arma::pow(1.0-samp_cov*p_joint,n_samples));
        
    // Set inf from divsion by zero to zero
    cs_corr.elem(arma::find(cs_corr==arma::datum::inf)).zeros();
        
    arma::mat tmp = samp_cov*p_joint % arma::log(samp_cov*p_joint + 1e-16) % cs_corr;
    h_joints[ij] = -arma::accu(tmp);
  }

  //Rcout << "done" << std::endl;
    
  // Compute mutual information
  //Rcout << "Computing mutual information...";
  arma::mat mim = arma::mat(n_genes, n_genes, arma::fill::zeros);
  int ij(0);
  for(int i(0); i<n_genes; ++i)
  {
    for(int j(i); j<n_genes; ++j)
    {
      mim(i,j) = h_marginals  [i] + h_marginals[j] - h_joints[ij];
      mim(j,i) = mim(i,j);
      ++ij;
    }
  }
  //Rcout << "done" << std::endl;

  return mim;
}


// Shrinkage estimator
// [[Rcpp::export]]
arma::mat mim_shrink_cpp(NumericMatrix disc_expr_data, int n_cores)
{
  arma::Mat<int> data = arma::conv_to<arma::Mat<int> >::from(R2armaMat_num(disc_expr_data));
  const int n_genes(data.n_cols), n_samples(data.n_rows);
  const int n_pairs = get_n_gene_pairs(n_genes);
  
  //Rcout << "There are " << n_genes << " genes and "
        //<< n_samples << " samples" << std::endl;
  
  // Change from R indexing to C++ indexing
  data -= 1;
  
  // Compute marginal entropies
  //Rcout << "Computing marginal entropies...";
  std::vector<double> h_marginals(n_genes);
  for(int j; j<n_genes; ++j)
  {
    // Compute the shrinkage intensity lambda
    arma::vec p_marginal = get_emp_marg_dist(data.col(j));
    //Rcout << "j=" << j << std::endl;
    //Rcout << "p_marg = \n" << p_marginal.t() << std::endl;
    double n_bins = (double)p_marginal.n_elem;
    double lambda_numer = 1.0 - arma::accu(arma::pow(p_marginal, 2.0));
    double lambda_denom = (double)(n_samples-1) *
      arma::accu(arma::pow(1.0/n_bins - p_marginal, 2.0));
    //Rcout << "lambda = " << lambda_num << " \ " << lambda_den << std::endl;
    
    double lambda;
    if(lambda_denom==0)
      lambda = 0;
    else
      lambda = lambda_numer/lambda_denom;

    // Lambda must be between 0 and 1
    if(lambda < 0.0)
      lambda = 0.0;
    else if(lambda > 1)
      lambda = 1.0;

    // Note: the target distribution is uniform
    arma::vec p_marg_shrink = lambda * (1.0/n_bins) + (1.0-lambda) * p_marginal;
    //Rcout << "p_marg_shrink = " << p_marg_shrink.t() << std::endl;
    h_marginals[j] = get_marginal_ml_entropy(p_marg_shrink);
    
  }
  //Rcout << "done" << std::endl;
  
  // Compute joint entropies in parallel
  //Rcout << "Computing joint entropies...";
  std::vector<double> h_joints(n_pairs);
  
  // Number of cores to use  int ij(0);
  //const int n_cores = sysconf(_SC_NPROCESSORS_ONLN);
  omp_set_num_threads(n_cores);

  const std::vector<std::pair <int,int> > ij_pairs = get_ij_list(n_genes);

  #pragma omp parallel for shared(data, h_joints)
  for(int ij=0; ij<n_pairs; ++ij)
  {
    std::pair <int,int> ij_pair = ij_pairs[ij];
    int i(ij_pair.first), j(ij_pair.second);
    arma::mat p_joint = get_emp_joint_dist(data.col(i), data.col(j));
    //Rcout << "p_joint = \n" << p_joint << std::endl;
    double n_bins = (double)p_joint.n_elem;
    double lambda_numer = 1.0 - arma::accu(arma::pow(p_joint, 2.0));
    double lambda_denom = (double)(n_samples-1) *
          arma::accu(arma::pow(1.0/n_bins - p_joint, 2.0));
        
    double lambda;
    if(lambda_denom==0)
      lambda = 0;
    else
      lambda = lambda_numer/lambda_denom;
        
    // Lambda must be between 0 and 1
    if(lambda < 0.0)
      lambda = 0.0;
    else if(lambda > 1)
      lambda = 1.0;
        
    //Rcout << "lambda = " << lambda << std::endl;
        
    // Note: the target distribution is uniform
    arma::mat p_joint_shrink = lambda * (1.0/n_bins) + (1.0-lambda) * p_joint;
    //Rcout << "p_j_shrink = \n" << p_joint_shrink << std::endl;
    h_joints[ij] = get_joint_ml_entropy(p_joint_shrink);
  }
  //Rcout << "done" << std::endl;
  
  // Compute mutual information
  //Rcout << "Computing mutual information...";
  arma::mat mim = arma::mat(n_genes, n_genes, arma::fill::zeros);
  int ij(0);
  for(int i(0); i<n_genes; ++i)
  {
    for(int j(i); j<n_genes; ++j)
    {
      if(false)
      {
        Rcout << "i=" << i << " j=" << j << std::endl;
        Rcout << "h.i = " << h_marginals[i] << " h.j = " << h_marginals[j] << std::endl;
        Rcout << "h.ij = " << h_joints[ij] << std::endl << std::endl;
      }

      mim(i,j) = h_marginals[i] + h_marginals[j] - h_joints[ij];
      mim(j,i) = mim(i,j);
      ++ij;
    }
  }
  //Rcout << "done" << std::endl;
  
  return mim;
}
