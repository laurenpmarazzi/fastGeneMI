// ----------------------------------------------------------------------------------
//  Implementations of the Mutual Information estimators using Maximum likelihood,
//  Miller-Madow, Chao-Shen, Jack-knififed (unvalidated) and Shrinkage entropy 
//  estimation (problem in joint entropy computation)
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
  const int n_pairs = (int)(((double)n_genes/2.0) * (double)(n_genes+1));

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
  const arma::Mat<int> idx_lookup = get_idx_lookup_mat(n_genes);
  arma::mat p_joint;
  std::vector<double> h_joints(n_pairs); // Change to armadillo vector for returning to R
  int i, j;
  
  // Number of cores to use
  // const int n_cores = sysconf(_SC_NPROCESSORS_ONLN);
  omp_set_num_threads(n_cores);
  
  //const int n_threads = omp_get_num_threads();
  //Rcout << "There are " << n_threads << " threads\n";
  
  #pragma omp parallel for shared(data, h_joints) private(i,j,p_joint) schedule(dynamic) default(none) collapse(2)
  for(i=0; i<n_genes; ++i)
  {
    for(j=0; j<n_genes; ++j)
    {
      if(i<=j)
      {
        int idx = idx_lookup(i,j);
        p_joint = get_emp_joint_dist(data.col(i), data.col(j));
        h_joints[idx] = get_joint_ml_entropy(p_joint);
      }
    }
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
  
  return mim;
}


// Mutual information using maximum likelihood entropy estimate and
// Miller-Madow bias correction

// [[Rcpp::export]]
arma::mat mim_MM_cpp(NumericMatrix disc_expr_data, int n_cores)
{
  arma::Mat<int> data = arma::conv_to<arma::Mat<int> >::from(R2armaMat_num(disc_expr_data));
  const int n_genes(data.n_cols), n_samples(data.n_rows);
  const int n_pairs = (int)(((double)n_genes/2.0) * (double)(n_genes+1));

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
  arma::mat p_joint;
  const arma::Mat<int> idx_lookup = get_idx_lookup_mat(n_genes);
  std::vector<double> h_joints(n_pairs);
  int i, j;
  
  // Number of cores to use
  //const int n_cores = sysconf(_SC_NPROCESSORS_ONLN);
  omp_set_num_threads(n_cores);
  
  #pragma omp parallel for shared(data,h_joints) private(i,j,p_joint) schedule(dynamic) default(none) collapse(2)
  for(i=0; i<n_genes; ++i)
  {
    for(j=0; j<n_genes; ++j)
    {
      if(i<=j)
      {
        // Compute the Miller-Madow correction to the entropy
        int idx = idx_lookup(i,j);
        p_joint = get_emp_joint_dist(data.col(i), data.col(j));
        int nonzero_bins = arma::size(arma::find(p_joint))(0);
        double mm_corr = (double)(nonzero_bins - 1)/(2.0 * (double)n_samples);
        
        h_joints[idx] = get_joint_ml_entropy(p_joint) + mm_corr;
      }
      
    }
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
  
  return mim;
}


// Chao-Shen Estimator
// [[Rcpp::export]]
arma::mat mim_CS_cpp(NumericMatrix disc_expr_data, int n_cores)
{
  arma::Mat<int> data = arma::conv_to<arma::Mat<int> >::from(R2armaMat_num(disc_expr_data));
  const int n_genes(data.n_cols), n_samples(data.n_rows);
  const int n_pairs = (int)(((double)n_genes/2.0) * (double)(n_genes+1));

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
  const arma::Mat<int> idx_lookup = get_idx_lookup_mat(n_genes);
  arma::mat p_joint;
  std::vector<double> h_joints(n_pairs);
  int i, j;
  
  // Number of cores to use
  //const int n_cores = sysconf(_SC_NPROCESSORS_ONLN);
  omp_set_num_threads(n_cores);
  
  #pragma omp parallel for shared(data,h_joints) private(i,j,p_joint) schedule(dynamic) default(none) collapse(2)
  for(i=0; i<n_genes; ++i)
  {
    for(j=0; j<n_genes; ++j)
    {
      if(i<=j)
      {
        // Get the number of bins with a single count
        int idx = idx_lookup(i,j);
        p_joint = get_emp_joint_dist(data.col(i), data.col(j));
        int sing_count_bins = arma::size(arma::find(p_joint==1.0/(double)n_samples))(0);
        double samp_cov = 1.0 - (double)sing_count_bins/(double)n_samples;
        arma::mat cs_corr= 1.0/(1.0 - arma::pow(1.0-samp_cov*p_joint,n_samples));
        
        // Set inf from divsion by zero to zero
        cs_corr.elem(arma::find(cs_corr==arma::datum::inf)).zeros();
        
        arma::mat tmp = samp_cov*p_joint % arma::log(samp_cov*p_joint + 1e-16) % cs_corr;
        h_joints[idx] = -arma::accu(tmp);
      }
      
    }
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
    
  return mim;
}


// Shrinkage estimator
// [[Rcpp::export]]
arma::mat mim_shrink_cpp(NumericMatrix disc_expr_data, int n_cores)
{
  arma::Mat<int> data = arma::conv_to<arma::Mat<int> >::from(R2armaMat_num(disc_expr_data));
  const int n_genes(data.n_cols), n_samples(data.n_rows);
  const int n_pairs = (int)(((double)n_genes/2.0) * (double)(n_genes+1));
  
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
  arma::mat p_joint, p_joint_shrink;
  int i, j;
  const arma::Mat<int> idx_lookup = get_idx_lookup_mat(n_genes);
  
  // Number of cores to use  int ij(0);
  //const int n_cores = sysconf(_SC_NPROCESSORS_ONLN);
  omp_set_num_threads(n_cores);
  
  #pragma omp parallel for shared(data,h_joints) private(i,j,p_joint,p_joint_shrink) schedule(dynamic) default(none) collapse(2)
  for(i=0; i<n_genes; ++i)
  {
    for(j=0; j<n_genes; ++j)
    {
      // Rcout << "i = " << i << " j = " << j << std::endl;
      if(i<=j)
      {
        // Compute the shrinkage intensity lambda
        int idx = idx_lookup(i,j);
        p_joint = get_emp_joint_dist(data.col(i), data.col(j));
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
        p_joint_shrink = lambda * (1.0/n_bins) + (1.0-lambda) * p_joint;
        //Rcout << "p_j_shrink = \n" << p_joint_shrink << std::endl;
        h_joints[idx] = get_joint_ml_entropy(p_joint_shrink);
      }
    }
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
