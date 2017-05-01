// ----------------------------------------------------------------------------------
//  Implementation of the B-spline Mutual information estimation method
// ----------------------------------------------------------------------------------

#include "fastGeneMI.h"
//#include "util_funcs.cpp"

using namespace Rcpp;

// Speed ideas: memoisation, move knot vector out of class (as it always the same)

// ----------------------------------------------------------------------------------
//  Class that returns bspline coefficients
// ----------------------------------------------------------------------------------

// Constructor
Bspline::Bspline(const int k_, const int n_bins_, const arma::vec& x_)
{
  k = k_;
  n_bins = n_bins_;
  x = x_;
  n_pts = x.n_elem;
  t = get_knots();
  tmax = t.max();
  
  // Transform data to range of knot vector and trim maximum values
  z = ((x-x.min())/(x.max()-x.min())) * tmax;
  z.replace(tmax, tmax-1e-12);
}

void Bspline::print() const
{
  Rcout << "k=" << k << "\t #bins=" << n_bins << "\t#pts" << n_pts << std::endl;
  Rcout << "t^T=" << arma::trans(t) << std::endl;
  Rcout << "x^T=" << arma::trans(x) << std::endl;
  Rcout << "z^T=" << arma::trans(z) << std::endl;
}

// Get the knot vector - take this out of class and pass to constructor
arma::vec Bspline::get_knots() const
{
  int n_knots(k + n_bins);
  arma::vec t(n_knots, arma::fill::zeros);
  
  // Intermediate (non-zero) knots
  for(int i(k); i<n_bins; ++i)
  {
    t(i) = (double)i - (double)k + 1.0;
  }
  
  // End values
  arma::vec ones = arma::vec(k, arma::fill::ones);
  t.rows(n_knots-k, n_knots-1) = ((double)n_bins - (double)k + 1.0) * ones;
  return t;
}

// Return zero-th order basis function
arma::vec Bspline::basis0(const double zi) const
{
  arma::vec res(n_bins+k-1, arma::fill::zeros);
  int loc = arma::find(zi >= t).max();
  res(loc) = 1.0;
  return res;
}

arma::vec Bspline::basis(const double zi, const int p) const
{
  //Rcout << "Calling basis with p=" << p << " at point zi=" << zi << std::endl;
  arma::vec basis_p_minus_1;
  // 0-th order basis function
  if(p==0)
    return basis0(zi);
  
  // Higher orders
  else
    basis_p_minus_1 = basis(zi, p-1);
  
  arma::vec first_term_num = zi - t.rows(0, t.n_elem-p-1);
  arma::vec first_term_denom = t.rows(p, t.n_elem-1) - t.rows(0, t.n_elem-p-1);
  
  arma::vec second_term_num = t.rows(p+1, t.n_elem-1) - zi;
  arma::vec second_term_denom = t.rows(p+1, t.n_elem-1) - t.rows(1, t.n_elem-p-1);
  
  arma::vec first_term = first_term_num/first_term_denom;
  arma::vec second_term = second_term_num/second_term_denom;
  
  // Set inf to zero and remove nan
  first_term.elem(arma::find(arma::abs(first_term)==arma::datum::inf)).zeros();
  second_term.elem(arma::find(arma::abs(second_term)==arma::datum::inf)).zeros();
  first_term.replace(arma::datum::nan, 0.0);
  second_term.replace(arma::datum::nan, 0.0);
  
  return first_term.rows(0, first_term.n_elem-2) %
  basis_p_minus_1.rows(0, basis_p_minus_1.n_elem-2) +
  second_term % basis_p_minus_1.rows(1, basis_p_minus_1.n_elem-1);
}

// Evaluate the B-spline basis functions at each point and return
// the result as a sparse matrix
arma::sp_mat Bspline::compute_bspl_coeffs()
{
  int ct(0);
  arma::mat res = arma::mat(n_pts, n_bins, arma::fill::zeros);
  for(int i(0); i<n_pts; ++i)
  {
    res.row(i) = basis(z[i], k-1).t();
  }
  // Zero small values
  res.elem(arma::find(res<1e-10)).zeros();
  return arma::sp_mat(res);
}

// Wrapper for computing bspline coefficients
arma::sp_mat get_bspl_coeffs(const arma::vec& x, const int n_bins, const int k)
{
  Bspline bspl = Bspline(k, n_bins, x);
  //bspl.print();
  return bspl.compute_bspl_coeffs();
}


// Compute the mutual information matrix of expression data using the B-spline method using 
// splines of order k
// [[Rcpp::export]]
arma::mat mim_bspline_cpp(NumericMatrix expr_data, const int k, const int n_bins, const int n_cores)
{
  // Convert arguments to Armadillo objects
  const arma::mat data = R2armaMat_num(expr_data); // Continuous data does not need reindexing R->C++
  const int n_genes(data.n_cols), n_samples(data.n_rows);
  const int n_pairs = (int)(((double)n_genes/2.0) * (double)(n_genes+1));
  
  // Compute bspline coefficients in parallel
  int j;
  
  // Number of cores to use
  //const int n_cores = sysconf(_SC_NPROCESSORS_ONLN);
  omp_set_num_threads(n_cores);
  
  std::vector<arma::sp_mat> bspl_coeffs(n_genes);
  #pragma omp parallel for shared(bspl_coeffs) private(j) schedule(auto) default(none)
  for(j=0; j<n_genes; ++j)
  {
    bspl_coeffs[j] = get_bspl_coeffs(data.col(j), n_bins, k);
  }
  
  //Rcout << "Got bspline coefficients\n";
  
  // Compute marginal probabilities and entropies
  arma::vec p_marginal;
  std::vector<double> h_marginals(n_genes);
  for(int j(0); j<n_genes; ++j)
  {
    p_marginal = arma::vec(arma::trans(arma::sum(bspl_coeffs[j], 0)/(double)n_samples));
    h_marginals[j] = -arma::sum(p_marginal % arma::log(p_marginal + 1e-16));
  }
  
  //Rcout << "Got marginal probabilities and entropies\n";
  
  // Compute joint probabilities and entropies in parallel
  arma::mat p_joint;
  const arma::Mat<int> idx_lookup = get_idx_lookup_mat(n_genes);
  std::vector<double> h_joints(n_pairs);
  int i;
  
  #pragma omp parallel for shared(h_joints,bspl_coeffs) private(i,j,p_joint) schedule(dynamic) default(none) collapse(2)
  for(i=0; i<n_genes; ++i)
  {
    for(j=0; j<n_genes; ++j)
    {
      if(i<=j)
      {
        int idx = idx_lookup(i,j);
        p_joint = (bspl_coeffs[i].t() * bspl_coeffs[j])/(double)n_samples;
        h_joints[idx] = -arma::accu(p_joint % arma::log(p_joint + 1e-16));
      }
    }
  }
  
  //Rcout << "Got joint probabilities and entropies\n";
  
  // Compute mutual information matrix
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
  
  return mim;
  
}

