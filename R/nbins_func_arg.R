# Called when the nbins argument to a discrete estimator is a function. Function
# must take the data column as an argument and return a integer number of bins

get.disc.data <- function(expr.data, discretisation, n.bins) {
  # Default number of bins is (int)n_samples^(1/3) if not using 
  # Bayesian Blocks
  if(discretisation=="bb")
    disc.data <- disc_dataset_bb_cpp(expr.data)
  
  # Equal width or equal frequency bins
  else {
    if(as.integer(n.bins)==n.bins)
      disc.data <- infotheo::discretize(expr.data, disc = discretisation,
                                        nbins = n.bins)
    else
      stop(paste(discretisation,
                 "discretisation requires an integer number of bins"))
  }
  
  
  return(disc.data)
}

