# functions to sample observations 

# 1. Mixture of correlated Gaussian 
# INPUT 
# sigma - an array of covariance matrix of each cluster 
# mu - a matrix of the mean of each cluster 
# coord - a matrix of coordinates corresponding to each cluster 
# n - number of observations
# d - dimensions
# prob - vector of probability to be in each cluster 
# the rest of the coordinates are i.i.d. N(0, 1)
# OUTPUT
# a matrix with n rows and d columns, each row corresponding to one observation
sample_gaussian_mixture <- function(n, d, mu, sigma, coord, prob){
  # cluster assignment 
  class <- rmultinom(1, size = n, prob = prob)
  class_sum <- c(0, cumsum(class))
  # draw from this Gaussian distribution
  X <- matrix(rnorm(n*d, 0, 1), n, d)
  for(i in 1:nrow(mu)){
    # large covariance matrix
    sigma_large <- diag(1, nrow = d, ncol = d)
    sigma_large[coord[i, ], coord[i,]] <- sigma[, ,i]
    R <- chol(sigma_large)
    X[(class_sum[i] + 1):class_sum[i+1], ] <- X[(class_sum[i] + 1):class_sum[i+1], ] %*% R 
    X[(class_sum[i] + 1) :class_sum[i+1], ] <- t(t(X[(class_sum[i] + 1):class_sum[i+1], ] ) + as.vector(mu[i,]))
  }
  # assign column names
  colnames(X) <- sapply(1:d, function(x) paste0("X", x))
  X <- as.data.frame(X)
  
  return(X)
}

# 2. High-dimensional example 
sample_hd <- function(n, d){
  X <- matrix(NA, n, d)
  
  sigma_1 <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, byrow = T)
  sigma_2 <- matrix(c(1, 0, 0, 1), 2, 2)
  
  mu_1 <- c(-1.5, 0.6)
  mu_2 <- c(2, -1.5)
  
  R1 <- chol(sigma_1)
  R2 <- chol(sigma_2)
  
  # first two dimensions 
  x1 <- matrix(rnorm(n*2, 0, 1), n, 2)
  X[,1:2] <- t( apply(x1, 1, function(t) if(sample(1:2, 1, prob = c(0.4, 0.6)) == 1){ t%*% R1 + mu_1
  }else{ t%*% R2 + mu_2}))
  
  # rest of the dimensions
  classes <- matrix(sample(1:2, n*(d-2),replace = T), n, d-2)
  X[,-c(1:2)] <- matrix(rnorm(n*(d-2), 0, 1), n, d-2) / classes / 2 + 1 + 2.5 * classes
  
  # assign column names
  colnames(X) <- sapply(1:d, function(x) paste0("X", x))
  X <- as.data.frame(X)
}















