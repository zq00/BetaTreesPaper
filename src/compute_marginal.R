# 1. Estimate marginal distributions 
# INPUT
# X - data matrix
# var_names - variables for which we construct the histogram 
# lower, upper - a named vector containing lower and upper bounds of the data. 
#     This will be the range of data values +/- 1 if not provided 
estimate_marginal <- function(X, var_names = colnames(X), lower = NULL, upper = NULL){
  Y <- select(X, all_of(var_names))
  lower_bound <- numeric(length(var_names))
  upper_bound <- numeric(length(var_names))
  if(is.null(lower)){
    lower <- apply(Y, 2, min) - 1
    names(lower) <- var_names
  }
  if(is.null(upper)){
    upper <- apply(Y, 2, max) + 1
    names(upper) <- var_names
  }
  # estimate marginal CDF 
  est_density <- list()
  est_cdf <- list()
  est_icdf <- list()
  
  for(var in var_names){
    x <- pull(Y, var)
    est_density[[var]] <- density(x, n = 2^12)
    est_cdf[[var]] <- CDF(est_density[[var]])
    est_icdf[[var]] <- inverse(est_cdf[[var]], lower = lower[var], upper = upper[var])
  }
  
  # compute the lower and the upper end 
  X_transformed <- matrix(NA, nrow = nrow(Y), ncol = ncol(Y))
  colnames(X_transformed) <- colnames(Y)
  for(var in var_names){
    X_transformed[,var] <- est_cdf[[var]](pull(Y, var))
  }
  
  X_transformed <- as_tibble(X_transformed)
  return(
    list(X_transformed = X_transformed,
         est_density = est_density,
         est_cdf = est_cdf,
         est_icdf = est_icdf,
         lower_bound = lower_bound,
         upper_bound = upper_bound)
  )
}


# 2. Compute original histogram from the transformed histogram 
# INPUT
# X - data matrix
# hist_transformed - adaptive histogram for tilde(X) = F_hat(X)
# var_names - variables for which we construct the histogram 
# est_icdf - estimated inverse CDF 
# lower, upper - a named vector containing lower and upper bounds of the data. 
#     This will be the range of data values +/- 1 if not provided 
compute_histogram_original <- function(X, 
                                       hist_transformed, var_names, 
                                       est_icdf, 
                                       lower = NULL, upper = NULL){
  d <- length(var_names)
  
  Y <- select(X, all_of(var_names))
  lower_bound <- numeric(length(var_names))
  upper_bound <- numeric(length(var_names))
  
  if(is.null(lower)){
    lower <- apply(Y, 2, min) - 1
    names(lower) <- var_names
  }
  if(is.null(upper)){
    upper <- apply(Y, 2, max) + 1
    names(upper) <- var_names
  }
  
  colnames(hist_transformed) <- c(sapply(var_names, function(t) paste0(t, "_l")), 
                                  sapply(var_names, function(t) paste0(t, "_u")),
                                  "density", "lower_density", "upper_density", "n", "level")
  
  hist_original <- hist_transformed
  for(var in var_names){
    hist_original[,paste0(var, "_l")] <- sapply(hist_transformed[,paste0(var, "_l")], 
                                                function(t) {
                                                  if(t == 0) {return(lower[var])}
                                                  else if(t==1) {return(upper[var])} 
                                                  else {return(est_icdf[[var]](t))}
                                                })
    
    hist_original[,paste0(var, "_u")] <- sapply(hist_transformed[,paste0(var, "_u")], 
                                                function(t) { 
                                                  if(t == 0) {return(lower[var])}
                                                  else if(t==1) {return(upper[var])} 
                                                  else {return(est_icdf[[var]](t))}
                                                })
  }
  
  volume_transformed <- apply(hist_transformed[,(d+1):(2*d)] - hist_transformed[,1:d], 1, prod) 
  volume_original <- apply(hist_original[,(d+1):(2*d)] - hist_original[,1:d], 1, prod) 
  
  hist_original[,(2*d+1)] <- volume_transformed * hist_transformed[,(2*d+1)] / volume_original
  hist_original[,(2*d+2)] <- volume_transformed * hist_transformed[,(2*d+2)] / volume_original
  hist_original[,(2*d+3)] <- volume_transformed * hist_transformed[,(2*d+3)] / volume_original
  
  return(hist_original)
}








