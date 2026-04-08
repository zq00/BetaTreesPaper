# Functions to compute summary statistics of the histogram

# 1. Compute the true probability mass in each region  
# INPUT
# rect - a matrix with 2d columns, the first and last d columns correspond to the lower and upper end of histogram regions 
# distribution - "gaussian_mixture" or "hd"
# ... - additional arguments for the Gaussian mixture distribution 
# OUTPUT
# probability mass in each region in the rect 
compute_probability_mass <- function(rect, distribution, ...){
  d <- ncol(rect) / 2
  if(distribution == "gaussian_mixture"){
    true_log_prob <- matrix(NA, nrow = nrow(rect), ncol = nrow(mu))
    
    for(j in 1:nrow(mu)){
      coord_j <- coord[j, ]
      
      if(d > ncol(mu)){
        coord_mj <- (1:d)[-coord[j, ]]
      }
      
      for(i in 1:nrow(rect)){
        val <- pmvnorm(lower = rect[i, coord_j], upper = rect[i, (coord_j + d)],
                       mean = mu[j, coord_j], sigma = sigma[,,j])[1]
        if(d > ncol(mu)){
          true_log_prob[i, j] <- ifelse(val <= 0, -Inf, log(val)) + 
            sum(log(pnorm(rect[i, (coord_mj + d)], mu[j, coord_mj], sd = 1) - pnorm(rect[i, coord_mj], mu[j, coord_mj], sd = 1)))
        }else{
          true_log_prob[i, j] <- ifelse(val <= 0, -Inf, log(val))
        }
      }
    }
    
    return(colSums(t(exp(true_log_prob)) * prob, na.rm = T))
  }
  if(distribution == "hd"){
    
    mu_1 <- c(-1.5, 0.6)
    mu_2 <- c(2, -1.5)
    sigma_1 <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, byrow = T)
    sigma_2 <- matrix(c(1, 0, 0, 1), 2, 2)
    
    return(apply(rect, 1, function(t){
      # first two dimensions
      p1 <- pmvnorm(lower = t[1:2], upper = t[(d+1):(d+2)], mean = mu_1, sigma = sigma_1)[1] * 0.4 + 
        pmvnorm(lower = t[1:2], upper = t[(d+1):(d+2)], mean = mu_2, sigma = sigma_2)[1] * 0.6
      # last dimensions
      p2 <- 0.5 * (pnorm(t[(d+3):(2*d)], mean = 3.5, sd = 0.5) - pnorm(t[3:d], mean = 3.5, sd = 0.5)) +
        0.5 * (pnorm(t[(d+3):(2*d)],mean = 6, sd = 0.25) - pnorm(t[3:d], mean = 6, sd = 0.25))
      p1 * prod(p2)
    }))
  }
}


# 2. Plot the average of true density within each region versus histogram estimates 
# INPUT
# true_probs - true probability within each region 
# est_avg_density - estimated average density within each region 
# lower_bound, upper_bound - lower and upper bound of average density within each region
# rect - a matrix with 2d columns, the first and last d columns correspond to the lower and upper end of histogram regions 
# Note: (rect, est_avg_density, lower_bound, upper_bound) are the first 2d + 3 columns of the Beta-tree output 
# OUTPUT
# a grapph of the estimated average density within each region versus true average density 
plot_density <- function(true_probs, est_avg_density, lower_bound, upper_bound, rect){
  nregion <- length(true_probs)
  
  d <- ncol(rect) / 2
  volume <- apply(rect[,(d+1):(2*d)] - rect[,1:d], 1, prod) 
  
  
  lower <- apply(rect[, 1:d], 2, min)
  upper <- apply(rect[, (d+1):(2*d)], 2, max)
  
  volume_standardized <- volume / prod(upper - lower) # divide volume by the total 
  
  true_standardized <- true_probs / volume_standardized
  est_standardized <- est_avg_density * volume / volume_standardized
  lower_standardized <- lower_bound * volume / volume_standardized
  upper_standardized <- upper_bound * volume / volume_standardized
  
  order <- sort(true_standardized, index.return = T)$ix
  ggplot() + 
    geom_segment(aes(x = 1:nregion, xend = 1:nregion, 
                     y = log(lower_standardized[order]), yend = log(upper_standardized[order])), color = "grey80") + 
    geom_point(aes(x = 1:nregion, y = log(true_standardized[order])), size = 0.2) + 
    xlab("Region index") + 
    ylab("Average density (log)") + 
    theme_classic() + 
    theme(text = element_text(size = 18))
}

# 3. Compute width of CI
# INPUTS
# d - data dimension
# B - number of repetitions
# file_path - location to store intermediate results 
# ... - additional inputs
get_summary_stat <- function(d, B = 20, file_path, ... ){
  
  data <- tibble(.rows = 0)
  
  n_vals <- get("n_val")
  alpha_vals <- get("alpha_val")
  
  for(i in 1:length(n_vals)){
    for(j in 1:length(alpha_vals)){
      # cat("\n n = ", n_vals[i], ", alpha = ", alpha_vals[j], ":\n")
      
      nbox <- numeric(B)
      ci_width <- matrix(NA, nrow = B, ncol = 3)
      lhs <- numeric(B)
      for(b in 1:B){
        X <- sample_gaussian(n_vals[i], d, Sigma)
        
        hist <- BuildHist(X, alpha = alpha_vals[j], 
                          method = "weighted_bonferroni", plot = FALSE,
                          bounded = TRUE, option = "ndat", qt = rep(1, d)) 
       
        hist <- as.matrix(hist, ncol = (2*d + 5))
        nbox[b] <- nrow(hist)
        # calculate summary statistics statistics 
        volume <- apply(hist[,(d+1):(2*d),drop = F] - hist[,1:d,drop = F], 1, prod)
        Fn <- (hist[,(2*d+ 4)] + 1 ) / n_vals[i]
        G_up <- hist[,(2*d+ 3)] * volume
        G_low <- hist[,(2*d+ 3)] * volume
        lhs_up <- sqrt(n_vals[i]) * abs(G_up - Fn) / sqrt(G_up * (1-G_up))
        lhs_low <- sqrt(n_vals[i]) * abs(G_low - Fn) / sqrt(G_low * (1-G_low))
        
        ci_width[b, ] <- quantile(hist[,(2*d + 3)]* volume - hist[,(2*d + 2)]* volume, c(0.25, 0.5, 0.75)) 
        lhs[b] <- max(c(lhs_up, lhs_low))
      }
      # store data to a table 
      if(!is.null(file_path)){
        write.table(cbind(nbox, ci_width), 
                    file = paste0(file_path, method, "_", d, "_", n_vals[i], "_", alpha_vals[j], ".txt"))
        
      }
  
      new_tab <- tibble(n = n_vals[i], alpha = alpha_vals[j],
                        ci_width_25 = mean(ci_width[,1]),
                        ci_width_25_sd = sd(ci_width[,1]),
                        ci_width_50 = mean(ci_width[,2]),
                        ci_width_50_sd = sd(ci_width[,2]),
                        ci_width_75 = mean(ci_width[,3]),
                        ci_width_75_sd = sd(ci_width[,3]),
                        nbox = mean(nbox),
                        lhs = mean(lhs))
      data <- rbind(data, new_tab)
    }
  }
  
  return(data)
}


# 4. Evaluate the mode hunting algorithm 
# INPUTS 
# n - sample size
# d - data dimension 
# m - number of clusters 
# method - "exact" or "approximate"
# nrep - number of repetitions
# ndim - dimension of the mixture components, 4 by default 
# ... - additional parameters: 
#   approximate algorithm: L - length of path; B - number of paths
#   exact algorithm: cutoff - cutoff if path length 
# OUTPUTS
# alg_time - computation time
# alg_distance - a list of distance matrix of each iteration
summary_mode_finding <- function(n, d, m, method, nrep, ndim = 4, ...){
  sigma <- array(NA, dim = c(ndim, ndim, m))
  for(i in 1:m){
    sigma[,,i] <- diag(rep(1, ndim))
  }
  
  coord <- t(replicate(m, sample(1:d, size = ndim, replace = F))) 
  prob <- rep(1/m, m)
  
  curr_vertices <- expand_grid(c(-1, 1), c(-1, 1))
  colnames(curr_vertices) <- c("X1", "X2")
  for(i in 2:(d-1)){
    curr_vertices <- expand_grid(curr_vertices,  c(-1, 1))
    colnames(curr_vertices) <- sapply(1:(i+1), function(t) paste0("X", t))
  }
  which_vertices <- sample(1:(2^d), m, replace = F)
  
  mu <- matrix(0, nrow = m, ncol = d)
  for(i in 1:m){
    mu[i,  coord[i, ]] <-  curr_vertices[which_vertices[i], coord[i, ]] %>% unlist()
  }
  # make mu larger
  mu <- mu * 8
  
  cat("number of duplicate distince cluster centers = ", sum(duplicated(mu)), "\n")
  while(sum(duplicated(mu)) > 0){
    ind_duplicate <- which(duplicated(mu) == T)
    which_vertices <- sample(1:(2^d), length(ind_duplicate), replace = F)
    for(i in 1:length(ind_duplicate))
    mu[ind_duplicate[i], ] <- curr_vertices[which_vertices[i], coord[ind_duplicate[i], ]] %>% unlist()
  }
  
  alg_time <- numeric(nrep)
  alg_distance <- list()
  n_region <- numeric(nrep)
  
  for(i in 1:nrep){
    cat(i, ":")
    X <- sample_gaussian_mixture(n, d, mu, sigma, coord, prob)
    hist_adaptive <- build_adaptive_histogram(as.matrix(X), thresh_marginal = 0.1 / d,
                                              thresh_interaction = qgamma(shape = d - 1,rate = 1, p = 1 - 0.1), 
                                              alpha = 0.1, 
                                              method = "weighted_bonferroni")$hist
    n_region[i] <- nrow(hist_adaptive)
    
    start_time <- Sys.time()
    if(method == "exact"){
      args <- list(...)
      if("cutoff" %in% names(args)){
        est_modes <- FindModes(hist_adaptive, d, cutoff = args$cutoff)
      }else{
        est_modes <- FindModes(hist_adaptive, d)
      }
      
    }else if(method == "approximate"){
      args <- list(...)
      est_modes <- FindModesApproximate(hist = hist_adaptive, d = d, L = args$L, B = args$B)
    }
    end_time <- Sys.time()
    alg_time[i] <- difftime(end_time, start_time, units = "secs")
    
    distance <- compute_distance(est_modes$mode, hist_adaptive, mu)
    alg_distance[[i]] <- distance
    
    cat(alg_time[i], "\n")
  }
  
  return(
    list(alg_time = alg_time,
         alg_distance = alg_distance,
         mu = mu,
         n_region = n_region)
  )
}

# 5. Compute distance between the estimated modes and the center of the Gaussian mixture 
# INPUTS
# est_modes - indices of estimated modes
# est_hist - estimated histogram
# mu - true mode
# OUTPUTS
# A matrix of distance between the estimated mode and true mode 
# each row correspond to one true mode
compute_distance <- function(est_modes, est_hist, mu){
  ncluster <- nrow(mu)
  d <- ncol(mu)
  nmodes <- length(est_modes)
  
  distance <- matrix(0, nrow = ncluster, ncol = nmodes)
  for(i in 1:ncluster){
    for(j in 1:nmodes){
      # calculate the point in the region with the smallest distance to the true mode 
      loc <- pmax(pmin(mu[i, ], est_hist[est_modes[j], (d+1):(2*d)]),  est_hist[est_modes[j], 1:d])
      
      distance[i, j] <- sqrt(sum((mu[i, ]  - loc)^2))
    }
  }
  distance
}












