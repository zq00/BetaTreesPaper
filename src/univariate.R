# Functions for plotting confidence intervals for two univariate distributions

# sample data
sample_gaussian_1d <- function(n){
  matrix(rnorm(n, 0, 1), nrow = n, ncol = 1)
}

sample_gaussian <- function(n, d, sigma = diag(d)){
  R <- chol(sigma)
  
  matrix(rnorm(n * d, 0, 1), nrow = n, ncol = d) %*% R
}

sample_harp <- function(n){
  class <- rmultinom(1, n, prob = rep(0.2, 5))
  
  X <- c(rnorm(class[1,1],0, 0.5), rnorm(class[2,1],5, 1), rnorm(class[3,1],15, 2), 
         rnorm(class[4,1],30, 4), rnorm(class[5,1],60, 8))
  matrix(X, nrow = n, ncol = 1)
}

# compute true densities
density_gaussian_1d <- tibble(
  t = seq(-2.5, 2.5, length.out = 500), 
  density = dnorm(t, 0, 1)
)

density_harp <- tibble(
  t = seq(0, 85, length.out = 10000), 
  density = 0.2 * dnorm(t,0, 0.5) + 0.2 * dnorm(t,5, 1) + 0.2 * dnorm(t,15, 2) +
    0.2 * dnorm(t,30, 4) + 0.2 * dnorm(t,60, 8)
)

# plot univariate histogram
plot_univariate <- function(n, distribution, alpha){
  if(distribution == "gaussian"){
    X <- sample_gaussian_1d(n)
  }else if(distribution == "harp"){
    X <- sample_harp(n)
  }
  
  hist <- BuildHist(X, alpha = alpha, method = "weighted_bonferroni", plot = F) 
  colnames(hist) <- c("left", "right", "density", "lower", "upper", "nobs", "depth")
  hist <- hist %>% 
    as_tibble() %>% 
    mutate(ci_width = upper - lower,
           width = right - left)
  
  g_hist <- ggplot(hist) + 
    geom_rect(aes(xmin=left, xmax = right, 
                  ymin = 0, ymax = density), 
              fill = "grey80", color = "white") + 
    geom_segment(aes(x = (left + right) / 2, 
                     xend = (left + right)/2, 
                     y = lower, yend = upper),color="grey20",
                 arrow = arrow(angle = 90, ends = "both", 
                               length = unit(0.03, "inches"))) + 
    
    xlab("X") + 
    ylab("Density") + 
    theme_classic() + 
    theme(text = element_text(size = 18)) 
  
  if(distribution == "gaussian"){
    g_hist <- g_hist + geom_line(aes(x = t, y = density), data = density_gaussian_1d, color = "black", linewidth = 0.6) 
  }else if(distribution == "harp"){
    g_hist <-  g_hist + geom_line(aes(x = t, y = density), data = density_harp, color = "black", linewidth = 0.6)  
  }
  
  list(hist = hist, g_hist = g_hist)
}
