# functions to compute the adaptive Beta-tree histogram 

# 1. Compute partition dimension 
# INPUT
# x - observations inside a region
# thresh_marginal, thresh_interaction - threshold on the p-value (marginal) and sum(-log(p-value)) (interaction)
# lower, upper - lower and upper bounds of a region
# OUTPUT 
# below_marginal - coordinates rejected in marginal test 
# below_interaction - coordinates rejected in interaction test (NA if some coordinates are rejected in the marginal test)
# p - coordinate for partitioning 
get_partition_dim <- function(x, thresh_marginal,  thresh_interaction, lower, upper){ 
  d <- ncol(x)
  
  # marginal test
  p_marginal <- numeric(d)
  for(i in 1:d){
    p_marginal[i] <- ad.test(x[,i], null = "punif", min = lower[i], max = upper[i])$p
  }
  # if the smallest p value is less than thresh_marginal, sample one coordinates from rejections
  if(min(p_marginal) < thresh_marginal){
    p <- sample(which(p_marginal < thresh_marginal), size = 1)
    return(
      list(
        below_marginal = which(p_marginal < thresh_marginal),
        below_interaction = NA,
        p = p)
    )
  }
  
  # test pairwise interactions if no rejection
  mlog_p_interatction <- matrix(NA, d, d)
  for(i in 1:(d-1)){
    for(j in (i+1):d){
      # calculate number of obs in each quadrant
      mid <- (lower[c(i, j)] + upper[c(i, j)]) / 2
      n00 <- sum(x[,i] < mid[1] & x[,j] > mid[2]) # upper left
      n01 <- sum(x[,i] > mid[1] & x[,j] > mid[2])   # upper right
      n11 <- sum(x[,i] > mid[1] & x[,j] < mid[2])  # lower right
      n10 <- sum(x[,i] < mid[1] & x[,j] < mid[2])  # lower left 
      dat <- matrix(c(n00, n01, n10, n11), 2, 2, byrow = T)
      mlog_p_interatction[i, j] <- mlog_p_interatction[j, i] <-  -log(fisher.test(dat)$p)
    }
  }
  # combine interaction p-values
  mlog_p_interaction_combined <- rowSums(mlog_p_interatction, na.rm = T)
  
  if(max(mlog_p_interaction_combined) > thresh_interaction){
    p <- (sample(which(mlog_p_interaction_combined > thresh_interaction), size = 1))
    
    return(list(below_marginal = NA, 
                below_interaction = which(mlog_p_interaction_combined > thresh_interaction),
                p = p)
    )
  }else{ # if cannot detect interaction, randomly pick a coordinate and split along the median
    p <- sample(1:d, 1)
    return(
      list(below_marginal = NA, 
           below_interaction = NA,
           p = p
      )
    )
  }
}

# 2. Adaptive Beta-tree histogram 
# the adaptive Beta-tree histogram uses a bounded histogram that removes the observations at the tail of each direction by default 
# INPUT 
# X - data 
# thresh_marginal, thresh_interaction - threshold on the p-value (marginal) and sum(-log(p-value)) (interaction)
# alpha - significance level
# method - method to adjust for multiple testing, "bonferroni" or "weighted_bonferroni"

build_adaptive_histogram <- function(X, thresh_marginal, thresh_interaction, 
                                     alpha = 0.1, method = "weighted_bonferroni"){
  nd <- NA # No. regions in each level
  split_dir <- NA # record the dimension used for partitioning at each step  
  d <- ncol(X)
  n <- nrow(X)
  
  # use the extreme obs. as bounds 
  lower <- numeric(d)
  upper <- numeric(d)
  for(i in 1:d){
    Xi <- X[,i]
    ind_min <- which.min(Xi)
    ind_max <- which.max(Xi)
    # set bounds
    lower[i] <- Xi[ind_min]
    upper[i] <- Xi[ind_max]
    # remove obs.
    X <- X[-c(ind_min, ind_max), ]
  }
  
  # put the bounds here
  rootnode <- list(
    leftchild = NULL,
    rightchild = NULL,
    ndat = nrow(X),   # number of observ. in this node
    depth = 0,
    low = lower,   # vector of lower bound on region
    up  =  upper,
    lower = NA,  # lower confidence bound for average density
    upper = NA,  # upper confidence bound
    leaf  = FALSE, # node is leaf or not
    bounded = TRUE, 
    dir = NA
  )
  
  # add a new node in the k-d tree
  add_node <- function(x, node, d, thresh_marginal, thresh_interaction){
    
    if (is.na(nd[node$depth + 1])) {
      nd[node$depth + 1] <<- 1L
    } else {
      nd[node$depth + 1] <<- nd[node$depth + 1] + 1L
    }
    
    if (node$ndat < 4*log(n)) { node$leaf = TRUE   # node is leaf, return it marked as such
    } else {  # split this node:
      leftnode = node; rightnode = node   # initialize by taking info from parent
      
      # partition dimension
      dir <- get_partition_dim(x, lower = node$low, upper = node$up, 
                               thresh_marginal = thresh_marginal, 
                               thresh_interaction = thresh_interaction)
      node$dir <- dir$p
      p <- dir$p
      qt <- 0.5 
      # debug: print partition direction
      # cat(dir$below_marginal, "   ;", dir$below_interaction, "\n")
      split_dir <<- c(split_dir, dir$p)
      
      depth = node$depth + 1
      leftnode$depth = depth; rightnode$depth = depth
      
      x = x[order(x[,p]),, drop = F] # sort according to the pth coordinate
      m = node$ndat  # number of data in x
      
      leftnode$ndat = ceiling(qt * m) - 1; rightnode$ndat = m-ceiling(qt * m) 
      xleft = x[1:(ceiling(qt * m) - 1),,drop = F]; xright = x[(ceiling(qt * m)+1):m,, drop = F]
      
      leftnode$up[p] <- x[ceiling(qt * m),p, drop = F]
      rightnode$low[p] <- x[ceiling(qt * m) ,p, drop = F]
      
      node$leftchild = add_node(xleft,leftnode, d, thresh_marginal, thresh_interaction )
      node$rightchild = add_node(xright,rightnode, d, thresh_marginal, thresh_interaction)
    }
    return(node)
  }
  
  kdtree <- add_node(X, rootnode,  d, thresh_marginal = thresh_marginal, 
                     thresh_interaction = thresh_interaction)
  
  ahat <- ConfLevel(nd, alpha, method)
  kdtree <- SetBounds(kdtree, ahat, n)
  B <- matrix(nrow = 0, ncol = (2*d + 5))         # matrix that will hold selected regions
  B <- SelectNodes(kdtree, B, ahat, n)
  
  return(list(split_dir = split_dir, 
              hist = B,
              lower = lower,
              upper = upper))
}









