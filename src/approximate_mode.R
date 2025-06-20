# Functions for the approximate mode hunting algorithm 

# 1. An approximate algorithm to identify modes 
# INPUTS
# hist - estimated histogram 
#        Note: the (2d+1) : (2d+3) columns store the estimated density level 
# d - data dimension
# L - path length 
# B - number of paths 
FindModesApproximate <- function(hist, d, L, B){
  cluster <- list() # cluster stores the regions where reaching those nodes is equivalent to reaching the mode
  
  density <- hist[,(2*d+1)] # empirical density of each region
  ci <- hist[, (2*d+2):(2*d+3)] # lower and upper confidence bounds of each region
  node_order <- order(density, decreasing = T) # order nodes by decreasing order of empirical density
  
  # compute adjacency matrix
  adj <- compute_adjacency_mat(hist, d)
  # calculate the neighbor of the mode with highest density
  # compute a graph from the adj
  g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected")
  rm(adj) # the adjacency matrix is a large matrix object that takes up lots of memory when the number of region is large 
  
  mode <- node_order[1]
  if(nrow(hist) == 1) {
    cat("The histogram contains a single region!")
    return(list(mode = c(1),
                hist = hist))
  }
  
  cluster[[1]] <- find_neighbors(mode[1], g, ci)
  candidates<- setdiff(node_order[-1], cluster[[1]])
  
  cat("Total # nodes = ", nrow(hist), "\n")
  cat("# regions in 1st cluster = ", length(cluster[[1]]), "\n")
  cat("# candidate modes = ", length(candidates), "\n")
  
  nmodes <- 1
  while(length(candidates) > 0){
    connected <- connected_to_modes(candidates[1], mode, g, cluster, ci, L, B)
    
    if(length(candidates) %% 100 == 0) cat("\n",connected$nn, "\n")
    if(connected$connected == "unconnected"){
      nmodes <- nmodes + 1
      mode <- c(mode, candidates[1])
      
      cluster[[nmodes]] <- find_neighbors(mode[nmodes], g, ci)
      candidates<- setdiff(candidates[-1], cluster[[nmodes]])
      
      cat("# regions in", nmodes, "-th cluster = ", length(cluster[[nmodes]]), "\n")
      cat("# candidate modes = ", length(candidates), "\n")
    }else{
      candidates <- candidates[-1]
    }
  }
  return(
    list(mode = mode,
         hist = hist,
         cluster = cluster 
    )
  )
}

# 2. find regions that are "neighbors" of a mode, i.e., the upper CI are higher than the lower CI of the mode
# INPUTS
# mode - index of the mode 
# g - an adjacency graph 
# ci - lower and upper confidence bounds 
# OUTPUT
# a vector of "neighbors" of a mode
find_neighbors <- function(mode, g, ci){
 
  x <- neighbors <-  as.vector(neighbors(g, mode))
  x <- x[which(ci[x, 2] > ci[mode, 1])]
  
  while(length(x) > 0){
    new_nbd <- unlist(sapply(x, function(t) as.vector(neighbors(g, t))))
    
    x <- new_nbd[which(ci[new_nbd, 2] > ci[mode, 1])]
    x <- setdiff(x, neighbors)
    neighbors <- c(neighbors, new_nbd) %>% unique() 
  }
  return(c(neighbors, mode))
}

# 3. Use Monte Carlo to verify whether a candidate region is connected to any current modes
# INPUTS
# candidate_mode - index of candidate mode 
# curr_mode - indices of current modes 
# g - adjacency graph 
# cluster - neighborhood of current modes
# ci - confidence bounds
# L - path length 
# B - number of paths 
# OUTPUTS
# "connected" if the candidate mode is connected to any of the current modes and "unconnected" otherwise 
connected_to_modes <- function(candidate_mode, curr_mode, g, cluster, ci, L, B){
  connected <- FALSE
  nmode <- length(curr_mode)
  nn <- numeric(nmode) # counts the number of paths that do not reach any current mode
  for(b in 1:B){
    # generate a random walk on the graph starting from the candidate mode 
    path <- random_walk(g, candidate_mode, steps = L)
    
    for(j in 1:nmode){
      intersections <- intersect(path, cluster[[j]]) # intersection with clusters corresponding to each mode
      
      if(length(intersections) == 0) {
        nn[j] <- nn[j] + 1; 
        break; 
      }else{
        ind <- min(which(path %in% cluster[[j]]))
        
        val <- sapply(path[1:ind], function(t) ci[t,2] < ci[curr_mode[j],1] & ci[t,2] < ci[candidate_mode,1]) # T if the CI of a rectangle along the path dips below the CI of both candidate and current mode
        if(all(!val)){ # all values are F, i.e., no region along the path has lower density compared to i and j
          connected <- TRUE
          break; 
        }
      }
    }
    if(connected == TRUE) break; 
  }
  
  
  if(connected == TRUE){
    return(list(connected = "connected",
                nn = nn))
  }else{
    return(list(connected = "unconnected",
                nn = nn))
  }
}











