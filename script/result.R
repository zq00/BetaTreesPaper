# report simulation results 
library(tidyverse)
# 1. Changing parameters of the approximate algorithm 

# input file location 
file_loc <- "/work/pi_qianzhao_umass_edu/Beta_trees/approx_params/result/"
# parameters
L_vals <- c(15, 50, 100)
B_vals <- c(2000, 10^4, 5*10^4)

# read file 
result <- NULL
for(i in 1:3){
  for(j in 1:3){
    dat <- read.table(file = paste0(file_loc, "L", L_vals[i], "_B", B_vals[j], ".txt"))
    colnames(dat) <- c("Avg time", "tp", "fp", "y")
    cat("L = ", L_vals[i], ", B = ", B_vals[j], "\n")
    cat(colnames(dat), "\n")
    cat("Avg = ", round(colMeans(dat), 2), "\n")
    cat("sd = ", round(apply(dat, 2, sd) / sqrt(nrow(dat)), 2), "\n")
    
    new_val <- c(L_vals[i], B_vals[j], round(colMeans(dat), 2), round(apply(dat, 2, sd) / sqrt(nrow(dat)), 2) )
    result <- rbind(result, new_val)
  }
}
colnames(result) <- c("L", "B",  "avg_time", "TP", "FP", "Y", "sd(time)", "sd (TP)", "sd (FP)", "sd (Y)")
# write.table(result,paste0(file_loc, "approximate_sim.txt") , sep = ",")

# 2. changing number of true modes 

# input file location 
file_loc <- "/work/pi_qianzhao_umass_edu/Beta_trees/approx_sim/result/"
# parameters
m_vals <- c(5, 10, 15)
n_vals <- c(50000, 10^5, 2*10^5)
d_vals <- c(4, 8, 12)

result <- NULL
# read file 
for(i in 1:3){
  for(j in 1:3){
    for(k in 1:3){
      dat <- read.table(file = paste0(file_loc, "m", m_vals[i], "_n", n_vals[j],"_d", d_vals[k], ".txt"))
      colnames(dat) <- c("Avg time", "tp", "fp", "y")
      cat("m = ", m_vals[i], ", n = ", n_vals[j], ", d = ", d_vals[k],"\n")
      cat(colnames(dat), "\n")
      cat("Avg = ", round(colMeans(dat), 2), "\n")
      cat("sd = ", round(apply(dat, 2, sd) / sqrt(nrow(dat)), 2), "\n")
      new_val <- c(m_vals[i], n_vals[j], d_vals[k],round(colMeans(dat), 2), round(apply(dat, 2, sd) / sqrt(nrow(dat)), 2) )
      result <- rbind(result, new_val)
    }
  }
}
result <- as_tibble(result)
colnames(result) <- c("m", "n", "d", "avg_time", "TP", "FP", "Y", "sd(time)", "sd (TP)", "sd (FP)", "sd (Y)")
# write.table(result,paste0(file_loc, "approximate_sim.txt") , sep = ",")

# 3. computation cost of the exact algorithm 

file_loc <- "/work/pi_qianzhao_umass_edu/Beta_trees/exact/result/"
# parameters
n_vals <- c(2000, 5000)
cutoff <- c(4, 6)

# read file 
for(i in 1:2){
  for(j in 1:2){
      dat <- read.table(file = paste0(file_loc, "n", n_vals[i], "_c", cutoff[j], ".txt"))
      colnames(dat) <- c("Avg time", "tp", "fp", "y")
      cat("n = ", n_vals[i], ", cutoff = ", cutoff[j], "\n")
      cat(colnames(dat), "\n")
      cat("Avg = ", round(colMeans(dat), 2), "\n")
      cat("sd = ", round(apply(dat, 2, sd) / sqrt(nrow(dat)), 2), "\n")
  }
}

# 4. width of confidence intervals 

file_loc <- "/work/pi_qianzhao_umass_edu/Beta_trees/width/result/"
# parameters
n_val <- c(2000, 5000, 10^4, 10^5, 5*10^5)
alpha_val <- c(0.01, 0.05, 0.1, 0.2)
d_val <- c(1, 4, 8, 16)

# read file 
for(i in 1:4){
  d <- d_val[i]
  
  data <- tibble(.rows = 0)
  
  for(j in 1:4){
    for(k in 1:4){
      new_data <- read.table(file = paste0(file_loc, "d", d, "_n", n_val[j],
                                           "_a", alpha_val[k] ,".txt"))
      
      # compute summary statistics and write to table 
      new_tab <- tibble(n = n_val[j], alpha = alpha_val[k],
                        ci_width_25 = mean(new_data[,2]),
                        ci_width_25_sd = sd(new_data[,2]),
                        ci_width_50 = mean(new_data[,3]),
                        ci_width_50_sd = sd(new_data[,3]),
                        ci_width_75 = mean(new_data[,4]),
                        ci_width_75_sd = sd(new_data[,4]),
                        nbox = mean(new_data[,1])
      )
      data <- rbind(data, new_tab)
    }
  }
  write.table(data, file = paste0(file_loc,"standard_", d, ".txt"))
}


































