# Evaluate how choices of parameters affect speed and accuracy of the approximate algorithm
args <- commandArgs(trailingOnly = TRUE)
L <- as.numeric(args[1])
B = as.numeric(args[2])

cat("L = ", L,"\n B = ",  B , "\n")

library(tidyverse)
library(mvtnorm)
library(spatstat)
library(igraph)
library(tidyr)
library(goftest)

# directory for storage 
store_dir <- "/work/pi_qianzhao_umass_edu/Beta_trees/approx_params/result/"

# load functions
src_dir <- "/home/qianzhao_umass_edu/Research/Beta_trees/Paper/V2/src"
src_file_names <- list.files(src_dir, full.names = T)
for(f_path in src_file_names){source(f_path)}

# load functions in the package directory 
package_dir <- "/home/qianzhao_umass_edu/Research/Beta_trees/Package/R_04252025/"
file_names <- list.files(package_dir, full.names = T)
for(f_path in file_names){source(f_path)}

n <- 10^5
d <- 4
m <- 10

nrep <- 100

result <- summary_mode_finding(n = n, d = d, m = m, ndim = 4,
                               method = "approximate", 
                               nrep = nrep, 
                               L = L, B = B)

tp <- numeric(nrep)
fp <- numeric(nrep)
y <- numeric(nrep)

for(l in 1:nrep){
  x <- colSums(result$alg_distance[[l]] < 2)
  tp[l] <- sum(x > 0)
  fp[l] <- sum(x == 0)
  y[l] <- mean(rowSums(result$alg_distance[[l]] < 2))
}

# store results 
data <- cbind(result$alg_time, tp, fp, y)
write.table(data, file = paste0(store_dir, "L", L, "_B", B, ".txt"), sep = "\t")

cat("L = ", L, ",", "B = ", B, "\n")

cat("Average time = ", mean(result$alg_time), "; sd (time) = ", sd(result$alg_time) / sqrt(nrep) , "\n")
cat("Average TP = ", mean(tp), "; sd (TP) = ", sd(tp) / sqrt(nrep), "\n")
cat("Average FP = ", mean(fp), "; sd (FP) = ", sd(fp) / sqrt(nrep), "\n")
cat("Average # modes within a true mode = ", mean(y), "; sd (FP) = ", sd(y) / sqrt(nrep), "\n")


