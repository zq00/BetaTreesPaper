This repository contains code to reproduce results from the paper "Beta-Trees: Multivariate Histograms With Confidence Statements"

Please install the R package BetaTree at the link: https://github.com/zq00/BetaTree. The R package includes code to compute the Beta-tree histogram, adaptive Beta-tree histogram, visualize two-dimensional marginal and conditional distributions, and both exact and approximate mode hunting algorithms. 

- **Visualization.qmd** reproduces results in Section 4 Summarizing data using the Beta-tree histogram.

- **ModeHunting.qmd** reproduces results in Section 5 Multivariate mode hunting

- **EMT.qmd** reproduces results in Section 6.1 Differentiating normal cells from EMT transitioning cells and in Appendix Section E. 

- **GvHD.qmd** reproduces results in Section 6.2 Identifying cell populations indicative of GvHD.

- **src** folder contains the source code for the simulations
  * compute_marginal.R used in Section 4 Example 4 transform the data to F-hat(Xj), F-hat is the estimated marginal distribution, construct the adaptive Beta-tree histogram for the transformed data, and convert the histogram to the scale of the original data. 
  * sample_obs.R sample observations from various distributions used in the simulations
  * summary.R computes the probability mass in each region in the examples shown in the paper; computes summary table of simulation results of the approximate mode hunting algorithm 
  * univariate.R generates data and plot of univariate examples
 
- **script** folder contains simulations that requires multiple repetitions that needs to be run on a cluster
  * exact_time.R reproduces results in Appendix D.3
  * parameter_approximate.R reproduces results in Appendix D.2.3 
  * sim_approximate reproduces results in Appendix D.2.2
  * result.R summarizes results into data tables in Appendix Table 1-5.
 

