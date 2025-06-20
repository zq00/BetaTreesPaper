This repository contains code to reproduce results from the paper "Beta-Trees: Multivariate Histograms With Confidence Statements"

Please install the R package BetaTree at the link: https://github.com/zq00/BetaTree

- **Visualization.qmd** reproduces results in Section 4 Summarizing data using the Beta-tree histogram.
  * Code in section Example 1 Beta-tree histogram with confidence intervals for univariate distributions reproduces Figure 1 in Appendix C.1
  * Code in section Example 3: Bivariate Gaussian distribution reproduces Figure 2 in Appendix C.2
  * Code in section Example 5: Ten-dimensional mixture of Gaussians reproduces Figure 3 in Appendix C.4 

- **ModeHunting.qmd** reproduces results in Section 5 Multivariate mode hunting

- **EMT.qmd** reproduces results in Section 6.1 Differentiating normal cells from EMT transitioning cells and in Appendix Section E. 

- **GvHD.qmd** reproduces results in Section 6.2 Identifying cell populations indicative of GvHD

- **src** folder contains the source code for the simulations
  * adaptive.R constructs adaptive Beta tree histogram
  * approximate_mode.R implements the approximate mode finding algorithm
  * sample_obs.R sample observations from various distributions used in the simulations
  * summar.R generates data for visualization e.g., calculating marginal and conditional distribution from the multivariate histogram and visualize the distribution
  * univariate.R generates data and plot of univariate examples
 
- **script** folder contains simulations that requires multiple repetitions that needs to be run on a cluster
  * exact_time.R reproduces results in Appendix D.3
  * parameter_approximate.R reproduces results in Appendix D.2.3 
  * sim_approximate reproduces results in Appendix D.2.2
  * result.R summarizes results into data tables in Appendix Table 1-5. 
