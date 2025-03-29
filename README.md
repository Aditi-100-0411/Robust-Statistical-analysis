# RobustStatsBattle 
This project investigates the robustness of the L₂-Norm Ratio Estimator (LNRE) compared to the Minimum Density Power Divergence Estimator (MDPDE) in the presence of contaminated data. The study employs a synthetic dataset with 30% adversarial outliers (sampled from 
N
(
15
,
1
)
N(15,1)) to evaluate the resilience of both estimators in recovering the true parameters of the inlier distribution (
N
(
0
,
1
)
N(0,1)).

Using kernel-smoothed density estimates, LNRE optimizes a divergence measure based on the ratio of empirical and model-implied densities, while MDPDE minimizes a density power divergence. Results demonstrate that LNRE provides more stable estimates under contamination, with mean and variance estimates closer to the true inlier parameters.

The findings highlight LNRE’s potential as a robust alternative to conventional methods in outlier-prone settings, offering improved reliability for statistical inference in noisy environments.
