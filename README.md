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



###############CODE#####################

========================================================================================================================================================
========================================================================================================================================================


# Load necessary library
       library(stats)
 
      # Set seed for reproducibility
      set.seed(1235)
      
      # Step 1: Generate Binomial Data (n=1, p=0.3) for 50 points
      binomial_data <- rbinom(50, 1,  0.70)#this means 70% are 1 so 30% percent contamination
      binomial_data
      
    
      # Step 2: Introduce Outliers (Replace 0s with N(5,1) Samples)
      contaminated_data <- ifelse(binomial_data == 0, rnorm(sum(binomial_data == 0), mean = 15, sd = 1),  rnorm(sum(binomial_data == 1), mean = 0, sd = 1))
      hist(contaminated_data,freq= F)
      mean(contaminated_data)
      # Step 3: Convert to Data Frame
      contaminated_df <- data.frame(Original = binomial_data, Contaminated = contaminated_data)
  
      # Step 4: Visualization - Histogram of Contaminated vs Original Data
      hist(contaminated_data, col = "blue", main = "Contaminated Data vs. Original Binomial Data",
           xlab = "Value", ylab = "Frequency", breaks = 10, probability = TRUE)
      hist(binomial_data, col = rgb(1,0,0,0.5), add = TRUE, probability = TRUE)
      
      legend("topright", legend = c("Contaminated Data", "Original Binomial Data"),
             fill = c("blue", rgb(1,0,0,0.5)))
      
      # Step 5: Display the first few rows
      head(contaminated_df)
      
      
      
      
      
      
      
      
      
           
           # Silverman's Rule of Thumb for Bandwidth Selection
            h <- 1.06 * sd(contaminated_data) * length(contaminated_data)^(-1/5)
            
              # Define Gaussian kernel function
              gaussian_kernel <- function(x, y, h) {
                  dnorm((x - y) / h) / h
               }
           
              # Define empirical density estimate (g_star) using Gaussian Kernel
              g_star <- function(x) {
                  sapply(x, function(x_i) mean(sapply(contaminated_data, function(y) gaussian_kernel(x_i, y, h))))
                }
           
              # Ensure g_star is vectorized for plotting
              g_star_vec <- Vectorize(g_star)
          plot(g_star)
            
              # Define f_star (smoothed model density, assuming normal distribution)
              f_star <- function(x, mu) {
                  sapply(x, function(x_i) {
                      integral <- integrate(function(y) gaussian_kernel(x_i, y, h) * dnorm(y, mean = mu, sd = 1),
                                              lower = -Inf, upper = Inf, stop.on.error = FALSE)
                      if (!is.finite(integral$value)) return(0)  # Prevent non-finite values
                      return(integral$value)
                    })
               }
           
              # Ensure f_star is vectorized for plotting
              f_star_vec <- Vectorize(function(x) f_star(x, mu = mean(contaminated_data)))
           plot(f_star_vec)
            # Generate x values for plotting
              x_vals <- seq(min(contaminated_data), max(contaminated_data), length.out = 100)
            
              # Compute function values
              g_star_vals <- g_star(x_vals)
           f_star_vals <- f_star(x_vals, mean(contaminated_data))
            
             # Plot the results
              plot(x_vals, g_star_vals, type = "l", col = "blue", lwd = 2, ylim = range(c(g_star_vals, f_star_vals)),
                          main = "Kernel Smoothed Density Estimations", xlab = "x", ylab = "Density")
           lines(x_vals, f_star_vals, col = "red", lwd = 2)
            legend("topright", legend = c("g* (Kernel Density Estimate)", "f* (Model Smoothed Density)"),
                           col = c("blue", "red"), lwd = 2)
           
             # Set given alpha and beta values
              alpha <- 0.6
           beta <- 0.1
            
              # Define the LNRE function to be maximized
              product1 <- function(x, mu) {
                  g_star(x)^beta * f_star(x, mu)^(alpha - beta)
                }
            
              product2 <- function(x, mu) {
                  f_star(x, mu)^alpha
                }
           
             # Ensure `integrate()` works properly and returns single values
              integral1 <- function(mu) {
                  result <- integrate(function(x) product1(x, mu), lower = -Inf, upper = Inf, stop.on.error = FALSE)
                 
                  return(result$value)
                }
           
              integral2 <- function(mu) {
                  result <- integrate(function(x) product2(x, mu), lower = -Inf, upper = Inf, stop.on.error = FALSE)
                 
                  return(result$value)
                }
           
              # Define LNRE function properly
              lnre_value <- function(mu) {
                 int1 <- integral1(mu)
                  int2 <- integral2(mu)
                 
                    
                    return(-(alpha / (beta * (alpha - beta))) * log(int1) + (1 / beta) * log(int2))
                }
           
              # Vectorized LNRE function for plotting
              lnre_value_vec <- Vectorize(lnre_value)
            
              # Define range for mu and compute LNRE values
              mu_vals <- seq(-3, 3, length.out = 50)
           lnre_results <- sapply(mu_vals, lnre_value)
            
              # Plot the LNRE function
              plot(mu_vals, lnre_results, type = "l", col = "purple", lwd = 2,
                          xlab = "mu", ylab = "LNRE Value", main = "Plot of LNRE Function")
           
              # Optimize LNRE function
              result <- optim(par = mean(contaminated_data), fn =lnre_value, method = "L-BFGS-B",
                                   lower = -3, upper = 3)
           
             # Display optimized parameters
             cat("Optimized Mean:", result$par, "\n")
             cat("Maximum LNRE Value:", -result$value, "\n")
             
             
             #######applying the same data for mdpde measure#########
              
                # Load necessary libraries
                library(pracma)  # For numerical integration
             
                # Compute the density power divergence function H_alpha_n (Eq 3.4)
               H_alpha_n <- function(theta, data, alpha) {
                   term1 <- integrate(function(z) (dnorm(z, mean = theta, sd = 1))^(1 + alpha), -Inf, Inf)$value
                    term2 <- mean((dnorm(data, mean = theta, sd = 1))^alpha)
                    return(term1 - (1 + 1 / alpha) * term2)
                  }
             
               
                # Optimization function to get MDPDE (Minimization of H_alpha_n)
                MDPDE <- function(data, alpha) {
                    result <- optim(par = mean(contaminated_data), fn = function(theta) H_alpha_n(theta, data, alpha), 
                                       method = "L-BFGS-B", lower = -Inf, upper = Inf)
                    return(result$par)  # Extract optimal parameter
                  }
              
               # Load contaminated data (Assuming `contaminated_data` is already available)
                data <- contaminated_data
             mean(contaminated_data)
                # Choose alpha (Robustness parameter)
               alpha <- 0.5  # Adjust based on robustness needs
             
                # Compute MDPDE estimate
                theta_hat <- MDPDE(data, alpha)
             
                # Print the estimated parameter
                cat("Sample Mean", mean(contaminated_data),"\n")
                cat("Estimated Parameter (MDPDE) =", theta_hat, "\n")
                cat("Optimized Mean(LNRE):", result$par, "\n")
                cat("Maximum LNRE Value:", -result$value, "\n")
                
                
                
                
                # Set seed for reproducibility
                set.seed(1235)
                
                # Step 1: Generate Binomial Data (n=1, p=0.3) for 50 points
                binomial_data <- rbinom(50, 1, 0.70) # 70% are 1 (30% contamination)
                binomial_data
                
                # Step 2: Introduce Outliers (Replace 0s with N(15,1) Samples)
                contaminated_data <- ifelse(binomial_data == 0, 
                                            rnorm(sum(binomial_data == 0), mean = 15, sd = 1), 
                                            rnorm(sum(binomial_data == 1), mean = 0, sd = 1))
                
                hist(contaminated_data, freq = F)
                mean(contaminated_data)
                sd(contaminated_data) # Original SD of contaminated data
                
                # Silverman's Rule of Thumb for Bandwidth Selection
                h <- 1.06 * sd(contaminated_data) * length(contaminated_data)^(-1/5)
                
                # Define Gaussian kernel function
                gaussian_kernel <- function(x, y, h) {
                  dnorm((x - y) / h) / h
                }
                
                # Define empirical density estimate (g_star) using Gaussian Kernel
                g_star <- function(x) {
                  sapply(x, function(x_i) mean(sapply(contaminated_data, function(y) gaussian_kernel(x_i, y, h))))
                }
                
                # Define f_star (smoothed model density, assuming normal distribution)
                f_star <- function(x, sigma) {
                  sapply(x, function(x_i) {
                    integral <- integrate(function(y) gaussian_kernel(x_i, y, h) * dnorm(y, mean = 0, sd = sigma),
                                          lower = -Inf, upper = Inf, stop.on.error = FALSE)
                    if (!is.finite(integral$value)) return(0)
                    return(integral$value)
                  })
                }
                
                # Vectorize functions
                g_star_vec <- Vectorize(g_star)
                f_star_vec <- Vectorize(function(x) f_star(x, sigma = sd(contaminated_data)))
                
                # Generate x values for plotting
                x_vals <- seq(min(contaminated_data), max(contaminated_data), length.out = 100)
                
                # Compute function values
                g_star_vals <- g_star(x_vals)
                f_star_vals <- f_star(x_vals, sd(contaminated_data))
                
                # Plot the results
                plot(x_vals, g_star_vals, type = "l", col = "blue", lwd = 2, ylim = range(c(g_star_vals, f_star_vals)),
                     main = "Kernel Smoothed Density Estimations", xlab = "x", ylab = "Density")
                lines(x_vals, f_star_vals, col = "red", lwd = 2)
                legend("topright", legend = c("g* (Kernel Density Estimate)", "f* (Model Smoothed Density)"),
                       col = c("blue", "red"), lwd = 2)
                
                # Set given alpha and beta values
                alpha <- 0.6
                beta <- 0.1
                
                # Define the LNRE function for sigma estimation
                product1_sigma <- function(x, sigma) {
                  g_star(x)^beta * f_star(x, sigma)^(alpha - beta)
                }
                
                product2_sigma <- function(x, sigma) {
                  f_star(x, sigma)^alpha
                }
                
                integral1_sigma <- function(sigma) {
                  result <- integrate(function(x) product1_sigma(x, sigma), lower = -Inf, upper = Inf, stop.on.error = FALSE)
                  return(result$value)
                }
                
                integral2_sigma <- function(sigma) {
                  result <- integrate(function(x) product2_sigma(x, sigma), lower = -Inf, upper = Inf, stop.on.error = FALSE)
                  return(result$value)
                }
                
                lnre_value_sigma <- function(sigma) {
                  int1 <- integral1_sigma(sigma)
                  int2 <- integral2_sigma(sigma)
                  return(-(alpha / (beta * (alpha - beta))) * log(int1) + (1 / beta) * log(int2))
                }
                
                # Vectorized LNRE function for sigma
                lnre_value_sigma_vec <- Vectorize(lnre_value_sigma)
                
                # Define range for sigma and compute LNRE values
                sigma_vals <- seq(0.5, 5, length.out = 50)
                lnre_results_sigma <- sapply(sigma_vals, lnre_value_sigma)
                
                # Plot the LNRE function for sigma
                plot(sigma_vals, lnre_results_sigma, type = "l", col = "purple", lwd = 2,
                     xlab = "sigma", ylab = "LNRE Value", main = "LNRE Function for Sigma")
                
                # Optimize LNRE function for sigma
                result_sigma <- optim(par = sd(contaminated_data), fn = lnre_value_sigma, method = "L-BFGS-B",
                                      lower = 0.1, upper = 5)
                
                # Display optimized parameters for sigma
                cat("Original Sample SD:", sd(contaminated_data), "\n")
                cat("Optimized Sigma (LNRE):", result_sigma$par, "\n")
                cat("Minimum LNRE Value:", result_sigma$value, "\n")
                
                ####### MDPDE for Sigma Estimation #########
                
                # Compute the density power divergence function for sigma
                H_alpha_n_sigma <- function(sigma, data, alpha) {
                  term1 <- integrate(function(z) (dnorm(z, mean = 0, sd = sigma))^(1 + alpha), -Inf, Inf)$value
                  term2 <- mean((dnorm(data, mean = 0, sd = sigma))^alpha)
                  return(term1 - (1 + 1 / alpha) * term2)
                }
                
                # Optimization function to get MDPDE for sigma
                MDPDE_sigma <- function(data, alpha) {
                  result <- optim(par = sd(data), fn = function(sigma) H_alpha_n_sigma(sigma, data, alpha), 
                                  method = "L-BFGS-B", lower = 0.1, upper = 5)
                  return(result$par)
                }
                
                # Compute MDPDE estimate for sigma
                sigma_hat <- MDPDE_sigma(contaminated_data, alpha = 0.5)
                
                # Print the estimated parameters
                cat("\n=== Final Results ===\n")
                cat("Original Sample SD:", sd(contaminated_data), "\n")
                cat("Estimated Sigma (MDPDE):", sigma_hat, "\n")
                cat("Optimized Sigma (LNRE):", result_sigma$par, "\n")
                # Print the estimated parameter
                cat("Sample Mean", mean(contaminated_data),"\n")
                cat("Estimated Parameter (MDPDE) =", theta_hat, "\n")
                cat("Optimized Mean(LNRE):", result$par, "\n")
                cat("Maximum LNRE Value:", -result$value, "\n")
                
