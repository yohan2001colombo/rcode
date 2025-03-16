Krss <- function(response, A, B, n_perm = 1000) {
  library(kernlab)  # Kernel regression
  library(ggplot2)  # Visualization
  library(MASS)     # For truehist
  
  # Gaussian kernel function
  gaussian_kernel <- function(x, y, sigma = 1) {
    exp(- (x - y)^2 / (2 * sigma^2))
  }
  
  # Create kernel matrices for two categorical factors
  kernel_matrix <- function(factor_levels, sigma = 1) {
    n <- length(factor_levels)
    K <- matrix(0, n, n)
    
    for (i in 1:n) {
      for (j in 1:n) {
        K[i, j] <- gaussian_kernel(as.numeric(factor_levels[i]), as.numeric(factor_levels[j]), sigma)
      }
    }
    return(K)
  }
  
  # Compute kernel matrices
  K1 <- kernel_matrix(A)  # Main effect of A
  K2 <- kernel_matrix(B)  # Main effect of B
  K12 <- K1 * K2          # Interaction effect
  
  # Function to fit kernel ridge regression and compute RSS
  compute_rss <- function(K, response) {
    model <- ksvm(K, response, type = "nu-svr", kernel = "matrix", C = 1)
    predictions <- predict(model)
    rss <- sum((response - predictions)^2)
    return(rss)
  }
  
  # Permutation test function
  permutation_test <- function(K_full, K_reduced, response, n_perm) {
    # Fit full and reduced models
    rss_full <- compute_rss(K_full, response)
    rss_reduced <- compute_rss(K_reduced, response)
    observed_statistic <- rss_reduced - rss_full  # Test statistic
    
    # Permutation loop
    perm_statistics <- numeric(n_perm)
    for (i in 1:n_perm) {
      perm_response <- sample(response)
      rss_full_perm <- compute_rss(K_full, perm_response)
      rss_reduced_perm <- compute_rss(K_reduced, perm_response)
      perm_statistics[i] <- rss_reduced_perm - rss_full_perm
    }
    
    # Compute p-value
    p_value <- mean(perm_statistics >= observed_statistic)
    return(p_value)
  }
  
  # Compute p-values
  p_value_main1 <- permutation_test(K1 + K2, K2, response, n_perm)  # Test for A
  p_value_main2 <- permutation_test(K1 + K2, K1, response, n_perm)  # Test for B
  p_value_interaction <- permutation_test(K1 + K2 + K12, K1 + K2, response, n_perm)  # Test for interaction
  
  # Create significance labels
  Significance <- ifelse(c(p_value_main1, p_value_main2, p_value_interaction) <= 0.05, "***", "")
  
  # Create a data frame similar to the ANOVA output
  anova_results <- data.frame(
    Factor = c("FactorA", "FactorB", "Interaction"),
    P_value = c(p_value_main1, p_value_main2, p_value_interaction),
    Significance = Significance
  )
  
  # Print the table
  print(anova_results, row.names = FALSE)
}
###################### EX:1 ###################################
# Define the data
Location <- as.factor(c(rep("Olympia", 6), rep("Ventura", 6),
                        rep("Northampton", 6), rep("Burlington", 6)))
Tribe <- as.factor(c(rep(c("Jedi", "Sith"), 12)))
Midichlorians <- c(10, 4, 12, 5, 15, 4, 15, 9, 15, 11, 18, 12,
                   8, 13, 8, 15, 10, 17, 22, 22, 20, 22, 20, 25)

##Compare with existing parametric two-way ANOVA
# Perform the two-way ANOVA
anova_result <- aov(Midichlorians ~ Location * Tribe, data = data)
summary(anova_result)

# Perform krss
Krss(Midichlorians,Location,Tribe)

###Test robustness under various violations
# Function to simulate data with violations
simulate_data <- function(effect_size = 0.5, violation_type = "none") {
  Location <- as.factor(c(rep("Olympia", 6), rep("Ventura", 6),
                          rep("Northampton", 6), rep("Burlington", 6)))
  Tribe <- as.factor(c(rep(c("Jedi", "Sith"), 12)))
  Midichlorians <- rnorm(24, mean = 10, sd = 2)  # Baseline
  
  # Add effect size
  Midichlorians[Tribe == "Sith"] <- Midichlorians[Tribe == "Sith"] + effect_size
  
  # Introduce violations
  if (violation_type == "non_normal") {
    Midichlorians <- rexp(24, rate = 1)  # Exponential distribution (non-normal)
  } else if (violation_type == "heteroscedastic") {
    Midichlorians[Tribe == "Sith"] <- Midichlorians[Tribe == "Sith"] + rnorm(12, mean = 0, sd = 3)  # Heteroscedasticity
  }
  
  return(list(Location = Location, Tribe = Tribe, Midichlorians = Midichlorians))
}

# Set a random seed for reproducibility
set.seed(123)

# Test robustness under non-normality
non_normal_data <- simulate_data(violation_type = "non_normal")
cat("Results for non-normal data:\n")
Krss(non_normal_data$Midichlorians, non_normal_data$Location, non_normal_data$Tribe)

# Test robustness under heteroscedasticity
hetero_data <- simulate_data(violation_type = "heteroscedastic")
cat("\nResults for heteroscedastic data:\n")
Krss(hetero_data$Midichlorians, hetero_data$Location, hetero_data$Tribe)

# Test with original data (no violations)
original_data <- simulate_data(violation_type = "none")
cat("\nResults for original data (no violations):\n")
Krss(original_data$Midichlorians, original_data$Location, original_data$Tribe)


#####Power Analyis 
# Function to simulate data with a given effect size
simulate_data <- function(n, effect_size) {
  Location <- as.factor(rep(c("Olympia", "Ventura", "Northampton", "Burlington"), each = n / 4))
  Tribe <- as.factor(rep(c("Jedi", "Sith"), each = n / 2))
  Midichlorians <- rnorm(n, mean = 10, sd = 2)  # Baseline
  Midichlorians[Tribe == "Sith"] <- Midichlorians[Tribe == "Sith"] + effect_size  # Add effect
  return(list(Location = Location, Tribe = Tribe, Midichlorians = Midichlorians))
}

# Power analysis function
power_analysis <- function(effect_size, n_range, n_sim = 100, n_perm = 100) {
  power_results <- data.frame(n = n_range, power = numeric(length(n_range)))
  
  for (i in 1:length(n_range)) {
    n <- n_range[i]
    significant <- numeric(n_sim)
    
    for (j in 1:n_sim) {
      data <- simulate_data(n, effect_size)
      krss_result <- Krss(data$Midichlorians, data$Location, data$Tribe, n_perm = n_perm)
      significant[j] <- krss_result$P_value[2] <= 0.05  # Check significance for FactorB (Tribe)
    }
    
    power_results$power[i] <- mean(significant)
  }
  
  return(power_results)
}

# Set parameters
effect_size <- 2  # Hypothesized effect size
n_range <- seq(20, 100, by = 20)  # Sample sizes to test
n_sim <- 50  # Number of simulations per sample size
n_perm <- 50  # Number of permutations per simulation

# Run power analysis
set.seed(123)  # For reproducibility
power_results <- power_analysis(effect_size, n_range, n_sim, n_perm)

# Print results
print(power_results)

# Plot power vs sample size
plot(power_results$n, power_results$power, type = "b", 
     xlab = "Sample Size", ylab = "Power", 
     main = "Power Analysis for Krss Function")

# Load the pwr package
library(pwr)

# Function to compute power for parametric two-way ANOVA
power_parametric_anova <- function(effect_size, n_range, k_A = 4, k_B = 2, alpha = 0.05) {
  power_results <- data.frame(n = n_range, power = numeric(length(n_range)))
  
  for (i in 1:length(n_range)) {
    n <- n_range[i]
    # Total sample size
    N <- n * k_A * k_B
    
    # Degrees of freedom for FactorA, FactorB, and Interaction
    df_A <- k_A - 1
    df_B <- k_B - 1
    df_AB <- (k_A - 1) * (k_B - 1)
    
    # Effect size (Cohen's f) for two-way ANOVA
    f <- effect_size
    
    # Compute power using pwr.anova.test (for one-way ANOVA, adjusted for two-way)
    # We use the degrees of freedom for the factor of interest (e.g., FactorA)
    power_results$power[i] <- pwr.anova.test(k = k_A, n = n, f = f, sig.level = alpha)$power
  }
  
  return(power_results)
}

# Set parameters
effect_size <- 0.5  # Hypothesized effect size (Cohen's f)
n_range <- seq(20, 100, by = 20)  # Sample sizes to test
k_A <- 4  # Number of levels for FactorA (Location)
k_B <- 2  # Number of levels for FactorB (Tribe)
alpha <- 0.05  # Significance level

# Run power analysis
power_results_parametric <- power_parametric_anova(effect_size, n_range, k_A, k_B, alpha)

# Print results
print(power_results_parametric)

# Plot power vs sample size
plot(power_results_parametric$n, power_results_parametric$power, type = "b", 
     xlab = "Sample Size", ylab = "Power", 
     main = "Power Analysis for Parametric Two-Way ANOVA")


