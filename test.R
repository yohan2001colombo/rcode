install.packages("devtools")
devtools::install_github("yohan2001colombo/nrprtwanov")

detach("package:npranov", unload=TRUE, force=TRUE)
remove.packages("npranov")

library(npranov)

set.seed(123)
soil_type <- factor(rep(c("Sandy", "Clay", "Loamy"), each = 30))  # Soil type (3 levels)
watering_frequency <- factor(rep(c("Low", "Medium", "High"), times = 30))  # Watering frequency (3 levels)
plant_height <- rnorm(n_obs, mean = 20, sd = 5) +
  as.numeric(soil_type) * 2 +
  as.numeric(watering_frequency) * 3 +
  as.numeric(soil_type) * as.numeric(watering_frequency) * 1.5 +
  rnorm(n_obs, sd = 5)  # Increased noise
# Create a data frame
data <- data.frame(soil_type, watering_frequency, plant_height)

krss(plant_height,soil_type,watering_frequency)

interaction_plot(data=data,soil_type,watering_frequency,plant_height)
main_effect_boxplot(data=data, soil_type,plant_height, "Main Effect of Soil Type")
main_effect_boxplot(data, watering_frequency,plant_height, "Main Effect of Tribe")

# Install and load your package
install.packages("devtools")
devtools::install_github("yohan2001colombo/nrprtwanov")
library(npranov)

########################################################################

# Set seed for reproducibility
set.seed(123)

# Define simulation parameters
n_simulations <- 1000  # Number of simulations
n_obs <- 90            # Number of observations per simulation
soil_levels <- c("Sandy", "Clay", "Loamy")  # Soil type levels
watering_levels <- c("Low", "Medium", "High")  # Watering frequency levels

# Initialize storage for results
results <- data.frame(
  simulation = 1:n_simulations,
  p_value_soil = numeric(n_simulations),
  p_value_watering = numeric(n_simulations),
  p_value_interaction = numeric(n_simulations)
)

# Run simulations
for (i in 1:n_simulations) {
  # Simulate data with smaller effects and more noise
  soil_type <- factor(rep(soil_levels, each = n_obs / length(soil_levels)))
  watering_frequency <- factor(rep(watering_levels, times = n_obs / length(watering_levels)))
  plant_height <- rnorm(n_obs, mean = 20, sd = 5) +
    as.numeric(soil_type) * 0.5 +  # Smaller main effect of soil type
    as.numeric(watering_frequency) * 0.75 +  # Smaller main effect of watering frequency
    as.numeric(soil_type) * as.numeric(watering_frequency) * 0.25 +  # Smaller interaction effect
    rnorm(n_obs, sd = 5)  # Increased noise

  # Create a data frame
  data <- data.frame(soil_type, watering_frequency, plant_height)

  # Apply your package's function
  krss_result <- krss(plant_height, soil_type, watering_frequency)

  # Store results
  results$p_value_soil[i] <- krss_result$P_value[1]
  results$p_value_watering[i] <- krss_result$P_value[2]
  results$p_value_interaction[i] <- krss_result$P_value[3]
}

# Analyze simulation results
# Calculate statistical power (proportion of significant results at alpha = 0.05)
power_soil <- mean(results$p_value_soil < 0.05)
power_watering <- mean(results$p_value_watering < 0.05)
power_interaction <- mean(results$p_value_interaction < 0.05)

# Print power results
cat("Power for Soil Type:", power_soil, "\n")
cat("Power for Watering Frequency:", power_watering, "\n")
cat("Power for Interaction:", power_interaction, "\n")

# Visualize p-value distributions
library(ggplot2)

# Plot p-values for Soil Type
ggplot(results, aes(x = p_value_soil)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black") +
  labs(title = "P-value Distribution for Soil Type", x = "P-value", y = "Frequency") +
  theme_minimal()

# Plot p-values for Watering Frequency
ggplot(results, aes(x = p_value_watering)) +
  geom_histogram(binwidth = 0.05, fill = "lightgreen", color = "black") +
  labs(title = "P-value Distribution for Watering Frequency", x = "P-value", y = "Frequency") +
  theme_minimal()

# Plot p-values for Interaction
ggplot(results, aes(x = p_value_interaction)) +
  geom_histogram(binwidth = 0.05, fill = "salmon", color = "black") +
  labs(title = "P-value Distribution for Interaction", x = "P-value", y = "Frequency") +
  theme_minimal()
