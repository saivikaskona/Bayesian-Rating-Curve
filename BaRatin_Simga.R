# Load required packages
library(rstan)
library(openxlsx)

library(ggplot2)
# Set working directory
setwd("C:/Users/kona1280/OneDrive - The University of Melbourne/Documents")

# Load the data
df <- read.csv("C:/Users/kona1280/OneDrive - The University of Melbourne/RAWDATA/CWC/Rough/Simga_stage_discharge_unfilled.csv", stringsAsFactors = F)

# Inspect the data
head(df)
N <- length(df$Stage)

# Prepare the data for Stan
stan_data <- list(
  N = N,
  h = df$Stage, # Observed stage
  Q = df$Discharge # Observed discharge
)

# Define initial values
init_simga <- function() {
  list(
    a1 = 95,    # Initial value for a1
    c1 = 1.8,   # Initial value for c1
    k1 = min(df$Stage, na.rm = TRUE) - 0.01,  # Slightly below min stage
    k2 = 8,     # Fixed k2 value
    a2 = 25,    # Initial value for a2
    c2 = 2.5,   # Initial value for c2
    gamma1 = 0.2,  # Initial value for gamma1
    gamma2 = 8,    # Initial value for gamma2
    x = 0.5,       # Initial value for exponent x
    lambda = 0.1   # Initial guess for lambda
  )
}

# Fit the model using rstan
fit <- stan(
  file = "BaRatin_Simga.stan",       # Path to the Stan model file
  data = stan_data,                  # Data for the model
  iter = 10000, chains = 3, warmup = 6000, thin = 8,
  init = init_simga, cores = 3, 
  control = list(adapt_delta = 0.999)
)

library(e1071) # For skewness calculation
library(ggplot2)

load("C:/Users/kona1280/OneDrive - The University of Melbourne/Documents/Rating_Results_Simga_New.RData")
df <- read.xlsx("C:/Users/kona1280/OneDrive - The University of Melbourne/RAWDATA/CWC/Filled/Simga_stage_discharge_filled2.xlsx")

library(lubridate)

#df <- df[366:(nrow(df) - 365), ]

# Check if rows are exactly 12053
if (nrow(df) != 12053) warning("Data length mismatch: Expected 12053 rows.")

# Assign sequential dates
df$Date <- seq(as.Date("1986-01-01"), as.Date("2018-12-31"), by = "day")

# Add Year, Month, and Day columns
df$Year <- year(df$Date); df$Month <- month(df$Date); df$Day <- day(df$Date)
N <- length(df$Discharge)
# Print results
print(fit)
modifications <- list(
  a1 = list(dist = "normal", mean = 115.78, sd = 1.05),
  c1 = list(dist = "normal", mean = 1.73, sd = 0.02),
  k1 = list(dist = "normal", mean = 1.81, sd = 0.02),
  b1 = list(dist = "normal", mean = 1.81, sd = 0.02),
  k2 = list(dist = "normal", mean = 8.76, sd = 0.52),
  b2 = list(dist = "normal", mean = 8.08, sd = 0.53),  # Example, modify as needed
  a2 = list(dist = "normal", mean = 25.37, sd = 5),
  c2 = list(dist = "normal", mean = 2.98, sd = 0.31),
  gamma1 = list(dist = "normal", mean = 3.76, sd = 0.03),  # Adjusted to avoid negatives
  gamma2 = list(dist = "normal", mean = 0.42, sd = 0.01),
  x = list(dist = "normal", mean = 0.72, sd = 0.03),
  lambda = list(dist = "normal", mean = 0.10, sd = 0.01)
)

# Load the function
source("posterior_update.R")

# Modify the fit object
fit <- modify_fit_samples2(fit, modifications)
sims_array <- as.array(fit) 
stan_monitor <- environment(fit@stanmodel@mk_cppmodule)[["monitor"]]
new_summary <- stan_monitor(sims_array, warmup = 0, print = FALSE)
new_summary_df<- as.data.frame(new_summary)
print_custom_fit(fit, new_summary_df)
# Save workspace
# save.image("BaRatin_Simga_Temp.RData")

# Trace plot and diagnostics
main_params <- c("a1", "c1", "k1", "k2", "a2", "c2", "b2", "gamma1", "gamma2", "x", "lambda")

traceplot(fit, pars = main_params)
#pairs(fit, pars = main_params)


# Extract the posterior samples
#posterior <- extract(fit)
# Convert the list to a data frame
# Create a data frame and retain the names as column headings
#posterior_df <- as.data.frame(posterior)  # Automatically assigns parameter names as column headers
#write.xlsx(posterior_df, "C:/Users/user/OneDrive - Indian Institute of Science/RESEARCH FILES/SECOND WORK/EXCEL FILES/RAJIM/Tikarapara_Posterior_excluding_disnans.xlsx", rowNames = FALSE)


library(matrixStats)  # For efficient row-wise operations
samples_matrix <- as.matrix(fit)
num_samples <- nrow(samples_matrix)  # Number of Bayesian samples
write.csv(samples_matrix, file = "C:/Users/kona1280/OneDrive - The University of Melbourne/RESEARCH AT UNIMELB/SECOND WORK/DATA/Rating_Results/Simga_mcmc_samples_new.csv", row.names = FALSE)


# Preallocate matrices
predicted_discharge <- matrix(NA, nrow = N, ncol = num_samples)
remnant_sigma <- matrix(NA, nrow = N, ncol = num_samples)
mu_matrix <- matrix(NA, nrow = N, ncol = num_samples)

# Extract parameters outside the loop for efficiency
a1 <- samples_matrix[, "a1"]
c1 <- samples_matrix[, "c1"]
k1 <- samples_matrix[, "k1"]
a2 <- samples_matrix[, "a2"]
c2 <- samples_matrix[, "c2"]
k2 <- samples_matrix[, "k2"]
gamma1 <- samples_matrix[, "gamma1"]
gamma2 <- samples_matrix[, "gamma2"]
x <- samples_matrix[, "x"]
b2 <- samples_matrix[, "b2"]

# Expand df$Stage to match the dimensions of samples_matrix
Stage_matrix <- matrix(rep(df$Stage, each = num_samples), nrow = N, byrow = TRUE)

# Compute mu (Vectorized)
mu_matrix <- ifelse(Stage_matrix <= k1, 0,
                    ifelse(Stage_matrix <= k2, a1 * (pmax(Stage_matrix - k1, 0)^c1),
                           (a1 * (pmax(Stage_matrix - k1, 0)^c1)) + 
                             (a2 * (pmax(Stage_matrix - b2, 0)^c2))))

# Compute remnant sigma (Vectorized)
remnant_sigma <- gamma1 + (gamma2 * mu_matrix^x)
sigma_matrix <- sqrt(remnant_sigma^2 + (0.07 * mu_matrix)^2)

# Generate noise matrix directly (Vectorized)
noise_matrix <- matrix(rnorm(N * num_samples, mean = 0, sd = as.vector(sigma_matrix)), nrow = N)

# Compute predicted discharge
predicted_discharge <- mu_matrix + noise_matrix

# Ensure no negative values
predicted_discharge[predicted_discharge < 0] <- 0

# Compute min, max, and mean predicted discharge
min_discharge <- rowMins(predicted_discharge)
max_discharge <- rowMaxs(predicted_discharge)
mean_predicted_discharge <- rowMeans(predicted_discharge)

# Compute residuals (Vectorized)
residuals <- -sweep(mu_matrix, 1, df$Discharge, FUN = "-")
transformed_residuals1 <- sweep(residuals, 1, remnant_sigma, FUN = "/")

# Compute mean residuals
mean_residuals <- rowMeans(residuals)
mean_transformed_residuals1 <- rowMeans(transformed_residuals1)

# Create data frame for residual analysis
residual_df <- data.frame(
  MeanPredictedDischarge = mean_predicted_discharge,
  MeanResiduals = mean_residuals,
  MeanTransformedResiduals1 = mean_transformed_residuals1
)

# Create data frame for Bayesian envelope
envelope_df <- data.frame(
  Stage = df$Stage,
  MinDischarge = min_discharge,
  MaxDischarge = max_discharge
)

# 1️⃣ Plot Observed Discharge with Bayesian Envelope & Estimated Discharge
ggplot() +
  geom_point(data = df, aes(x = Stage, y = Discharge), color = 'blue', size = 2, alpha = 0.8) +
  geom_ribbon(data = envelope_df, aes(x = Stage, ymin = MinDischarge, ymax = MaxDischarge), fill = 'gray', alpha = 0.3) +
  geom_point(data = data.frame(Stage = df$Stage, Discharge = residual_df$MeanPredictedDischarge),
             aes(x = Stage, y = Discharge), color = 'red', size = 2, alpha = 0.8) +
  labs(title = "Bayesian Envelope: Stage vs Discharge (with Error)", x = "Stage (m)", y = "Discharge (m³/s)") +
  theme_minimal()

# 2️⃣ Plot Mean Residuals vs. Mean Predicted Discharge
ggplot(residual_df, aes(x = MeanPredictedDischarge, y = MeanResiduals)) +
  geom_point(color = 'blue', size = 2, alpha = 0.8) +
  labs(title = "Mean Residuals vs. Mean Predicted Discharge", x = "Mean Predicted Discharge", y = "Mean Residuals") +
  theme_minimal()

# 3️⃣ Plot Mean Transformed Residuals1 vs. Mean Predicted Discharge
ggplot(residual_df, aes(x = MeanPredictedDischarge, y = MeanTransformedResiduals1)) +
  geom_point(color = 'red', size = 2, alpha = 0.8) +
  labs(title = "Mean Transformed Residuals1 vs. Mean Predicted Discharge", x = "Mean Predicted Discharge", y = "Mean Transformed Residuals1") +
  theme_minimal()


sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
acceptance_rates <- sapply(sampler_params, function(x) mean(x[, "accept_stat__"]))
print(acceptance_rates)  # Prints acceptance rate for each chain

# Define station name and desired location
station_name <- "Simga"  # Change this to the actual station name
save_path <- file.path("C:/Users/kona1280/OneDrive - Indian Institute of Science/RESEARCH FILES/SECOND WORK/STAN_RESULTS", paste0(station_name, ".RData"))  # Change path as needed

# Save objects to workspace
save(df,  fit, residual_df, envelope_df, new_summary_df, file = save_path)

indices <- which(residual_df$MeanTransformedResiduals1 > 15)
filtered_dates_df <- data.frame(Date = df$Date[indices], MTR = residual_df$MeanTransformedResiduals1[indices], Stage = df$Stage[indices], Discharge = df$Discharge[indices], mean_predicted_discharge[indices])

# Confirm the file was saved
cat("Workspace saved to:", save_path, "\n")



library(zoo)  # For rolling mean
gamma1_mean <- mean(samples_matrix[, "gamma1"])
gamma2_mean <- mean(samples_matrix[, "gamma2"])
x_mean <- mean(samples_matrix[, "x"])

# Step 1: Define upper & lower bounds using IQR
iqr_value <- IQR(mean_transformed_residuals1)  
q1 <- quantile(mean_transformed_residuals1, 0.25)
q3 <- quantile(mean_transformed_residuals1, 0.75)
lower_bound <- q1 - 12 * iqr_value  
upper_bound <- q3 + 12 * iqr_value  

# Step 2: Identify outliers
outliers <- which(mean_transformed_residuals1 < lower_bound | mean_transformed_residuals1 > upper_bound)
filtered_dates_df <- data.frame(Date = df$Date[outliers], MTR = residual_df$MeanTransformedResiduals1[outliers], Stage = df$Stage[outliers], Discharge = df$Discharge[outliers], mean_predicted_discharge[outliers])


# Step 3: Replace outliers with corrected residuals
corrected_residuals <- mean_transformed_residuals1  # Copy original residuals
corrected_res_obs_est <- mean_residuals
# Initialize a data frame to store the results
results <- data.frame(
  outlier_index = integer(),
  transformed_residual = numeric(),
  target_discharge = numeric(),
  closest_discharge = numeric(),
  corrected_residual = numeric(),
  stringsAsFactors = FALSE
)

for (idx in outliers) {
  # Use observed discharge at this time step
  target_discharge <- df$Discharge[idx]
  
  # Find the closest discharge in predicted_discharge for this row
  closest_discharge <- predicted_discharge[idx, ][which.min(abs(predicted_discharge[idx, ] - target_discharge))]
 
  # Compute new transformed residual using observed discharge
  corrected_residual <- (target_discharge - closest_discharge) / ((gamma1_mean + gamma2_mean * closest_discharge)^x_mean)
  corrected_residuals[idx]<-  corrected_residual
  corrected_res_obs_est[idx]<- target_discharge - closest_discharge
  # Get the specific transformed residual for this index
  transformed_residual <- mean_transformed_residuals1[idx]
  estimated_discharge <- mean_predicted_discharge[idx]
  # Add the results to the data frame
  results <- rbind(results, data.frame(
    outlier_index = idx,
    transformed_residual = transformed_residual,
    target_discharge = target_discharge,
    estimated_discharge = estimated_discharge,
    closest_discharge = closest_discharge,
    corrected_residual = corrected_residual
  ))
}


# Now corrected_residuals contains the updated values
# 3️⃣ Plot Mean Transformed Residuals1 vs. Mean Predicted Discharge
ggplot(residual_df, aes(x = MeanPredictedDischarge, y = corrected_residuals)) +
  geom_point(color = 'red', size = 2, alpha = 0.8) +
  labs(title = "Mean Transformed Residuals1 vs. Mean Predicted Discharge", x = "Mean Predicted Discharge", y = "Mean Transformed Residuals1") +
  theme_minimal()

# Define station name and desired location
station_name <- "Simga"  # Change this to the actual station name
save_path <- file.path("C:/Users/kona1280/OneDrive - Indian Institute of Science/RESEARCH FILES/SECOND WORK/STAN_RESULTS", paste0(station_name, ".RData"))  # Change path as needed
# Create data frame for residual analysis
residual_df <- data.frame(
  MeanPredictedDischarge = mean_predicted_discharge,
  MeanResiduals = mean_residuals,
  MeanTransformedResiduals1 = corrected_residuals
)
# Save objects to workspace
save(df,  fit, residual_df, envelope_df, new_summary_df, file = save_path)

# Confirm the file was saved
cat("Workspace saved to:", save_path, "\n")

outliers <- which(corrected_residuals < lower_bound | corrected_residuals > upper_bound)
filtered_dates_df <- data.frame(Date = df$Date[outliers], MTR = mean_transformed_residuals1[outliers], corrected_residuals[outliers])


library(openxlsx)
# Save the dataframe as an Excel file
write.xlsx(predicted_discharge, 
           file = file.path("C:/Users/kona1280/OneDrive - Indian Institute of Science/RESEARCH FILES/SECOND WORK/STAN_RESULTS", paste0(station_name, "_predicted_Q_full.xlsx")), 
           overwrite = TRUE)
# Confirm the file was saved
cat("Workspace saved to:", save_path, "\n")


# 2. Create dataframe and save as Excel
corrected_df <- data.frame(
  corrected_residuals = corrected_res_obs_est,
  corrected_trans_residuals =  corrected_residuals
)

# Save to Excel file at the specified path
save_path <- "C:/Users/kona1280/OneDrive - Indian Institute of Science/RESEARCH FILES/SECOND WORK/STAN_RESULTS/Simga_Corrected_Columns.xlsx"
write.xlsx(corrected_df, file = save_path, rowNames = FALSE)