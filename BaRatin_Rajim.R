library(rstan)
library(ggplot2)
setwd("C:/Users/kona1280/OneDrive - The University of Melbourne/Documents")
# Load the data
df <- read.csv("C:/Users/kona1280/OneDrive - The University of Melbourne/RAWDATA/CWC/Rough/Rajim_stage_discharge_unfilled.csv", stringsAsFactors = F)

# Inspect the data
head(df)
N <- length(df$Discharge)
#weights <- df$Discharge / sum(df$Discharge)
#weights <- rank(df$Discharge) / sum(rank(df$Discharge)) * length(df$Discharge)


# Prepare the data for Stan
stan_data <- list(
  N = N,
  h = df$Stage, # Observed stage
  Q = df$Discharge, # Observed discharge
  weights = weights
)

init_1 <- function() {
  list(
    a1 = rnorm(1, 70, 5),  # Initial value for a1
    c1 = 2.7,  # Initial value for c1
    k1 = rnorm(1, 1.02, 0.05),   # Initial value for k1
    gamma1 = runif(1, 0.1, 10),
    gamma2 = runif(1, 0.1, 10),  # Initial value for gamma2
    x = rnorm(1, 0.5, 0.04) # Initial value for exponent x
  )
}  

# Fit the Stan model
fit <- stan(file = "BaRatin_Rajim.stan", data = stan_data,
            iter = 6000, chains = 1, warmup =5000, thin = 2,
            init = init_1, cores = 6,  control = list(adapt_delta = 0.99, ))
library(e1071) # 
load("C:/Users/kona1280/OneDrive - The University of Melbourne/Documents/Rating_Results_Rajim_New3.RData")
df <- read.csv("C:/Users/kona1280/OneDrive - The University of Melbourne/RAWDATA/CWC/Filled/Filled/Rajim_stage_discharge_filled.csv", stringsAsFactors = F)
library(lubridate)
df <- df[366:(nrow(df) - 365), ]

# Check if rows are exactly 12053
if (nrow(df) != 12053) warning("Data length mismatch: Expected 12053 rows.")

# Assign sequential dates
df$Date <- seq(as.Date("1986-01-01"), as.Date("2018-12-31"), by = "day")

# Add Year, Month, and Day columns
df$Year <- year(df$Date); df$Month <- month(df$Date); df$Day <- day(df$Date)
library(readxl)  # Load package to read Excel files

# Set discharge to 0 where stage < 1.13
df$Discharge[df$Stage < 1.13] <- 0
# Inspect the data
head(df)
N <- length(df$Discharge)
#weights <- df$Discharge / sum(df$Discharge)
#weights <- rank(df$Discharge) / sum(rank(df$Discharge)) * length(df$Discharge)


# Prepare the data for Stan
stan_data <- list(
  N = N,
  h = df$Stage, # Observed stage
  Q = df$Discharge, # Observed discharge
  weights = weights
)

# Print the results
print(fit)
modifications <- list(
  a1 = list(dist = "normal", mean = 101.02, sd = 0.3),
  c1 = list(dist = "normal", mean = 2.43, sd = 0.06),
  k1 = list(dist = "normal", mean = 1.13, sd = 0.06),
  gamma1 = list(dist = "normal", mean = 0.12, sd = 0.01),  # Adjusted to avoid negatives
  gamma2 = list(dist = "normal", mean =2.65, sd = 0.07),
  x = list(dist = "normal", mean = 0.54, sd = 0.07)
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




main_params <- c("a1", "c1", "k1","gamma1", "gamma2", "x")    

traceplot(fit, pars = main_params)
samples_matrix <- as.matrix(fit)
write.csv(samples_matrix, file = "C:/Users/kona1280/OneDrive - The University of Melbourne/RESEARCH AT UNIMELB/SECOND WORK/DATA/Rating_Results/Rajim_mcmc_samples_new.csv", row.names = FALSE)
num_samples <- nrow(samples_matrix)
# Pre-allocate matrices for efficiency
predicted_discharge <- matrix(NA, nrow = N, ncol = num_samples)
remnant_sigma <- matrix(NA, nrow = N, ncol = num_samples)
mu <- matrix(NA, nrow = N, ncol = num_samples)


library(matrixStats)  # For efficient row-wise operations

num_samples <- nrow(samples_matrix)  # Number of Bayesian samples
# Preallocate matrices
predicted_discharge <- matrix(NA, nrow = N, ncol = num_samples)
remnant_sigma <- matrix(NA, nrow = N, ncol = num_samples)
mu_matrix <- matrix(NA, nrow = N, ncol = num_samples)

# Extract parameters outside the loop for efficiency
a1 <- samples_matrix[, "a1"]
c1 <- samples_matrix[, "c1"]
k1 <- samples_matrix[, "k1"]
gamma1 <- samples_matrix[, "gamma1"] 
gamma2 <- samples_matrix[, "gamma2"]
x <- samples_matrix[, "x"]

# Expand df$Stage to match the dimensions of samples_matrix
Stage_matrix <- matrix(rep(df$Stage, each = num_samples), nrow = N, byrow = TRUE)

# Compute mu (Vectorized)
mu_matrix <- ifelse(Stage_matrix <= k1, 0, a1 * (pmax(Stage_matrix - k1, 0)^c1))

# Compute remnant sigma (Vectorized)
remnant_sigma <- gamma1 + (gamma2 * mu_matrix^x)
sigma_matrix <- sqrt(remnant_sigma^2 + (0.07 * mu_matrix)^2)

# Generate noise matrix (Vectorized)
noise_matrix <- matrix(rnorm(N * num_samples, mean = 0, sd = as.vector(sigma_matrix)), nrow = N)

# Apply max noise adjustment to avoid negative discharge
noise_matrix <- pmax(noise_matrix, -mu_matrix)

# Compute predicted discharge
predicted_discharge <- mu_matrix + noise_matrix

# Compute min, max, and mean predicted discharge
min_discharge <- rowMins(predicted_discharge)
max_discharge <- rowMaxs(predicted_discharge)
mean_predicted_discharge <- rowMeans(predicted_discharge)

# Compute residuals (Vectorized)
residuals <- -sweep(mu_matrix, 1, df$Discharge, FUN = "-")

# Calculate skewness of residuals for each sample
residual_skewness <- apply(residuals, 2, skewness)

# Identify the sample with the minimum residual skewness
min_skewness_index <- which.min(residual_skewness)

# Use the sample with the minimum residual skewness to predict discharge
selected_predicted_discharge <- mu_matrix[, min_skewness_index]

# Compute transformed residuals (Vectorized)
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
  # Observed discharge (blue points)
  geom_point(
    data = df, 
    aes(x = Stage, y = Discharge), 
    color = 'blue', 
    size = 2.5,          # Increased from 2
    alpha = 0.8
  ) +
  
  # Bayesian envelope (gray ribbon)
  geom_ribbon(
    data = envelope_df, 
    aes(x = Stage, ymin = MinDischarge, ymax = MaxDischarge), 
    fill = 'gray', 
    alpha = 0.3
  ) +
  
  # Estimated discharge (red points)
  geom_point(
    data = data.frame(
      Stage = df$Stage, 
      Discharge = residual_df$MeanPredictedDischarge
    ),
    aes(x = Stage, y = Discharge), 
    color = 'red', 
    size = 2,          # Increased from 2
    alpha = 0.8
  ) +
  
  # Labels and titles
  labs(
    title = "Rajim",
    x = "Stage (m)", 
    y = "Discharge (m³/s)"
  ) +
  
  # Theme modifications
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 14),  # X-axis title
    axis.title.y = element_text(size = 14),  # Y-axis title
    axis.text = element_text(size = 12),     # Axis numbers
    plot.title = element_text(
      size = 16, 
      face = "bold",         # Bold title
      margin = margin(b = 15)  # Add space below title
    )
  )


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


# Define station name and desired location
station_name <- "Rajim"  # Change this to the actual station name
save_path <- file.path("C:/Users/kona1280/OneDrive - Indian Institute of Science/RESEARCH FILES/SECOND WORK/STAN_RESULTS", paste0(station_name, ".RData"))  # Change path as needed

# Save objects to workspace
save(df,  fit, residual_df, envelope_df, new_summary_df, file = save_path)

# Confirm the file was saved
cat("Workspace saved to:", save_path, "\n")

indices <- which(residual_df$MeanTransformedResiduals1 > 20)
filtered_dates_df <- data.frame(Date = df$Date[indices], MTR = residual_df$MeanTransformedResiduals1[indices], Stage = df$Stage[indices], Discharge = df$Discharge[indices], mean_predicted_discharge[indices])


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
indices <- which(residual_df$MeanTransformedResiduals1 > 20)
filtered_dates_df <- data.frame(Date = df$Date[indices], MTR = residual_df$MeanTransformedResiduals1[indices], Stage = df$Stage[indices], Discharge = df$Discharge[indices], mean_predicted_discharge[indices])


# Step 3: Replace outliers with corrected residuals
 corrected_trans_residuals <- mean_transformed_residuals1  # Copy original residuals
 corrected_residuals <- mean_residuals
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
  corrected_residuals[idx] <- target_discharge - closest_discharge
  # Compute new transformed residual using observed discharge
  corrected_trans_residual <- (target_discharge - closest_discharge) / ((gamma1_mean + gamma2_mean * closest_discharge)^x_mean)
   corrected_trans_residuals[idx]<-  corrected_trans_residual
  # Get the specific transformed residual for this index
  transformed_residual <- mean_transformed_residuals1[idx]
  estimated_discharge <- mean_predicted_discharge[idx]
  corrected_trans_residual <-  corrected_trans_residuals[idx] 
  # Add the results to the data frame
  results <- rbind(results, data.frame(
    outlier_index = idx,
    transformed_residual = transformed_residual,
    target_discharge = target_discharge,
    estimated_discharge = estimated_discharge,
    closest_discharge = closest_discharge,
    corrected_trans_residual = corrected_trans_residual
  ))
}


# Now  corrected_trans_residuals contains the updated values
# 3️⃣ Plot Mean Transformed Residuals1 vs. Mean Predicted Discharge
ggplot(residual_df, aes(x = MeanPredictedDischarge, y =  corrected_trans_residuals)) +
  geom_point(color = 'red', size = 2, alpha = 0.8) +
  labs(title = "Mean Transformed Residuals1 vs. Mean Predicted Discharge", x = "Mean Predicted Discharge", y = "Mean Transformed Residuals1") +
  theme_minimal()

# Define station name and desired location
station_name <- "Rajim"  # Change this to the actual station name
save_path <- file.path("C:/Users/kona1280/OneDrive - Indian Institute of Science/RESEARCH FILES/SECOND WORK/STAN_RESULTS", paste0(station_name, ".RData"))  # Change path as needed
# Create data frame for residual analysis
residual_df <- data.frame(
  MeanPredictedDischarge = mean_predicted_discharge,
  MeanResiduals = mean_residuals,
  MeanTransformedResiduals1 =  corrected_trans_residuals
)
# Save objects to workspace
save(df,  fit, residual_df, envelope_df, new_summary_df, file = save_path)

library(openxlsx)
# Save the dataframe as an Excel file
write.xlsx(predicted_discharge, 
           file = file.path("C:/Users/kona1280/OneDrive - Indian Institute of Science/RESEARCH FILES/SECOND WORK/STAN_RESULTS", paste0(station_name, "_predicted_Q_full.xlsx")), 
           overwrite = TRUE)
# Confirm the file was saved
cat("Workspace saved to:", save_path, "\n")

# 2. Create dataframe and save as Excel
corrected_df <- data.frame(
  corrected_residuals = corrected_residuals,
  corrected_trans_residuals = corrected_trans_residuals
)

# Save to Excel file at the specified path
save_path <- "C:/Users/kona1280/OneDrive - Indian Institute of Science/RESEARCH FILES/SECOND WORK/STAN_RESULTS/Rajim_Corrected_Columns.xlsx"
write.xlsx(corrected_df, file = save_path, rowNames = FALSE)