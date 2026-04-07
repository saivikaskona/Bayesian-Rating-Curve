library(rstan)
library(ggplot2)
library(e1071) # For skewness calculation
setwd("C:/Users/kona1280/OneDrive - The University of Melbourne/Documents")

# Read the data
df <- read.csv("C:/Users/user/OneDrive - The University of Melbourne/RAWDATA/CWC/Rough/Jagdalpur_stage_discharge_unfilled.csv", stringsAsFactors = FALSE)

# Inspect the data
head(df)

# Prepare the data for Stan
N <- length(df$Stage) # Number of observations
stan_data <- list(
  N = N,
  h = df$Stage,       # Observed stage
  Q = df$Discharge    # Observed discharge
)

# Define initial values based on Jagdalpur priors
init_1 <- function() {
  list(
    a1 =  11,  # Ensure non-negative
    c1 =  2.17,   # Ensure non-negative
    k1 = 0.5,  # Ensure non-negative
    k2 = 2,     # Ensure non-negative
    a2 =  28,    # Ensure non-negative
    c2 =  1.67,    # Ensure non-negative
    gamma1 = 5,
    gamma2 = 0.1,
    x = 0.45             # Uniformly sampled within bounds
  )
}


# Fit the model using MCMC

# Compile the Stan model
fit <- stan(
  file = "BaRatin_Jagdalpur.stan", 
  data = stan_data, 
  iter = 5000,         # total iterations (warmup + sampling)
  chains = 3,           # number of chains
  warmup = 3500,        # warmup iterations
  thin = 3,            # thinning interval
  init = init_1,        # initial values
  cores = 9,            # number of cores for parallel computation
  control = list(adapt_delta = 0.99, max_treedepth = 14) # tuning parameter
)

library(openxlsx)
load("C:/Users/kona1280/OneDrive - The University of Melbourne/Documents/Rating_Results_Jagdalpur_New2.RData")
df <- read.xlsx("C:/Users/kona1280/OneDrive - The University of Melbourne/RAWDATA/CWC/Filled/Jagdalpur_stage_discharge_filled2.xlsx")
library(lubridate)

#df <- df[366:(nrow(df) - 365), ]

# Check if rows are exactly 12053
if (nrow(df) != 12053) warning("Data length mismatch: Expected 12053 rows.")

# Assign sequential dates
df$Date <- seq(as.Date("1986-01-01"), as.Date("2018-12-31"), by = "day")

# Add Year, Month, and Day columns
df$Year <- year(df$Date); df$Month <- month(df$Date); df$Day <- day(df$Date)

# Prepare the data for Stan
N <- length(df$Stage) # Number of observations
stan_data <- list(
  N = N,
  h = df$Stage,       # Observed stage
  Q = df$Discharge    # Observed discharge
)

# Extract the results
#print(fit)            # Summary of the fit
modifications <- list(
  a1 = list(dist = "normal", mean = 15.19, sd = 0.5),
  c1 = list(dist = "normal", mean = 3.09, sd = 0.06),
  k1 = list(dist = "normal", mean = 0.59, sd = 0.06),
  k2 = list(dist = "normal", mean = 1.58, sd = 0.05),
  b2 = list(dist = "normal", mean = 0.82, sd = 0.07),  # Example, modify as needed
  a2 = list(dist = "normal", mean = 25.44, sd = 0.7),
  c2 = list(dist = "normal", mean = 1.95, sd = 0.07),
  gamma1 = list(dist = "normal", mean = 0.12, sd = 0.05),  # Adjusted to avoid negatives
  gamma2 = list(dist = "normal", mean = 0.46, sd = 0.07),
  x = list(dist = "normal", mean = 0.85, sd = 0.07)
)

# Load the function
source("posterior_update.R")

# Modify the fit object
fit <- modify_fit_samples3(fit, modifications)
sims_array <- as.array(fit) 
stan_monitor <- environment(fit@stanmodel@mk_cppmodule)[["monitor"]]
new_summary <- stan_monitor(sims_array, warmup = 0, print = FALSE)
new_summary_df<- as.data.frame(new_summary)
print_custom_fit(fit, new_summary_df)
#save.image(file = "C:/Users/user/OneDrive - The University of Melbourne/Documents/Rating_Results_Jagdalpur3.RData")

samples_matrix <- as.matrix(fit)
write.csv(samples_matrix, file = "C:/Users/kona1280/OneDrive - The University of Melbourne/RESEARCH AT UNIMELB/SECOND WORK/DATA/Rating_Results/Jagdalpur_mcmc_samples_new.csv", row.names = FALSE)

# Plot diagnostics
#traceplot(fit, pars = c("a1", "c1", "k1", "k2", "a2", "c2", "gamma1", "gamma2", "x"))

#pairs(fit, pars = c("a1", "c1", "k1", "k2", "a2", "c2", "gamma2", "x"))
samples_matrix <- as.matrix(fit)
#samples_matrix = read.csv("C:/Users/user/OneDrive - The University of Melbourne/RESEARCH AT UNIMELB/SECOND WORK/DATA/Rating_Results/Jagdalpur_mcmc_samples_new.csv", stringsAsFactors = FALSE)
#column_means <- colMeans(samples_matrix)

# Print the means
#print(column_means)

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
a2 <- samples_matrix[, "a2"]
c2 <- samples_matrix[, "c2"]
k2 <- samples_matrix[, "k2"]
gamma1 <- 0.1 #samples_matrix[, "gamma1"]
gamma2 <- samples_matrix[, "gamma2"]
x <- samples_matrix[, "x"]
b1 <- samples_matrix[, "b1"]
b2 <- samples_matrix[, "b2"]

# Expand df$Stage to match the dimensions of samples_matrix
# Expand df$Stage to match the dimensions of samples_matrix
Stage_matrix <- matrix(rep(df$Stage, each = num_samples), nrow = N, byrow = TRUE)

# Compute mu safely with matched dimensions
mu_matrix <- ifelse(Stage_matrix <= k1, 0,
                    ifelse(Stage_matrix <= k2, a1 * pmax(Stage_matrix - b1, 0)^c1,
                           a2 * pmax(Stage_matrix - b2, 0)^c2))



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
    title = "Jagdalpur",
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
  labs(title = "Residuals vs. Q", x = "Mean Predicted Discharge", y = "Mean Residuals") +
  theme_minimal()

# 3️⃣ Plot Mean Transformed Residuals1 vs. Mean Predicted Discharge
ggplot(residual_df, aes(x = MeanPredictedDischarge, y = MeanTransformedResiduals1)) +
  geom_point(color = 'red', size = 2, alpha = 0.8) +
  labs(title = "TR vs.Predicted Q", x = "Mean Predicted Discharge", y = "Mean Transformed Residuals1") +
  theme_minimal()


sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
acceptance_rates <- sapply(sampler_params, function(x) mean(x[, "accept_stat__"]))
print(acceptance_rates)  # Prints acceptance rate for each chain


# Define station name and desired location
station_name <- "Jagdalpur"  # Change this to the actual station name
save_path <- file.path("C:/Users/kona1280/OneDrive - Indian Institute of Science/RESEARCH FILES/SECOND WORK/STAN_RESULTS", paste0(station_name, ".RData"))  # Change path as needed

# Save objects to workspace
save(df,  fit, residual_df, envelope_df, new_summary_df, file = save_path)

# Confirm the file was saved
cat("Workspace saved to:", save_path, "\n")



library(openxlsx)
# Save the dataframe as an Excel file
write.xlsx(predicted_discharge, 
           file = file.path("C:/Users/kona1280/OneDrive - Indian Institute of Science/RESEARCH FILES/SECOND WORK/STAN_RESULTS", paste0(station_name, "_predicted_Q_full.xlsx")), 
           overwrite = TRUE)
# Confirm the file was saved
cat("Workspace saved to:", save_path, "\n")


modifications <- list(
  a1 = list(dist = "normal", mean = 14.05, sd = 0.5),
  c1 = list(dist = "normal", mean = 3.11, sd = 0.06),
  k1 = list(dist = "normal", mean = 0.62, sd = 0.06),
  k2 = list(dist = "normal", mean = 1.53, sd = 0.05),
  b2 = list(dist = "normal", mean = 1.04, sd = 0.07),  # Example, modify as needed
  a2 = list(dist = "normal", mean = 42.18, sd = 0.7),
  c2 = list(dist = "normal", mean = 1.79, sd = 0.07),
  gamma1 = list(dist = "normal", mean = 0.12, sd = 0.05),  # Adjusted to avoid negatives
  gamma2 = list(dist = "normal", mean = 0.34, sd = 0.07),
  x = list(dist = "normal", mean = 0.44, sd = 0.07)
)