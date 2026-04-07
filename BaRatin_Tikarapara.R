library(rstan)
library(openxlsx)
library(e1071) 
library(ggplot2)
setwd("C:/Users/user/OneDrive - The University of Melbourne/Documents")
# Load the data
df <- read.csv("C:/Users/user/OneDrive - The University of Melbourne/RAWDATA/CWC/Rough/Tikarapara_stage_discharge_unfilled.csv", stringsAsFactors = F)

# Inspect the data
head(df)
N <- length(df$Discharge)

# Prepare the data for Stan
stan_data <- list(
  N = N,
  h = df$Stage, # Observed stage
  Q = df$Discharge # Observed discharge
)


    init_1 <- function() {
      list(
        a1 = 9.82,  # Initial value for a1
        c1 = 2.66,  # Initial value for c1
        k1 = 1.13,   # Initial value for k1\
        k2 = 5.5,   # Initial value for k2
        a2 = 152, # Initial value for a2
        c2 = 1.78, 
        gamma1 = 9.7,
        gamma2 = 0.34,  # Initial value for gamma2
        x = 0.65 # Initial value for exponent x
      )
    }  
  


# Fit the Stan model
fit <- stan(file = "BaRatin_Tikarapara3.stan", data = stan_data,
            iter =5000, chains = 3, warmup = 3000, thin = 4,init = init_1, cores = 15, control = list(adapt_delta = 0.99))

save.image("C:/Users/user/OneDrive - The University of Melbourne/Documents/Rating_Results_Tikarapara_New4.RData")




# Start from here
load("C:/Users/kona1280/OneDrive - The University of Melbourne/Documents/Rating_Results_Tikarapara_New4.RData")
# Print the results
print(fit)
# Load the data
df <- read.xlsx("C:/Users/kona1280/OneDrive - The University of Melbourne/RAWDATA/CWC/Filled/Tikarapara_stage_discharge_filled2.xlsx")
library(lubridate)

#df <- df[366:(nrow(df) - 365), ]

# Check if rows are exactly 12053
if (nrow(df) != 12053) warning("Data length mismatch: Expected 12053 rows.")

# Assign sequential dates
df$Date <- seq(as.Date("1986-01-01"), as.Date("2018-12-31"), by = "day")

# Add Year, Month, and Day columns
df$Year <- year(df$Date); df$Month <- month(df$Date); df$Day <- day(df$Date)
# Inspect the data
head(df)
N <- length(df$Discharge)

# Prepare the data for Stan
stan_data <- list(
  N = N,
  h = df$Stage, # Observed stage
  Q = df$Discharge # Observed discharge
)

# Trace plot and diagnostics
main_params <- c("a1", "c1", "b1", "b2", "a2", "c2","gamma1", "gamma2", "x")

traceplot(fit, pars = main_params)
samples_matrix <- as.matrix(fit)
#write.csv(samples_matrix, file = "C:/Users/kona1280/OneDrive - The University of Melbourne/RESEARCH AT UNIMELB/SECOND WORK/DATA/Rating_Results/Rajim_mcmc_samples_old.csv", row.names = FALSE)
num_samples <- nrow(samples_matrix)
# Pre-allocate matrices for efficiency
predicted_discharge <- matrix(NA, nrow = N, ncol = num_samples)
remnant_sigma <- matrix(NA, nrow = N, ncol = num_samples)
mu <- matrix(NA, nrow = N, ncol = num_samples)

# Extract the posterior samples
#posterior <- extract(fit)
# Convert the list to a data frame
# Create a data frame and retain the names as column headings
#posterior_df <- as.data.frame(posterior)  # Automatically assigns parameter names as column headers
#write.xlsx(posterior_df, "C:/Users/kona1280/OneDrive - Indian Institute of Science/RESEARCH FILES/SECOND WORK/EXCEL FILES/RAJIM/Tikarapara_Posterior_excluding_disnans.xlsx", rowNames = FALSE)

# Loop over Bayesian samples
for (i in 1:num_samples) {
  # Extract parameters from samples_matrix
  a1 <- samples_matrix[i, "a1"]
  c1 <- samples_matrix[i, "c1"]
  k1 <- samples_matrix[i, "k1"]
  a2 <- samples_matrix[i, "a2"]
  c2 <- samples_matrix[i, "c2"]
  k2 <- samples_matrix[i, "k2"]
  gamma1 <- samples_matrix[i, "gamma1"]
  gamma2 <- samples_matrix[i, "gamma2"]
  x <- samples_matrix[i, "x"]
  b1 <- samples_matrix[i, "b1"]
  b2 <- samples_matrix[i, "b2"]
  
  # Loop over time steps
  for (j in 1:N) {
    mu[j, i] <- ifelse(df$Stage[j] <= k1, 0,
                       ifelse(df$Stage[j] <= k2, a1 * (df$Stage[j] - b1)^c1,
                              a2 * (df$Stage[j] - b2)^c2))
    remnant_sigma[j, i] <- gamma1 + (gamma2 * mu[j, i]^x) # Use remnant sigma from Stan model
    sigma <- sqrt(remnant_sigma[j, i]^2 + (0.07 * mu[j, i])^2)
    noise <- max(rnorm(1, mean = 0, sd = sigma), -mu[j, i])
    predicted_discharge[j, i] <- mu[j, i] + noise
  }
}

min_discharge <- apply(predicted_discharge, 1, min)
max_discharge <- apply(predicted_discharge, 1, max)
mean_predicted_discharge <- rowMeans(predicted_discharge)
# Calculate the median predicted discharge for each time step
mean_predicted_discharge <- apply(predicted_discharge, 1, median)
residuals <- -(sweep(mu, 1, df$Discharge, FUN = "-"))
# Assuming residuals and df are already defined and loaded


# Calculate skewness of residuals for each sample
residual_skewness <- apply(residuals, 2, skewness)

# Identify the sample with the minimum residual skewness
min_skewness_index <- which.min(residual_skewness)

# Use the sample with the minimum residual skewness to predict the discharge
selected_predicted_discharge <- mu[, min_skewness_index]


# Compute transformed residuals
transformed_residuals1 <- sweep(residuals, 1, remnant_sigma, FUN = "/")

# Compute mean residuals per time step (averaging across Bayesian samples)
mean_residuals <- rowMeans(residuals)
mean_transformed_residuals1 <- rowMeans(transformed_residuals1)

# Create a data frame for plotting
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
  geom_point(data = data.frame(Stage = df$Stage, Discharge = selected_predicted_discharge),
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


sims_array <- as.array(fit) 
stan_monitor <- environment(fit@stanmodel@mk_cppmodule)[["monitor"]]
new_summary <- stan_monitor(sims_array, warmup = 0, print = FALSE)
new_summary_df<- as.data.frame(new_summary)
station_name <- "Tikarapara"  # Change this to the actual station name
save_path <- file.path("C:/Users/kona1280/OneDrive - Indian Institute of Science/RESEARCH FILES/SECOND WORK/STAN_RESULTS", paste0(station_name, ".RData"))  # Change path as needed

# Save objects to workspace
save(df,  fit, residual_df, envelope_df, new_summary_df, file = save_path)
cat("Workspace saved to:", save_path, "\n")

library(openxlsx)
# Save the dataframe as an Excel file
write.xlsx(predicted_discharge, 
           file = file.path("C:/Users/kona1280/OneDrive - Indian Institute of Science/RESEARCH FILES/SECOND WORK/STAN_RESULTS", paste0(station_name, "_predicted_Q_full.xlsx")), 
           overwrite = TRUE)
# Confirm the file was saved
cat("Workspace saved to:", save_path, "\n")