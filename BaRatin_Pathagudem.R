library(rstan)
library(ggplot2)
setwd("C:/Users/kona1280/OneDrive - The University of Melbourne/Documents")

# Load the data
df2 <- read.csv("C:/Users/kona1280/OneDrive - The University of Melbourne/RAWDATA/CWC/Rough/Pathagudem_stage_discharge_unfilled.csv", stringsAsFactors = F)

# Inspect the data
head(df)
N <- length(df$Discharge)

# Compute weights (higher flows get more weight)
#weights <- df$Discharge / sum(df$Discharge)

# Prepare the data for Stan
stan_data <- list(
  N = N,
  h = df$Stage,
  Q = df$Discharge
  #weights = weights  # Pass weights to Stan
)

init_1 <- function() {
  list(
    a1 = rnorm(1, 236.12, 50),
    c1 = rnorm(1, 1.67, 0.0375),
    k1 = rnorm(1, 0.5, 0.1),
    gamma1 = runif(1, 0.1, 1000),
    gamma2 = runif(1, 0.1, 100),
    x = rnorm(1, 0.7, 0.2)
  )
}

# Fit the Stan model
fit <- stan(file = "BaRatin_Pathagudem.stan", data = stan_data,
            iter = 3000, chains = 3, warmup = 1500, thin = 5,
            init = init_1, cores = 15, control = list(adapt_delta = 0.99))

library(openxlsx)
load("C:/Users/kona1280/OneDrive - The University of Melbourne/Documents/Rating_Results_Pathagudem_New2.RData")
df <- read.xlsx("C:/Users/kona1280/OneDrive - The University of Melbourne/RAWDATA/CWC/Filled/Pathagudem_stage_discharge_filled2.xlsx")
library(lubridate)

#df <- df[366:(nrow(df) - 365), ]

# Check if rows are exactly 12053
if (nrow(df) != 12053) warning("Data length mismatch: Expected 12053 rows.")

# Assign sequential dates
df$Date <- seq(as.Date("1986-01-01"), as.Date("2018-12-31"), by = "day")

# Add Year, Month, and Day columns
df$Year <- year(df$Date); df$Month <- month(df$Date); df$Day <- day(df$Date)

N <- length(df$Discharge)

# Compute weights (higher flows get more weight)
#weights <- df$Discharge / sum(df$Discharge)

# Prepare the data for Stan
stan_data <- list(
  N = N,
  h = df$Stage,
  Q = df$Discharge
  #weights = weights  # Pass weights to Stan
)

# Print the results
print(fit) 
sims_array <- as.array(fit) 
stan_monitor <- environment(fit@stanmodel@mk_cppmodule)[["monitor"]]
new_summary <- stan_monitor(sims_array, warmup = 0, print = FALSE)
new_summary_df<- as.data.frame(new_summary)



main_params <- c("a1", "c1", "k1" , "gamma1", "gamma2", "x")
traceplot(fit, pars = main_params)

samples_matrix <- as.matrix(fit)

# Generate predicted discharge values
num_samples <- nrow(samples_matrix)
write.csv(samples_matrix, file = "C:/Users/kona1280/OneDrive - The University of Melbourne/RESEARCH AT UNIMELB/SECOND WORK/DATA/Rating_Results/Pathagudem_mcmc_samples_new.csv", row.names = FALSE)


predicted_discharge <- matrix(NA, nrow = N, ncol = num_samples)
remnant_sigma <- matrix(NA, nrow = N, ncol = num_samples)
mu <- matrix(NA, nrow = N, ncol = num_samples)



for (i in 1:num_samples) {
  a1 <- samples_matrix[i, "a1"]
  c1 <- samples_matrix[i, "c1"]
  k1 <- samples_matrix[i, "k1"]
  gamma1 <- samples_matrix[i, "gamma1"]
  gamma2 <- samples_matrix[i, "gamma2"]
  x <- samples_matrix[i, "x"]
  
  for (j in 1:N) {
    mu[j, i] <- ifelse(df$Stage[j] <= k1, 0, a1 * (df$Stage[j] - k1)^c1)
    remnant_sigma[j, i] <- gamma1 + (gamma2 * mu[j, i]^x) # Use remnant sigma from Stan model
    sigma <- sqrt(remnant_sigma[j, i]^2 + (0.07 * mu[j, i])^2)
    noise <- max(rnorm(1, mean = 0, sd = sigma), -mu[j, i])
    predicted_discharge[j, i] <- mu[j, i] + noise
  }
}

# Compute min and max predicted discharge including error component
min_discharge <- apply(predicted_discharge, 1, min)
max_discharge <- apply(predicted_discharge, 1, max)
mean_predicted_discharge <- rowMeans(predicted_discharge)

# Compute residuals
residuals <- -(sweep(mu, 1, df$Discharge, FUN = "-"))
transformed_residuals1 <- sweep(residuals, 1, remnant_sigma, FUN = "/")

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
      Discharge = mean_predicted_discharge
    ),
    aes(x = Stage, y = Discharge), 
    color = 'red', 
    size = 2,          # Increased from 2
    alpha = 0.8
  ) +
  
  # Labels and titles
  labs(
    title = "Pathagudem",
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


# Define station name and desired location
station_name <- "Pathagudem"  # Change this to the actual station name
save_path <- file.path("C:/Users/kona1280/OneDrive - Indian Institute of Science/RESEARCH FILES/SECOND WORK/STAN_RESULTS", paste0(station_name, ".RData"))  # Change path as needed

# Save objects to workspace
#save(df,  fit, residual_df, envelope_df, new_summary_df, file = save_path)

# Confirm the file was saved
cat("Workspace saved to:", save_path, "\n")

library(openxlsx)
# Save the dataframe as an Excel file
write.xlsx(predicted_discharge, 
           file = file.path("C:/Users/kona1280/OneDrive - Indian Institute of Science/RESEARCH FILES/SECOND WORK/STAN_RESULTS", paste0(station_name, "_predicted_Q_full.xlsx")), 
           overwrite = TRUE)
# Confirm the file was saved
cat("Workspace saved to:", save_path, "\n")