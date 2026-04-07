data {
  int<lower=0> N; // Number of observations
  vector[N] h;    // Observed stage
  vector[N] Q;    // Observed discharge
}

parameters {
  real<lower=0> a1; // Scaling parameter for segment 1
  real<lower=0> c1; // Power law exponent for segment 1
  real<lower=0> k1; // Threshold parameter for segment 1
  real<lower=0> k2; // Threshold parameter for segment 2
  real<lower=0> a2; // Scaling parameter for segment 2
  real<lower=0> c2; // Power law exponent for segment 2
  real gamma1;       // Parameter for remnant sigma
  real<lower=0> gamma2; // Parameter for remnant sigma
  real<lower=0,upper=1> x; // Exponent for remnant sigma
}

transformed parameters {
  real b1 = k1; // b1 continuity condition: b1 = k1
  real b2 = k2 - pow((1/a2) * a1 * pow(k2 - b1, c1), 1/c2); // Continuity condition for b2
  vector[N] mu;  // Predicted discharge
  vector[N] remnant_sigma;  // Remnant sigma
  vector[N] sigma;  // Total sigma
  vector[N] tau;    // Precision (1/sigma^2)

  for (i in 1:N) {
    // Use if-else to select the correct mu[i] based on the value of h[i]
    if (h[i] <= k1) {
      mu[i] = 0;
    } else if (h[i] <= k2) {
      mu[i] = a1 * pow(fmax(0, h[i] - b1), c1);
    } else {
      mu[i] = a2 * pow(fmax(0, h[i] - b2), c2);
    }

    // Compute the remnant sigma[i] using gamma1, gamma2, and the exponent x
    remnant_sigma[i] = gamma1 + gamma2 * pow(mu[i], x);

    // New sigma[i] calculation as the square root of the sum of squares
    sigma[i] = sqrt(square(remnant_sigma[i]) + square(0.07 * mu[i]));

    // Precision as inverse of variance (tau[i] = 1/sigma[i]^2)
    tau[i] = 1 / square(sigma[i]);
  }
}

model {
  // Priors
  a1 ~ normal(10.63, 0.1);
  c1 ~ normal(2.67, 0.16);
  k1 ~ normal(1.5, 0.04);
  k2 ~ normal(7.0, 0.04);
  a2 ~ normal(150.0, 0.04);
  c2 ~ normal(2.0, 0.16);
  gamma1 ~ normal(0.0, 0.01);
  gamma2 ~ uniform(0, 100);
  x ~ normal(0.7, 0.04);

  // Likelihood
  for (i in 1:N) {
    Q[i] ~ normal(mu[i], tau[i]);
  }
}
