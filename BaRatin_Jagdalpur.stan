data {
  int<lower=0> N;       // Number of observations
  vector[N] h;          // Observed stage
  vector[N] Q;          // Observed discharge
}

parameters {
  real<lower=0> a1;     // Scaling parameter for segment 1
  real<lower=0> c1;     // Power law exponent for segment 1
  real<lower=0> k1;     // Threshold parameter for segment 1
  real<lower=0> k2;     // Threshold parameter for segment 2
  real<lower=0> a2;     // Scaling parameter for segment 2
  real<lower=0> c2;     // Power law exponent for segment 2
  real <lower=0> gamma1;          // Parameter for remnant sigma
  real <lower=0> gamma2;          // Parameter for remnant sigma
  real<lower=0> x;      // Exponent for remnant sigma
}

transformed parameters {
  real b1 = k1;         // Continuity condition for b1: b1 = k1
  real b2 = k2 - pow((1 / a2) * a1 * pow(k2 - b1, c1), 1 / c2); // Continuity condition for b2
}

model {
  vector[N] mu;             // Predicted discharge
  vector[N] remnant_sigma;  // Remnant sigma
  vector[N] sigma;          // Total sigma

  // Priors
  a1 ~ normal(13, 5);  // Prior for a1
  c1 ~ normal(2, 1);    // Prior for c1
  k1 ~ normal(0.5, 0.1);       // Prior for k1
  k2 ~ normal(1.7, 0.5);     // Prior for k2
  a2 ~ normal(30, 10);    // Prior for a2
  c2 ~ normal(1.7, 0.3);    // Prior for c2
  gamma1 ~ uniform(0.01, 50);  // Prior for gamma2
  gamma2 ~ normal(0.1, 0.1);  // Prior for gamma2
  x ~ normal(0.45, 0.1);       // Prior for exponent x

  // Calculate mu, remnant_sigma, and sigma
  for (i in 1:N) {
         mu[i] = if_else(h[i] <= k1, 0,
         if_else(h[i] <= k2, a1 * pow(fmax(0, h[i] - b1), c1),
         a2 * pow(fmax(0, h[i] - b2), c2)));


    remnant_sigma[i] = gamma1 + (gamma2 * pow(mu[i], x));
    sigma[i] = sqrt(remnant_sigma[i]^2 + pow(0.07 * mu[i], 2));
  }

  // Likelihood
  Q ~ normal(mu, sigma);
}
