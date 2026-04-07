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
  real<lower=0> gamma1;
  real<lower=0> gamma2; // Parameter for remnant sigma
  real<lower=0> x; // Exponent for remnant sigma
  real<lower=0> lambda; // Scaling factor for smooth transition
}

transformed parameters {
  real b1 = k1; // Continuity condition for segment 1

  // Adjusted continuity condition for segment 2
  real b2 = k2 - lambda * pow(fmax(1e-6, (a1 / a2) * pow(fmax(1e-6, k2 - b1), c1)), 1 / c2);
}

model {
  vector[N] mu;             // Predicted discharge
  vector[N] remnant_sigma;  // Remnant sigma
  vector[N] sigma;          // Total sigma

  // Priors based on values that worked
  a1 ~ normal(95, 10);
  c1 ~ normal(1.8, 0.25);
  k1 ~ normal(2, 0.5);
  k2 ~ normal(8, 1);
  a2 ~ normal(25, 5);
  c2 ~ normal(2.5, 0.5);
  gamma1 ~ normal(0.2, 0.1);
  gamma2 ~ normal(8, 3);
  x ~ normal(0.5, 0.1);

  // Prior for lambda
  lambda ~ normal(0.1, 0.05); // Encourages small values around 0.1

  // Vectorized discharge calculation
  for (i in 1:N) {
    real seg1 = a1 * pow(fmax(0, h[i] - b1), c1);
    real seg2 = a2 * pow(fmax(0, h[i] - b2), c2);

    if (h[i] <= k1) {
      mu[i] = 0;
    } else if (h[i] <= k2) {
      mu[i] = seg1;
    } else {
      mu[i] = seg1 + seg2; // Sum both segments when h > k2
    }

    // Remnant sigma and total sigma calculations
    remnant_sigma[i] = gamma1 + gamma2 * pow(mu[i], x);
    sigma[i] = sqrt(pow(remnant_sigma[i], 2) + pow(0.07 * mu[i], 2));
  }

  // Likelihood
  Q ~ normal(mu, sigma);
}
