data {
    int<lower=0> N;          // Number of observations
    vector[N] h;             // Observed stage
    vector[N] Q;             // Observed discharge
   //vector[N] weights;  // Added weights
}

parameters {
    real<lower=0> a1;        // Scaling parameter
    real<lower=0> c1;        // Power law exponent
    real<lower=0, upper=1.09> k1; // Threshold parameter
    real<lower=0> gamma1;             // Parameter for remnant sigma
    real<lower=0> gamma2;             // Parameter for remnant sigma
    real<lower=0> x;         // Exponent for remnant sigma
}

transformed parameters {
    real b1 = k1;            // Continuity condition for b1
}

model {
    vector[N] mu;            // Predicted discharge
    vector[N] remnant_sigma; // Remnant sigma
    vector[N] sigma;         // Total sigma

    // Priors
    a1 ~ normal(250, 50);      // Prior for a1
    c1 ~ normal(2.7, 0.05);     // Prior for c1
    k1 ~ normal(1.05, 1);      // Prior for k1
    gamma1 ~ uniform(0, 50);   // Prior for gamma2
    gamma2 ~ uniform(0.1, 50);   // Prior for gamma2
    x ~ normal(0.5, 0.3);        // Prior for exponent x

    // Calculate mu, remnant_sigma, and sigma
    for (i in 1:N) {
        if (h[i] <= k1) {
            mu[i] = 0;
        } else {
            mu[i] = a1 * pow((h[i] - b1), c1);
        }
        remnant_sigma[i] = gamma1 + (gamma2 * pow(mu[i], x)); // , 1
        sigma[i] = sqrt(remnant_sigma[i]^2 + pow(0.07 * mu[i], 2));
    }
    // Weighted Likelihood
 //  for (i in 1:N) {
  //  target += weights[i] * normal_lpdf(Q[i] | mu[i], sigma[i]);
  // }
    // Likelihood
    Q ~ normal(mu, sigma);
}
