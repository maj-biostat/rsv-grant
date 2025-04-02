// Log-logistic AFT model
data {
  int<lower=0> N;             // Number of observations
  int<lower=0> P;             // Number of predictors
  matrix[N, P] X;             // Predictor matrix X[, 1] is intercept
  vector<lower=0>[N] y;       // Observed survival times
  vector<lower=0, upper=1>[N] event;  // Event indicator (1=event, 0=censored)

  // int N_pred;
  // vector[N_pred] t_surv;    // time to predict survival at

  // prior
  vector[2] mu0_gamma;   // e.g. c(5, 0)
  vector[2] sd0_gamma;   // e.g. c(2, 2)

  real rho_shape;        // e.g. 0.5 to 1

}

parameters {
  vector[P] gamma;             // Regression coefficients for scale
  real<lower=0> shape;        // Shape parameter (b in the formula)
}

transformed parameters {
  // Location parameter (log-scale)
  vector[N] mu;

  mu = X * gamma;
}

model {
  // Priors - arbitrary at the moment
  target += normal_lpdf(gamma[1] | mu0_gamma[1], sd0_gamma[1]);
  target += normal_lpdf(gamma[2] | mu0_gamma[2], sd0_gamma[2]);
  target += exponential_lpdf(shape | rho_shape);

  // Likelihood
  for (i in 1:N) {
    if (event[i] == 1) {
      // For observed events, use the log-logistic density
      target += log(shape) - mu[i] + (shape - 1) * (log(y[i]) - mu[i]) -
                2 * log1p(pow(y[i] / exp(mu[i]), shape));
    } else {
      // For censored observations, use the log survival function
      target += -log1p(pow(y[i] / exp(mu[i]), shape));
    }
  }
}

generated quantities {
  vector[N_pred] surv0;
  vector[N_pred] surv1;

  real<lower=0, upper=1> p0_360;
  real<lower=0, upper=1> p1_360;

  real scale0;
  real scale1;

  // these equate to the median survival time
  scale0 = exp(gamma[1]);
  scale1 = exp(gamma[1] + gamma[2]);

  p0_360 = 1 - (1 / (1 + pow(360.0/scale0,  shape)));
  p1_360 = 1 - (1 / (1 + pow(360.0/scale1,  shape)));

  // for(i in 1:N_pred){
  //   surv0[i] =  1 / (1 + pow(t_surv[i]/scale0,  shape));
  //   surv1[i] =  1 / (1 + pow(t_surv[i]/scale1,  shape));
  // }


}


