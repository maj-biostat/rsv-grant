data {

  int n0;
  int y0;

  int nt;
  int yt;

  // historic controls
  int H;
  array[H] int yh;
  array[H] int nh;

  // prior hyperparameters:
  // distribution of control arm response based on current and historic data
  // central value (mean, sd)
  vector[2] pri_mu;
  // variation
  real<lower=0> pri_tau_rho;

  // hyperparameters for treatment arm
  vector[2] pri_pt;
}
transformed data {
}
parameters {
  // control arm parameters (distribution of historic controls)
  // central value for response in control arm
  real mu;
  // variation in distribution of controls
  real<lower=0> tau;
  // log-odds response for current and historic controls
  real gamma0;
  vector[H] gamma;
  // risk on intervention arm
  real<lower=0, upper=1> pt;
}
transformed parameters{
  real p0;
  vector[H] ph;
  real rd;
  real ve;

  p0 = inv_logit(gamma0);
  for(i in 1:H){
    ph[i] = inv_logit(gamma[i]);
  }

  // risk on treatment less risk on control
  rd = pt - p0;
  ve = 1 - (pt/p0);
}
model {
  target += normal_lpdf(mu | pri_mu[1], pri_mu[2]);
  target += exponential_lpdf(tau | pri_tau_rho);

  // log-odds response for current and historic controls arise from
  // common distribution
  target += normal_lpdf(gamma0 | mu, tau);
  target += normal_lpdf(gamma | mu, tau);

  // independently estimate pt
  target += beta_lpdf(pt | pri_pt[1], pri_pt[2]);

  // likelihood
  target += binomial_lpmf(y0 | n0, p0);
  for(i in 1:H){
    target += binomial_lpmf(yh[i] | nh[i], ph[i]);
  }
  target += binomial_lpmf(yt | nt, pt);
}


