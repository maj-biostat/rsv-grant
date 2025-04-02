data {

  int n0;
  int y0;

  int nt;
  int yt;

  // historic controls
  int H;
  array[H] int yh;
  array[H] int nh;

}
transformed data {
}
parameters {
  real mu;
  real<lower=0> tau;
  real phi;
  real gamma0;
  vector[H] gamma;
}
transformed parameters{
  real p0;
  real pt;
  vector[H] ph;
  real theta;

  p0 = inv_logit(gamma0);
  for(i in 1:H){
    ph[i] = inv_logit(gamma[i]);
  }
  pt = inv_logit(gamma0 + phi);
  theta = exp(phi);
}
model {
  target += logistic_lpdf(mu | 0, 1);
  target += exponential_lpdf(tau | 1);

  target += normal_lpdf(gamma0 | mu, tau);
  target += normal_lpdf(gamma | mu, tau);

  target += normal_lpdf(phi | 0, 2);

  // likelihood
  target += binomial_lpmf(y0 | n0, p0);
  for(i in 1:H){
    target += binomial_lpmf(yh[i] | nh[i], ph[i]);
  }
  target += binomial_lpmf(yt | nt, pt);
}


