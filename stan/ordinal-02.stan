data {
  int<lower=1> N;
  int<lower=1> K;
  int<lower=1> P; // num pars in linear predictor (no intercept)
  array[N] int<lower=1,upper=K> y;
  matrix[N,P] X;
  vector[N] wgt;
  //
  real pri_b_s;
  real pri_cuts_s;
}
parameters {
  // ordered[K-1] cuts;
  simplex[K] c;
  vector[P] b;
}
transformed parameters{
  ordered[K-1] cuts = logit(cumulative_sum(c[1:(K-1)]));
}
model {
  target += normal_lpdf(b | 0, pri_b_s);
  // target += normal_lpdf(cuts | 0, pri_cuts_s);
  target += dirichlet_lpdf(c | rep_vector(1, K));

  for(i in 1:N){
    target += wgt[i] * ordered_logistic_lpmf(y[i] | X[i, ] * b, cuts);
  }
}
generated quantities {
  vector[K] p0;
  vector[K] p1;
  vector[K] rd;

  p0[1] = 1 - inv_logit(- cuts[1]);
  p1[1] = 1 - inv_logit(b[1] - cuts[1]);
  p0[K] = inv_logit(- cuts[K-1]);
  p1[K] = inv_logit(b[1] - cuts[K-1]);

  for(k in 2:(K-1)){
    p0[k] = inv_logit(- cuts[k-1]) - inv_logit(- cuts[k]);
    p1[k] = inv_logit(b[1] - cuts[k-1]) - inv_logit(b[1] - cuts[k]);
  }

  for(k in 1:K){
    rd[k] = p1[k] - p0[k];
  }
}
