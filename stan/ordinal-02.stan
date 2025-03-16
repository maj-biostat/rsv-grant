data {
  int<lower=1> N;
  int<lower=1> K;
  int<lower=1> P; // num pars in linear predictor (no intercept)
  array[N] int<lower=1,upper=K> y;
  matrix[N,P] X;
  vector[N] wgt;
}
parameters {
  ordered[K-1] cuts;
  vector[P] b;
}
model {
  target += normal_lpdf(b | 0, 1);
  target += normal_lpdf(cuts | 0, 5);

  for(i in 1:N){
    target += wgt[i] * ordered_logistic_lpmf(y[i] | X[i, ] * b, cuts);
  }


}
generated quantities {
  real p_MA_rti_rsv_soc = inv_logit(cuts[3]);
  real p_MA_rti_rsv_new = inv_logit(cuts[3] + b[1]);
  real rd = p_MA_rti_rsv_new - p_MA_rti_rsv_soc;
}
