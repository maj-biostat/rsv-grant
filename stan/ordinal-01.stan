data {
  int<lower=1> N;
  int<lower=1> K;
  int<lower=1> P; // num pars in linear predictor (no intercept)
  array[N] int<lower=1,upper=K> y;
  matrix[N,P] X;

}
parameters {
  ordered[K-1] cuts;
  vector[P] b;
}
model {
  target += normal_lpdf(b | 0, 2);
  target += normal_lpdf(cuts | 0, 2);


  target += ordered_logistic_lpmf(y | X * b, cuts);
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
