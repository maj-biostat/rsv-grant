data {
  int N;
  array[N] int n;
  array[N] int y;
  matrix[N, 1] X;
}
transformed data {
}
parameters {
  real alpha;
  vector[1] b;
}
transformed parameters{
  real theta = exp(b[1]);
}
model {
  target += logistic_lpdf(alpha | 0, 1);
  target += normal_lpdf(b | 0, 2);
  target += binomial_logit_glm_lpmf(y | n, X, alpha, b);
}
generated quantities {
  // makes assumptions about order of data in n and y and X.
  real p0 = inv_logit(alpha);
  real p1 = inv_logit(alpha + b[1]);
  real rd = p1 - p0;
}
