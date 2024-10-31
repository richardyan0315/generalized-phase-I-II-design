data{
  int<lower=0> l;
  array[l] int<lower=0> n;
  array[l] real sum_qtp;
  array[l] real x;
  real a;
  real b;
  real<lower=0> c;
}

parameters{
  real alpha;
  real<lower=0> beta;
}

transformed parameters{
  vector[l] log_lik;
  for (i in 1:l){
    log_lik[i] = - sum_qtp[i] * log(1 + exp(-alpha - beta * x[i])) - 
      (n[i] - sum_qtp[i])* log(1 + exp(alpha + beta * x[i]));
  }
}

model{
  alpha ~ uniform(a, b);
  beta ~ uniform(0, c);
  target += uniform_lpdf(alpha| a, b);
  target += uniform_lpdf(beta|0, c);
  target += sum(log_lik);
}
