//
// https://experienced-sampler.netlify.app/post/stan-hierarchical-ar/#testing
//

// The input data is a vector 'y' of length 'N', we I individuals who we 
// measured T times for N (I*T) observations in total. We also have an indicator 
// variable that illustrates what individual (1,..., I) data belong to, and what 
// measurement occasion (1,..., T) the data was collected at. 
// Data is in long-format
data {
  int<lower=0> N;
  int<lower=0> I;
  int<lower=0> T;
  array[N] int<lower=1, upper=I> individual;
  array[N] int<lower=1, upper=T> time;
  vector[N] y;
}
transformed data {
  vector[N] y_std;
  real meanY;
  real sdY;
  meanY = mean(y);
  sdY = sd(y);
  y_std = (y - meanY) / sdY;
}
// The parameters accepted by the model. 
parameters {
  real alpha_hat;
  real<lower=0> alpha_scale;
  real<lower=-1, upper=1> beta_hat;
  real<lower=0> beta_scale;
  vector[I] alpha;
  vector<lower=-1, upper=1>[I] beta;
  real<lower=0> sigma;
}
// The model to be estimated. We model the output 'y_std[n]' to be normally 
// distributed with mean 'alpha[n] + beta[n] * y_c[n-1]' and standard deviation
// 'sigma'. We use the group-mean centered values of y as predictors so that 
// alpha gives us individual means instead of intercepts.
model {
  vector[N] y_c;
  
  alpha_hat ~ normal(0, 5);
  beta_hat ~ normal(0, .5);
  alpha_scale ~ normal(0, 1);
  beta_scale ~ normal(0, 1);
  sigma ~ normal(0, 2);
  
  for (i in 1 : I) {
    alpha[i] ~ normal(alpha_hat, alpha_scale);
    beta[i] ~ normal(beta_hat, beta_scale);
  }
  
  y_c[1] = y_std[1] - alpha[individual[1]];
  
  for (n in 2 : N) {
    y_c[n] = y_std[n] - alpha[individual[n]];
    if (time[n] > 1) 
      y_std[n] ~ normal(alpha[individual[n]]
                        + beta[individual[n]] * y_c[n - 1], sigma);
  }
}
generated quantities {
  vector[I] alphas_ind;
  real alpha_hat_raw;
  real<lower=0> alpha_scale_raw;
  real<lower=0> sigma_raw;
  alphas_ind = (sdY * alpha) + meanY;
  alpha_hat_raw = (sdY * alpha_hat) + meanY;
  alpha_scale_raw = sdY * alpha_scale;
  sigma_raw = sigma * sdY;
}


