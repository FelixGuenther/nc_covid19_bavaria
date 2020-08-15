data {
  // No days
  int T;
  // D
  int maxDelay;
  // Number cps delay dist
  int n_cp;
  // Number additional Covariate effects delay distribution
  int n_wextra;
  // observation triangle
  int rT[T, maxDelay+1];
  // Design matrix delay dist
  matrix[maxDelay+1, n_cp + n_wextra] W[T];
  // Lambda Random Walk Params
  // Prior Delay gammas (log-OR Haz)
  real mu_gamma[maxDelay];
  real sd_gamma[maxDelay];
  // Prior Effects Delay dist (Changepoints and Wextra)
  real eta_mu[n_cp + n_wextra];
  real eta_sd[n_cp + n_wextra];
  // Predictions up to predLag days before now
  int predLag;
}
parameters {
  real logLambda[T];
  real gamma[maxDelay];
  vector[n_cp + n_wextra] eta;
  real<lower=0> sd_logLambda;
  real<lower=0> phi_negbinom;
}
transformed parameters {
  real p[T, maxDelay+1];
  real haz[T, maxDelay+1];
  
  for(t in 1:T) {
    //delay distribution
    haz[t, 1] = inv_logit(gamma[1] + W[t,1] * eta);
    p[t, 1] = inv_logit(gamma[1] + W[t,1] * eta);
    for(d in 1:(maxDelay-1)) {
      haz[t, d+1] = inv_logit(gamma[d+1] + W[t,d+1] * eta);
      p[t, d+1] = (1-sum(p[t, 1:d])) * haz[t, d+1];
    }
    haz[t, maxDelay+1] = 1;
    p[t,maxDelay+1] = 1 - sum(p[t, 1:maxDelay]);
  }
}
model {
  // priors
  // logLambda
  logLambda[1] ~ normal(0,3);
  sd_logLambda ~ normal(0,5); // half normal due to constraint
  for (t in 2:T) {
    logLambda[t] ~ normal(logLambda[t-1], sd_logLambda);
  }
  // Delay dist
  //coefs for logit @ delay 0,..,maxDelay-1
  gamma ~ normal(mu_gamma, sd_gamma);
  // Changepoint effects
  eta ~ normal(eta_mu, eta_sd);
  //phi_negbinom ~ inv_gamma(0.01, 0.01);

  
  // Likelihood
  for(t in 1:T) {
    for (d in 0:min(T-t,maxDelay)) {
      rT[t,d+1] ~ neg_binomial_2(exp(logLambda[t])*p[t,d+1]+0.1, phi_negbinom);
    }
  }
}
generated quantities {
  int nT[maxDelay-predLag, maxDelay+1];
  real ntInf[maxDelay-predLag];
  for(t in (T-maxDelay+1):(T-predLag)) {
    for (d in 0:maxDelay) {
       if((t+d) <= T) {
       nT[t-(T-maxDelay), d+1] = rT[t, d+1];
       } else {
        nT[t-(T-maxDelay),d+1] = neg_binomial_2_rng(exp(logLambda[t])*p[t,d+1]+0.1, phi_negbinom);
       }
    }
    //ntInf[t-(T-maxDelay)] = sum(nT[t-(T-maxDelay),]) - (t-(T-maxDelay)) * 0.1;
    ntInf[t-(T-maxDelay)] = sum(nT[t-(T-maxDelay),]);
  }
}
