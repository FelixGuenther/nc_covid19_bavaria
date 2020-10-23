data {
  // No days
  int T;
  // D
  int maxDelay;
  // observation triangle
  int rT[T, maxDelay+1];
  // Prior Delay gammas (log-OR Haz)
  real mu_gamma[maxDelay];
  real sd_gamma[maxDelay];
  // Delay rW
  int rW_params;
  matrix[T, rW_params] rW_delay_des;
  // Predictions up to predLag days before now
  int predLag;
  // True data for scoring rules
  int nT_true[maxDelay-predLag, maxDelay+1];
}
parameters {
  //logLambda
  real logLambda[T];
  real<lower=0> sd_logLambda;
  //delay Distribution
  real gamma[maxDelay];
  real<lower=0> phi_negbinom;
  vector[rW_params] alpha_del;
  real<lower=0> sd_alpha_del;
}
transformed parameters {
  real p[T, maxDelay+1];
  real haz[T, maxDelay+1];
  
  haz[1, 1] = inv_logit(gamma[1]);
  p[1, 1] = inv_logit(gamma[1]);
  for(d in 1:(maxDelay-1)) {
    haz[1, d+1] = inv_logit(gamma[d+1]);
    p[1, d+1] = (1-sum(p[1, 1:d])) * haz[1, d+1];
  }
  haz[1, maxDelay+1] = 1;
  p[1,maxDelay+1] = 1 - sum(p[1, 1:maxDelay]);

  for(t in 2:T) {
    //delay distribution
    haz[t, 1] = inv_logit(gamma[1] + dot_product(rW_delay_des[t,], alpha_del));
    p[t, 1] = inv_logit(gamma[1] + dot_product(rW_delay_des[t,], alpha_del));
    for(d in 1:(maxDelay-1)) {
      haz[t, d+1] = inv_logit(gamma[d+1] + dot_product(rW_delay_des[t,], alpha_del));
      p[t, d+1] = (1-sum(p[t, 1:d])) * haz[t, d+1];
    }
    haz[t, maxDelay+1] = 1;
    p[t,maxDelay+1] = 1 - sum(p[t, 1:maxDelay]);
  }
}
model {
  // priors
  // logLambda
  logLambda[1] ~ normal(0,1);
  sd_logLambda ~ normal(0,0.5); // half normal due to constraint
  // Delay dist
  //coefs for logit @ delay 0,..,maxDelay-1
  gamma ~ normal(mu_gamma, sd_gamma);
  // First order random walk on delay dist
  alpha_del[1] ~ normal(0,.1);
  sd_alpha_del ~ normal(0,0.1);

  for (t in 2:T) {
    logLambda[t] ~ normal(logLambda[t-1], sd_logLambda);
  }
  for (i in 2:rW_params) {
    alpha_del[i] ~ normal(alpha_del[i-1], sd_alpha_del);
  }

  // Likelihood
  for(t in 1:T) {
    for (d in 0:min(T-t,maxDelay)) {
      rT[t,d+1] ~ neg_binomial_2(exp(logLambda[t])*p[t,d+1]+0.01, phi_negbinom);
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
        nT[t-(T-maxDelay),d+1] = neg_binomial_2_rng(exp(logLambda[t])*p[t,d+1]+0.01, phi_negbinom);
       }
    }
    ntInf[t-(T-maxDelay)] = sum(nT[t-(T-maxDelay),]) - round((t-(T-maxDelay)) * 0.01);
  }
}

