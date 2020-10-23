data {
  int T;
  int maxDelay;
  int rT[T, maxDelay+1];
  real mu_gamma[maxDelay];
  real sd_gamma[maxDelay];
  int predLag;
}
parameters {
  real logLambda[T];
  real gamma[maxDelay];
  real<lower=0> sd_logLambda;
}
transformed parameters {
  real p[T, maxDelay+1];
  real haz[T, maxDelay+1];
  
  for(t in 1:T) {
    //delay distribution
    haz[t, 1] = inv_logit(gamma[1]);
    p[t, 1] = inv_logit(gamma[1]);
    for(d in 1:(maxDelay-1)) {
      haz[t, d+1] = inv_logit(gamma[d+1]);
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
  
  // Likelihood
  for(t in 1:T) {
    for (d in 0:min(T-t,maxDelay)) {
      rT[t,d+1] ~ poisson(exp(logLambda[t])*p[t,d+1]+0.01);
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
         nT[t-(T-maxDelay),d+1] = poisson_rng(exp(logLambda[t])*p[t,d+1]+0.01);
       }
    }
    ntInf[t-(T-maxDelay)] = sum(nT[t-(T-maxDelay),]) - round((t-(T-maxDelay)) * 0.01);
  }
}

