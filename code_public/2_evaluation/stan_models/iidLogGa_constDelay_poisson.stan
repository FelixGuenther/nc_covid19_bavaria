data {
  int T;
  int maxDelay;
  int rT[T, maxDelay+1];
  real<lower=0> alpha_lambda;
  real<lower=0> beta_lambda;
  real mu_gamma[maxDelay];
  real sd_gamma[maxDelay];
  int predLag;
}
parameters {
  real<lower=0> lambda[T];
  real gamma[maxDelay];
}
transformed parameters {
  real p[T, maxDelay+1];
  real haz[T, maxDelay+1];
  real logLambda[T];
  
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
    //logLambda
    logLambda[t] = log(lambda[t]);
  }
}
model {
  // priors
  lambda ~ gamma(alpha_lambda, beta_lambda);
  //coefs for logit @ delay 0,..,maxDelay-1
  gamma ~ normal(mu_gamma, sd_gamma);
  
  // likelihood
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

