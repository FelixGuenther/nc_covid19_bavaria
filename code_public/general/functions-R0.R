## nc is a nowcast object as created by
sample_from_nc_posterior <- function(nc, n_sample = 100) {

  mcmc_samples <- attr(nc@delayCDF[["bayes.trunc.ddcp"]], "model")$mcmc_samples
  ## sample n_sample rows, each containing on time-series of Ntinf from each chain
  Ntinf <- map(mcmc_samples, ~{
    dn2 <- dimnames(.x)[[2]]
    ind_ntinf <- grep("NtInf", dn2)
    .x[sample(seq_len(nrow(.x)), n_sample), ind_ntinf]

  })

  do.call(rbind, Ntinf)

}

tidy_Rt <- function(R0.R) {

  data.frame(
    Date = as.Date(names(R0.R$R)),
    Rt = R0.R$R,
    Rt_lower = R0.R$conf.int$lower,
    Rt_upper = R0.R$conf.int$upper)
}

tidy_Rsimu <- function(R0.R, nc, start, end) {

  epochs <- epoch(nc)
  ind <- which(epochs >= start & epochs <= end)

  R_simu_df <- R0.R$R.simu %>% as.data.frame()
  R_simu_df <- R_simu_df[ind, , drop = FALSE]
  R_simu_df$Date <- as.Date(epochs[ind])
  R_simu_df

}


sample_from_nc_posterior2 <- function(nc, n_sample = 100) {

  mcmc_samples <- attr(nc@delayCDF[["bayes.trunc.ddcp"]], "model")$mcmc_samples

  # Extract the logLambda samples
  names <- colnames(mcmc_samples[[1]])
  idx <- grep("logLambda\\[[0-9]+\\]", names)
  logLambda <- do.call(rbind, mcmc_samples[, idx])

  idx_to_keep <- sample(seq_len(nrow(logLambda)), size = n_sample, replace = FALSE)
  logLambda_subset <- logLambda[idx_to_keep,]

  # Sample y_t ~ Po( exp(logLambda_t))
  ts <- apply(logLambda_subset, 2, function(log_lambda) rpois(length(log_lambda), lambda=exp(log_lambda)))

  #Done
  return(ts)

}
