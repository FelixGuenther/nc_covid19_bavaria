# Load function to estimate Rt
source("../../1_main_analysis/general/est.R0.TD.R") # customized version of R0::est.R0.TD

#' Create a list of objects that are required as data input for different nowcast models
#' in STAN
#' @param dat person-specific dataset with columns rep_date (reporting date), 
#' disease_start (symptom onset date)
#' @param now current day for estimating nowcast
#' @param begin_date date to start estimation
#' @para D maximum considered delay
#' @param predLag difference between `now` and the last desired prediction of 
#' number of disease onsets
#' @param mod_name name of the nowcast model (to return the required data objects). 
#' one of: "negBinom_rW_rW", "poisson_iid_consDel", "negBinom_rW_cp2W_wd", 
#' "negBinom_rW_rW_2d", "negBinom_rW_rW_wd", "negBinom_rW_cp2W", 
#' "poisson_rW_cp2W", "poisson_rW_const", "negBinom_rW_cptrue", 
#' "poisson_rW_cptrue"
#' @return list of preprocessed data objects/priors

prepare_data_list = function(dat,
                             now,
                             begin_date = ymd("2020-02-24"),
                             D = 21,
                             predLag = 2,
                             mod_name = "negBinom_rW_rW") {
  dat_full = dat
  dat = dat %>% filter(rep_date<= now)
  so_range <- range(begin_date, now)

  sts <- linelist2sts(
    as.data.frame(dat),
    dateCol="disease_start",
    aggregate.by = "1 day",
    dRange = so_range)

  dMin <- begin_date
  dMax <- max(dat$disease_start, na.rm = TRUE)

  # Reporting triangle
  cat("Building reporting triangle...\n")
  t0 <- dMin
  t02s <- seq(t0, now, by = "1 day")
  t02s_full = seq(t0, max(dat_full$rep_date)-1, by = "1 day")
  T <- length(t02s) - 1
  T_full = length(t02s_full) - 1 
  n <- matrix(NA, nrow = T + 1, ncol = T + 1, dimnames = list(as.character(t02s),
                                                              NULL))
  n_full = matrix(NA, nrow = T + 1, ncol = T_full + 1, dimnames = list(as.character(t02s),
                                                                  NULL))
  timeDelay <- function(d1, d2) {
    as.numeric(d2 - d1)
  }
  dat = dat %>% mutate(delay = rep_date - disease_start)
  dat_full = dat_full %>% mutate(delay = rep_date - disease_start)
  for (t in 0:T) {
    data.att <- dat[which(dat$disease_start == t02s[t+1]), ]
    data.att.full <- dat_full[which(dat_full$disease_start == t02s[t+1]), ]
    for (x in 0:(T - t)) {
      n[t + 1, x + 1] <- sum(data.att$delay == x)
    }
    for(x in 0:(T_full-t)) {
      n_full[t + 1, x + 1] <- sum(data.att.full$delay == x)
    }
  }
  cat("No. cases: ", sum(n, na.rm = TRUE), "\n")
  nLongDelay <- apply(n[, (D + 1) + seq_len(T - D), drop = FALSE],
                      1, sum, na.rm = TRUE)
  if (any(nLongDelay > 0)) {
    warning(paste(sum(nLongDelay), " cases with a delay longer than D=",
                  D, " days forced to have a delay of D days.\n", sep = ""))
    n <- n[, 1:(D + 1)]
    n[, (D + 1)] <- n[, (D + 1)] + nLongDelay
  } else {
    n <- n[, 1:(D + 1)]
  }
  # Replace unobserved dates with 0
  n[is.na(n)] <- 0
  
  # Full data
  nLongDelay_full <- apply(n_full[, (D + 1) + seq_len(T_full - D), drop = FALSE],
                      1, sum, na.rm = TRUE)
  if (any(nLongDelay_full > 0)) {
    n_full <- n_full[, 1:(D + 1)]
    n_full[, (D + 1)] <- n_full[, (D + 1)] + nLongDelay_full
    } else {
    n_full <- n_full[, 1:(D + 1)]
  }
  # Replace unobserved dates with 0
  n_full[is.na(n_full)] <- 0
  
  # Priors
  expectation_p_haz_plogisN <- function(mu, sigma, p_smaller) {
    f <- function(x) {
      (1-p_smaller) * plogis(x) * dnorm(x, mu, sigma)
    }
    int <- integrate(f, lower = -Inf, upper = Inf)
    int$value
  }
  # Variance of P(D>=d)*P(D=d|D>=d), with P(D=d|D>=d)=1/(1+exp(-x)) and X~N(mu, sigma^2)
  var_p_haz_plogisN <- function(mu, sigma, p_smaller) {
    E = expectation_p_haz_plogisN(mu, sigma, p_smaller)
    f <- function(x) {
      (((1-p_smaller) * plogis(x)) - E)^2 * dnorm(x, mu, sigma)
    }
    int <- integrate(f, lower = -Inf, upper = Inf)
    int$value
  }
  get_prior_gamma <- function(gd.prior.kappa = 1, D = 20) {
    alpha_i = gd.prior.kappa
    alpha_0 = (D+1) * gd.prior.kappa
    e_dir = alpha_i/alpha_0
    var_dir = (alpha_i*(alpha_0-alpha_i)) / (alpha_0^2*(alpha_0+1))

    sum_sq_diff <- function(theta, p_smaller) {
      mu <- theta[1]
      sigma <- theta[2]
      E_diff <- expectation_p_haz_plogisN(mu, sigma, p_smaller) - e_dir
      V_diff <- var_p_haz_plogisN(mu, sigma, p_smaller) -  var_dir
      E_diff^2 + V_diff^2
    }

    mu_gamma = rep(0, D)
    sigma_gamma = rep(0, D)
    for(i in 1:D) {
      theta_start = if(i == 1) {
        c(-3,1)
      } else {
        c(mu_gamma[i-1], sigma_gamma[i-1])
      }
      optim_res = optim(theta_start, function(theta) {
        sum_sq_diff(theta, p_smaller = (i-1)*e_dir)
      })
      mu_gamma[i] = ifelse(optim_res$convergence==0, optim_res$par[1], NA)
      sigma_gamma[i] = ifelse(optim_res$convergence==0, optim_res$par[2], NA)
    }
    data.frame(mu_gamma = mu_gamma, sigma_gamma = sigma_gamma)
  }

  # Prior for logit hazards
  gamma_prior_kappa = get_prior_gamma(gd.prior.kappa = 5, D=D)

  # Prior for iidLogGa
  exp = mean(observed(sts))
  var = 5*var(observed(sts))

  dslnex <- function(x, aim = c(exp, var)) {
    y = numeric(2)
    y[1] <- (x[1] * x[2]) -  aim[1]
    y[2] <- (x[1] * x[2]^2) - aim[2]
    y
  }

  x_start = rep(sqrt(exp), 2)
  res = nleqslv(x_start, dslnex, control=list(btol=.01))
  beta.lambda = 1/res$x[2]
  alpha.lambda = res$x[1]

  # W matrix
  ## IID Log Ga cps in delay dist
  ddChangepoint <- sort(seq(dMax, dMin, by = "-2 weeks")[-1])
  W <- array(NA, dim = c(T + 1, D + 1, length(ddChangepoint)), dimnames = list(as.character(t02s),
                                                                               paste("delay", 0:D, sep = ""),
                                                                               c(as.character(ddChangepoint))))
  for (t in 0:T) {
    for (i in 1:length(ddChangepoint)) {
      W[t + 1, ,i] <- pmax(0, as.numeric((t02s[t+1] + 0:D) - ddChangepoint[i]))
    }
  }
  
  # Wtrue matrix
  ## IID Log Ga cps in delay dist
  ddChangepoint_true <- c(lubridate::ymd("2020-03-01", 
                                         "2020-03-15",
                                         "2020-03-25",
                                         "2020-04-15", 
                                         "2020-04-30"))
  ddChangepoint_true = ddChangepoint_true[which(ddChangepoint_true<now)] 
  W_true<- array(NA, dim = c(T + 1, D + 1, length(ddChangepoint_true)), 
                 dimnames = list(as.character(t02s),
                                 paste("delay", 0:D, sep = ""),
                                 c(as.character(ddChangepoint_true))))
  for (t in 0:T) {
    for (i in 1:length(ddChangepoint_true)) {
      W_true[t + 1, ,i] <- pmax(0, as.numeric((t02s[t+1] + 0:D) - ddChangepoint_true[i]))
    }
  }

  ## Weekday W
  ## Additional part of the discrete survival model design matrix W = [ W_cp W_extra]
  # Make designmatrix for weekdays
  # Make a factor for each day of the week
  wdays <- levels(lubridate::wday(t02s, label=TRUE))[-1]
  # Create the extension of the W matrix
  Wextra <- array(NA, dim=c(length(t02s), D+1, length(wdays)),
                  dimnames=list(as.character(t02s), paste("delay",0:D,sep=""), wdays))

  # Replace easter holidays and may 1st with sunday indicator
  adjust_day = function(nc_dates_t, D) {
    if_else((nc_dates_t + 0:D) %in% c(ymd("2020-04-10"), 
                                      ymd("2020-04-13"), 
                                      ymd("2020-05-01"), 
                                      ymd("2020-06-01")), 
            ymd("2020-04-12"), nc_dates_t + 0:D)
  }

  # Loop over all times and lags
  for (t in seq_len(length(t02s))) {
    for (w in seq_len(length(wdays))) {
      #
      Wextra[t,, w] <- as.numeric(lubridate::wday(adjust_day(t02s[t], D), label=TRUE) == wdays[w])
    }
  }

  W_wd = array(NA, dim = c(T + 1, D + 1, length(ddChangepoint) + length(wdays)),
               dimnames = list(as.character(t02s),
                               paste("delay", 0:D, sep = ""),
                               c(as.character(ddChangepoint),
                                 as.character(wdays))))
  for (t in 0:T) {
    for (i in 1:length(ddChangepoint)) {
      W_wd[t + 1, ,i] = pmax(0, as.numeric((t02s[t+1] + 0:D) - ddChangepoint[i]))
    }
    W_wd[t+1, , (length(ddChangepoint)+1):(length(ddChangepoint) + length(wdays))] = Wextra[t+1,,]
  }

  ## W_wd 1-week cps
  ddChangepoint_1w = sort(seq(dMax, dMin, by = "-1 weeks")[-1])
  W_1w_wd = array(NA, dim = c(T + 1, D + 1, length(ddChangepoint_1w) + length(wdays)),
                  dimnames = list(as.character(t02s),
                                  paste("delay", 0:D, sep = ""),
                                  c(as.character(ddChangepoint_1w),
                                    as.character(wdays))))
  for (t in 0:T) {
    for (i in 1:length(ddChangepoint_1w)) {
      W_1w_wd[t + 1, ,i] = pmax(0, as.numeric((t02s[t+1] + 0:D) - ddChangepoint_1w[i]))
    }
    W_1w_wd[t+1, , (length(ddChangepoint_1w)+1):(length(ddChangepoint_1w) + length(wdays))] = Wextra[t+1,,]
  }
  
  # rW delay_dist
  rW_params = ceiling((T+1)/2) - 1
  rW_delay_des = matrix(nrow=(T+1), ncol = rW_params)
  for (i in (T+1):2) {
    if ((T+1) %% 2 == 0) {
    rW_delay_des[i,ceiling((i)/2)-1] = 1
    } else {
      rW_delay_des[i,ceiling((i+1)/2)-1] = 1
    }
  }
  rW_delay_des[is.na(rW_delay_des)] = 0
  
  ## True data
  nT_true=n_full[((T+1)-D+1):((T+1)-predLag),]
  # Combine all data/preprocessd objects into list
  if (mod_name == "negBinom_rW_rW") {
    dat_list = list(T = T+1,
                    maxDelay = D,
                    rT = n,
                    mu_gamma = gamma_prior_kappa$mu_gamma,
                    sd_gamma = gamma_prior_kappa$sigma_gamma,
                    predLag = predLag)
  } else if (mod_name == "poisson_iid_consDel") {
    dat_list = list(T = T+1, 
                    maxDelay = D,
                    rT = n, 
                    alpha_lambda = alpha.lambda,
                    beta_lambda = beta.lambda,
                    mu_gamma = gamma_prior_kappa$mu_gamma,
                    sd_gamma = gamma_prior_kappa$sigma_gamma,
                    predLag = predLag)
  } else if (mod_name == "negBinom_rW_cp2W_wd") {
    dat_list = list(T = T+1, 
                    maxDelay = D,
                    n_cp=length(ddChangepoint),
                    n_wextra = length(wdays),
                    rT = n, 
                    W = W_wd,
                    mu_gamma = gamma_prior_kappa$mu_gamma,
                    sd_gamma = gamma_prior_kappa$sigma_gamma,
                    eta_mu = rep(0, length(ddChangepoint) + length(wdays)),
                    eta_sd = c(rep(0.01, length(ddChangepoint)),
                                  rep(0.5, length(wdays))),
                    predLag = predLag)
  } else if (mod_name == "negBinom_rW_rW_2d") {
    dat_list = list(T = T+1,
                    maxDelay = D,
                    rT = n,
                    mu_gamma = gamma_prior_kappa$mu_gamma,
                    sd_gamma = gamma_prior_kappa$sigma_gamma,
                    rW_params = rW_params,
                    rW_delay_des = rW_delay_des,
                    predLag = predLag)
  } else if (mod_name == "negBinom_rW_rW_wd") {
    dat_list = list(T = T+1,
                    maxDelay = D,
                    rT = n,
                    mu_gamma = gamma_prior_kappa$mu_gamma,
                    sd_gamma = gamma_prior_kappa$sigma_gamma,
                    n_wextra = length(wdays),
                    W = Wextra,
                    eta_mu = rep(0, length(wdays)),
                    eta_sd = rep(0.5, length(wdays)),
                    predLag = predLag)
  } else if (mod_name == "negBinom_rW_cp2W") {
    dat_list = list(T = T+1, 
                    maxDelay = D,
                    n_cp=length(ddChangepoint),
                    rT = n, 
                    W = W,
                    mu_gamma = gamma_prior_kappa$mu_gamma,
                    sd_gamma = gamma_prior_kappa$sigma_gamma,
                    eta_mu = as.array(rep(0, length(ddChangepoint))),
                    eta_sd = as.array(rep(0.01, length(ddChangepoint))),
                    predLag = predLag)
  } else if (mod_name == "poisson_rW_cp2W") {
    dat_list = list(T = T+1, 
                    maxDelay = D,
                    n_cp=length(ddChangepoint),
                    rT = n, 
                    W = W,
                    mu_gamma = gamma_prior_kappa$mu_gamma,
                    sd_gamma = gamma_prior_kappa$sigma_gamma,
                    eta_mu = as.array(rep(0, length(ddChangepoint))),
                    eta_sd = as.array(rep(0.01, length(ddChangepoint))),
                    predLag = predLag)
  } else if (mod_name == "poisson_rW_const") {
    dat_list = list(T = T+1, 
                    maxDelay = D,
                    rT = n, 
                    mu_gamma = gamma_prior_kappa$mu_gamma,
                    sd_gamma = gamma_prior_kappa$sigma_gamma,
                    predLag = predLag)
  } else if (mod_name == "negBinom_rW_cptrue") {
    dat_list = list(T = T+1, 
                    maxDelay = D,
                    n_cp=length(ddChangepoint_true),
                    rT = n, 
                    W = W_true,
                    mu_gamma = gamma_prior_kappa$mu_gamma,
                    sd_gamma = gamma_prior_kappa$sigma_gamma,
                    eta_mu = as.array(rep(0, length(ddChangepoint_true))),
                    eta_sd = as.array(rep(0.01, length(ddChangepoint_true))),
                    predLag = predLag)
  } else if (mod_name == "poisson_rW_cptrue") {
    dat_list = list(T = T+1, 
                    maxDelay = D,
                    n_cp=length(ddChangepoint_true),
                    rT = n, 
                    W = W_true,
                    mu_gamma = gamma_prior_kappa$mu_gamma,
                    sd_gamma = gamma_prior_kappa$sigma_gamma,
                    eta_mu = as.array(rep(0, length(ddChangepoint_true))),
                    eta_sd = as.array(rep(0.01, length(ddChangepoint_true))),
                    predLag = predLag)
  } else {
    dat_list = NULL
  }
  dat_list = list_modify(dat_list, t02s = t02s, nT_true = nT_true)
  dat_list
}

# Function to summarise Nowcast results, estimate Rt and compute scoring rules
summarise_stan_fit = function(fit, dataDays, predLag, D, now, nT_true, data) {
  ntInf_est = rstan::extract(fit, "ntInf")$ntInf
  nc_dates = seq(now-D+1, now-predLag, by = 1)
  # Obtain posterior summary for the estimated number of disease onsets and
  # corresponding scoring rules
  ntInf = tibble(date = nc_dates,
                 med = apply(ntInf_est, 2, median),
                 q025 = apply(ntInf_est, 2, function(x) quantile(x, 0.025)),
                 q05 = apply(ntInf_est, 2, function(x) quantile(x, 0.05)),
                 q25 = apply(ntInf_est, 2, function(x) quantile(x, 0.25)),
                 q75 = apply(ntInf_est, 2, function(x) quantile(x, 0.75)),
                 q95 = apply(ntInf_est, 2, function(x) quantile(x, 0.95)),
                 q975 = apply(ntInf_est, 2, function(x) quantile(x, 0.975)),
                 log_score=logs_sample(rowSums(nT_true), t(ntInf_est)),
                 crps_score = crps_sample(rowSums(nT_true), t(ntInf_est)),
                 type = fit@model_name)
  
  p_est = rstan::extract(fit, "p")$p
  p = apply(p_est, c(2,3), function(x) median(x))
  colnames(p) = paste0("D", 0:(ncol(p)-1))
  p = tibble(date = seq(min(dataDays), max(dataDays), by = 1),
             as_tibble(p)) %>% filter(date <= max(dataDays)-predLag)
  # Estimate R(t)
  gt_ni = c(mean = 4.70, sd = 2.90)
  gt = R0::generation.time("lognormal", gt_ni)
  set.seed(131015)
  now_r = now - predLag
  time_2ndcase = ymd("2020-02-24")
  
  most_secondary_transmissions_occurred <- gt$time[min(which(cumsum(gt$GT)
                                                             >= 0.95))]
  # End of R(t) estimation
  end_Rt_estim <- now_r - most_secondary_transmissions_occurred + 1
  
  # Estimate R0 using the Wallinga & Teunis method.
  # sample time-series from posterior of now-cast object
  # each row is one smaple from posterior
  n_sample = 300
  mcmc_samples = sample(1:nrow(ntInf_est), size = n_sample, replace = FALSE)
  ntInf_obs = data %>% group_by(date = disease_start) %>% summarise(n_dis = n()) %>%
    right_join(tibble(date=seq(time_2ndcase, min(nc_dates)-1, by = "1 day")))
  
  posterior_samples = pmax(cbind(do.call(rbind, lapply(1:n_sample, function(x) unlist(ntInf_obs$n_dis))),
                                 ntInf_est[mcmc_samples,]),0)
  posterior_samples[is.na(posterior_samples)] = 0
  library(future)
  plan("multicore")
  options(mc.cores = 4)
  Rt_est_list <- furrr::future_map(
    .x = seq_len(nrow(posterior_samples)),
    .f = ~{
      ts <- posterior_samples[.x,]
      est.R0.TD( # use custom function that returns R.simu
        epid    = ts,
        GT      = gt,
        t       = seq(time_2ndcase, now-predLag, by = "1 day"),
        correct = TRUE,
        begin   = time_2ndcase,
        end     = end_Rt_estim,
        nsim    = 500)
    }
  )
  
  Rt_df = do.call(rbind, lapply(
    Rt_est_list,
    function(x) {
      as_tibble(x$R.simu) %>% mutate(Date = seq(time_2ndcase, now-predLag, 1)) %>%
        filter(Date<=end_Rt_estim) %>% pivot_longer(cols = starts_with("V"), values_to = "Rt")
    }))
  
  Rt_df_smry <- Rt_df %>%
    group_by(Date) %>%
    summarize(
      Rt_lower = quantile(Rt, .025),
      Rt_upper = quantile(Rt, .975),
      Rt = mean(Rt))
  # Save posterior information, estimated delay distribution and estimated Rt
  list(ntInf = ntInf, p=p, Rt_df_smry=Rt_df_smry)
}

#' Perform nowcast for a specific 'now' based on specific model
#' @param dat_list list of prepared data objects
#' @now current date for nowcast
#' @param begin_date date to start estimation
#' @para D maximum considered delay
#' @param predLag difference between `now` and the last desired prediction of 
#' number of disease onsets
#' @model stan model object
#' 
nowcast_stan_models = function(dat_list,
                               now,
                               begin_date = ymd("2020-02-24"),
                               D = 21,
                               predLag = 2,
                               model=NULL,
                               model_name_res=NULL,
                               data) {


  # Estimate model
  nc_mod = sampling(model,
                    data = dat_list,
                    iter = 4000,
                    chains = 4,
                    include = TRUE,
                    pars = c("ntInf", "p"),
                    cores = 4)
  res = tryCatch(summarise_stan_fit(nc_mod,
                                    dataDays = dat_list$t02s,
                                    predLag = predLag, now = now, D=D,
                                    nT_true=dat_list$nT_true, 
                                    data = data),
                 error = function(e) list(ntInf = tibble(date=NA, med=NA,
                                                         q025=NA, q05=NA, q25 = NA,
                                                         q75 = NA, q95=NA, q975 = NA,
                                                         log_score = NA,
                                                         crps_score = NA,
                                                         type = NA),
                                          p = tibble(date=NA, as_tibble(matrix(rep(NA, D+1), nrow = 1,
                                                                               dimnames = list(NULL, paste0("D", 0:D))))),
                                          Rt_df_smry = tibble(Date=NA, Rt_lower=NA, Rt_upper=NA, Rt=NA)))
  
  ntInf = res$ntInf %>% mutate(model = model_name_res)
  p = res$p %>% mutate(model = model_name_res)
  Rt = res$Rt_df_smry %>% mutate(model = model_name_res)
  
  list(ntInf = ntInf %>% mutate(now=now), p = p, Rt=Rt %>% mutate(now=now))
}


estimate_nowcast = function(now,
                            data,
                            begin_date = ymd("2020-02-24"),
                            D = 21,
                            predLag = 2,
                            model,
                            mod_name) {
  # dat_mod = data %>% filter(rep_date<= now)
  dat_list = prepare_data_list(dat=data,
                               now = now,
                               begin_date = begin_date,
                               D = D,
                               predLag = predLag, mod_name = mod_name)
  res = nowcast_stan_models(dat_list = dat_list,
                            now = now, 
                            begin_date = begin_date, D = D, predLag = predLag, 
                            model = model,
                            model_name_res = mod_name,
                            data = data)
  res
}

