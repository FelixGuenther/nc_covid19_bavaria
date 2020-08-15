#' Helper functions

#' Format results of STAN nowcast to (estimated ntInf and delay distribution)
#' @param fit fitted STAN model
#' @param dataDays dates of the nowcast days for formatting
#' @return list contraining tibbles of summaries of estimated ntInf and posterior median of delay distribution per day
summarise_stan_fit = function(fit, dataDays, predLag, dat_mod, D, now) {
  ntInf_est = rstan::extract(fit, "ntInf")$ntInf
  ntInf_est = tibble(date = seq(now-D+1, now - predLag, by = 1),
                 med = apply(ntInf_est, 2, median),
                 q025 = apply(ntInf_est, 2, function(x) quantile(x, 0.025)),
                 q05 = apply(ntInf_est, 2, function(x) quantile(x, 0.05)),
                 q25 = apply(ntInf_est, 2, function(x) quantile(x, 0.25)),
                 q75 = apply(ntInf_est, 2, function(x) quantile(x, 0.75)),
                 q95 = apply(ntInf_est, 2, function(x) quantile(x, 0.95)),
                 q975 = apply(ntInf_est, 2, function(x) quantile(x, 0.975)),
                 type = fit@model_name)
  observed = dat_mod %>% group_by(date = disease_start) %>% summarise(observed = n()) %>%
    filter(date >= min(dataDays))
  observed = left_join(tibble(date=seq(ymd("2020-02-24"), now-predLag, by = 1)), 
                       observed, by = "date")
  ntInf = right_join(ntInf_est, observed, by = "date")
  ntInf = ntInf %>% mutate(med = ifelse(date<=now-D, observed, med),
                           q025 = ifelse(date<=now-D, observed, q025),
                           q05 = ifelse(date<=now-D, observed, q05),
                           q25 = ifelse(date<=now-D, observed, q25),
                           q75 = ifelse(date<=now-D, observed, q75),
                           q95 = ifelse(date<=now-D, observed, q95),
                           q975 = ifelse(date<=now-D, observed, q975))
  p_est = rstan::extract(fit, "p")$p
  p = apply(p_est, c(2,3), function(x) median(x))
  colnames(p) = paste0("D", 0:(ncol(p)-1))
  p = tibble(date = seq(min(dataDays), max(dataDays), by = 1),
             as_tibble(p)) %>% filter(date <= max(dataDays)-predLag)
  list(ntInf = ntInf, p=p)
}


#' Function to prepare data for nowcasting in STAN
prepare_data_list = function(dat,
                             now,
                             begin_date = ymd("2020-02-24"),
                             D = 21,
                             predLag = 2) {
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
  T <- length(t02s) - 1
  
  n <- matrix(NA, nrow = T + 1, ncol = T + 1, dimnames = list(as.character(t02s),
                                                              NULL))
  timeDelay <- function(d1, d2) {
    as.numeric(d2 - d1)
  }
  dat = dat %>% mutate(delay = rep_date - disease_start)
  
  for (t in 0:T) {
    data.att <- dat[which(dat$disease_start == t02s[t+1]), ]
    for (x in 0:(T - t)) {
      n[t + 1, x + 1] <- sum(data.att$delay == x)
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
  ddChangepoint <- sort(seq(dMax, dMin, by = "-2 weeks")[-1])
  W <- array(NA, dim = c(T + 1, D + 1, length(ddChangepoint)), dimnames = list(as.character(t02s),
                                                                               paste("delay", 0:D, sep = ""),
                                                                               c(as.character(ddChangepoint))))
  for (t in 0:T) {
    for (i in 1:length(ddChangepoint)) {
      W[t + 1, ,i] <- pmax(0, as.numeric((t02s[t+1] + 0:D) - ddChangepoint[i]))
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
    if_else((nc_dates_t + 0:D) %in% c(ymd("2020-04-10"), ymd("2020-04-13"), ymd("2020-05-01"), 
                                      ymd("2020-06-01")), ymd("2020-04-12"), nc_dates_t + 0:D)
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
  # Combin all data/preprocessd objects into list
  list(T = T+1,
       maxDelay = D,
       n_cp = length(ddChangepoint),
       n_cp_1w = length(ddChangepoint_1w),
       n_wextra = length(wdays),
       rT = n,
       W = W,
       W_wd=W_wd,
       W_1w_wd = W_1w_wd,
       alpha_lambda = alpha.lambda,
       beta_lambda = beta.lambda,
       mu_gamma = gamma_prior_kappa$mu_gamma,
       sd_gamma = gamma_prior_kappa$sigma_gamma,
       eta_mu = rep(0, length(ddChangepoint)),
       eta_sd = c(rep(0.01, length(ddChangepoint))),
       eta_mu_wd = rep(0, length(ddChangepoint) + length(wdays)),
       eta_sd_wd = c(rep(0.01, length(ddChangepoint)),
                     rep(0.5, length(wdays))),
       eta_mu_1w_wd = rep(0, length(ddChangepoint_1w) + length(wdays)),
       eta_sd_1w_wd = c(rep(0.01, length(ddChangepoint_1w)),
                        rep(0.5, length(wdays))),
       predLag = predLag,
       t02s = t02s)
}

#' Estimate nowcast models
#' This function takes imputed data of reporting date and dates of disease onset and
#' estimates nocasts based STAN model
#' @param imputed_data imputed data based on function perform_imputation
#' @param data_date the date when the data was received
#' @param safePredictLag difference between `now` and the last desired prediction of newly infected cases, default to 1
#' @param save_results logical; should the results and plots of the nowcast be safed in a (new) folder under path
#' `../results/nowcast_estimates/data_date,add_folder_path/`, `../results/nowcast_plots/data_date,add_folder_path/`
#' @param add_folder_path see save_results
#' @return list including list of estimated nowcasts, nowcast results as tibble used for plotting, two plots of the nowcast results
#' and list with additional information on the nowcast (`now`, `safePredictLag`)
estimate_nowcasts = function(imputed_data,
                             data_date,
                             safePredictLag = 2,
                             save_results = FALSE,
                             add_folder_path = "",
                             D=21) {
  mod_randomWalk_cpDelay_wExtra_negBinom = stan_model("general/randomWalk_cpDelay_wExtra_negBinom.stan")
  ## Prepare folder to save results
  if (save_results & !dir.exists(paste0("../results/nowcast_estimates/",data_date, add_folder_path,"/"))) {
    dir.create(paste0("../results/nowcast_estimates/",data_date, add_folder_path,"/"))
  }

  # Select relevant columns for nowcast
  dat_mod = imputed_data %>%
    dplyr::select(rep_date = rep_date_reg, 
                  rep_date_weekday = rep_date_reg_weekday, 
                  disease_start = disease_start_imp)
  # Now is day before most current due to incomplete reporting on actual day
  now <- data_date - 1

  dat_mod = dat_mod %>% filter(rep_date <= now)
  prep_dat_list = prepare_data_list(dat = dat_mod,
                                    now = now,
                                    begin_date = ymd("2020-02-24"), 
                                    D = D, 
                                    predLag = safePredictLag)
  
  randomWalk_cpDelay_WD_dat = prep_dat_list[c("T", "maxDelay", "n_cp", "n_wextra",
                                         "rT", "W_wd", "mu_gamma", "sd_gamma",
                                         "eta_mu_wd", "eta_sd_wd", "predLag")]
  names(randomWalk_cpDelay_WD_dat)[c(6,9,10)] = c("W", "eta_mu", "eta_sd")
  
  nc_randomWalk_cpDelay_WD = sampling(mod_randomWalk_cpDelay_wExtra_negBinom,
                                       data = randomWalk_cpDelay_WD_dat,
                                       iter = 2500,
                                       chains = 4,
                                       include = TRUE,
                                       pars = c("ntInf", "p"),
                                       cores = 4)
  nc_summary = summarise_stan_fit(nc_randomWalk_cpDelay_WD, 
                                  dataDays = seq(ymd("2020-02-24"), now, by = 1), 
                                  predLag = safePredictLag, 
                                  D = D, 
                                  dat_mod = dat_mod, 
                                  now = now)
  
  
  
  

  res_list = list(nc_obj = nc_randomWalk_cpDelay_WD,
                  nc_smry = nc_summary,
       add_inf = list(now = now, safePredictLag = safePredictLag, D=D)
  )
  if (save_results) {
    save(res_list, file = paste0("../results/nowcast_estimates/", data_date, add_folder_path, "/nc_res_list.RData"))
  }
  res_list
}



#' Perform impuation of disease onset
#' This function takes edited individual-specific data with information on disease reporting date at the local and regional
#' health authorities (column name 'rep_date_local', 'rep_date_reg'),
#' reporting week ('rep_week_local','rep_week_reg'), weekday of reporting date ('rep_date_local_weekday', 'rep_date_reg_weekday'), 
#' an individuals age ('age') and the disease onset date (i.e. symptom onset date) ('disease_start') and delay 
#' (time between symptom onset and local reporting date) if available and 
#' estimates a Weibull GAMLSS for the expected delay. There exist two different options for the linear predictor of the GAMLSS model.
#' Afterwards missing disease onset days are imputed based on the fitted GAMLSS model.
#' @param dat_edit dataframe with columns rep_date_local, rep_date_reg, rep_date_local_weekday, rep_date_reg_weekday, 
#' age, disease_start and delay (NA if not available)
#' @param type string specifiying the GAMLSS imputation model, either  "week_weekday_age" for smooth effects of the reporting week and individiuals age
#' and a dummy effect for the reporting weekday, or "week" for only a smooth effect of the reporting week
#' @return list of imputation results, first entry: imputed data; second entry: GAMLSS model for imputation; third entry: ggplot illustrating
#' GAMLSS results; fourth entry: training data used for fitting the GAMLSS.

perform_imputation = function(dat_edit, type = "week_weekday_age") {
  # select training data with available delay
  # Select training data
  train <- dat_edit %>%
    dplyr::select(rep_date_local, rep_week_local, disease_start, delay, rep_date_local_weekday, age) %>%
    filter(delay >= 0 ) %>%
    mutate(delay_prep_log = if_else(delay == 0, 1e-2, delay),
           event_observed = 1)
  train_complete <- train %>% dplyr::select(delay_prep_log, rep_week_local, rep_date_local, 
                                            rep_date_local_weekday, age) %>%
    na.omit()
  ## Estimate imputation model
  if (type == "week_weekday_age") {
    impute_model <- gamlss::gamlss(
      delay_prep_log  ~ gamlss::cs(rep_week_local) + gamlss::cs(age) + rep_date_local_weekday,
      sigma.formula = ~ gamlss::cs(rep_week_local) + gamlss::cs(age) + rep_date_local_weekday,
      family = gamlss.dist::WEI2,
      data = train_complete,
      method = mixed(80,80)
    )} else if(type == "week") {
      impute_model <- gamlss::gamlss(
        delay_prep_log  ~ gamlss::cs(rep_week_local), 
        sigma.formula = ~ gamlss::cs(rep_week_local),
        family = gamlss.dist::WEI2,
        data = train_complete,
        method = mixed(80,80)
      )} else {
        stop("No valid imputation model specified")
      }
  ## Plot imputation results
  if (type == "week_weekday_age") {
    plot_imp_mod_dat <- data.frame(expand.grid(rep_week_local = unique(train$rep_week_local), age=seq(0,100,5),
                                               rep_date_local_weekday = unique(train$rep_date_local_weekday)))
    plot_imp_mod_dat = plot_imp_mod_dat %>% mutate(mu =exp(predict(impute_model, newdata = plot_imp_mod_dat, 
                                                what = "mu",
                                                data = train_complete)),
                                sigma = exp(predict(impute_model, newdata = plot_imp_mod_dat, 
                                                    what="sigma",
                                                    data = train_complete)))
    plot_imp = plot_imp_mod_dat %>% mutate(med = gamlss.dist:::qWEI2(0.5, mu=mu, sigma=sigma),
                                           q25 = gamlss.dist:::qWEI2(0.25, mu=mu, sigma=sigma),
                                           q75 = gamlss.dist:::qWEI2(0.75, mu=mu, sigma=sigma),
                                           rep_week_local = factor(rep_week_local, levels = sort(unique(rep_week_local)), 
                                                             labels = sort(unique(rep_week_local))),
                                           rep_date_local_weekday = factor(rep_date_local_weekday,
                                                                     levels = c("Monday", "Tuesday", "Wednesday",
                                                                                "Thursday", "Friday", "Saturday",
                                                                                "Sunday")))
    
    imp_median_weibull = ggplot(aes(age, med, col = rep_week_local), data = plot_imp) +
      geom_line() + facet_grid(.~rep_date_local_weekday) +
      scale_color_discrete(name="Week") +
      theme_bw() + theme(legend.position = "bottom") +
      guides(color = guide_legend(nrow = 1)) +
      xlab(label = "Age") + ylab("Delay") +
      coord_cartesian(xlim=c(20,90))
  } else {
    imp_median_weibull = NULL
  }
  ## Impute dat_edit
  if (type == "week_weekday_age") {
    dat_imputed <- dat_edit %>%
      mutate(
        # Get estimated person-specific location and scale
        mu = exp(predict(impute_model, newdata = dat_edit %>% dplyr::select(rep_week_local, age, 
                                                                            rep_date_local_weekday),
                         what = "mu", data = train_complete)),
        sigma = exp(predict(impute_model, newdata = dat_edit %>% dplyr::select(rep_week_local, age, 
                                                                               rep_date_local_weekday), 
                            what = "sigma", data = train_complete)),
        # Randomly sample from person-specific delay distribution
        delay_imputed = gamlss.dist::rWEI2(nrow(dat_edit), mu = mu, sigma = sigma),
        # Select imputed or reported delay
        delay_mod = ifelse(is.na(delay), delay_imputed, delay),
        # Get day of symptom onset assuming all cases reported at noon (i.e. delay <0.5 -> reporting at same 
        # day as symptom onset)
        disease_start_imp = ymd(rep_date_local - delay_mod + 0.5))
  } else {
    dat_imputed <- dat_edit %>%
      mutate(
        mu = exp(predict(impute_model, newdata = dat_edit %>% dplyr::select(rep_week_local), what = "mu", 
                         data = train_complete)),
        sigma = exp(predict(impute_model, newdata = dat_edit %>% dplyr::select(rep_week_local), what = "sigma", 
                            data = train_complete)),
        delay_imputed = gamlss.dist::rWEI2(nrow(dat_edit), mu = mu, sigma = sigma),
        delay_mod = ifelse(is.na(delay), delay_imputed, delay),
        disease_start_imp = ymd(rep_date_local - delay_mod + 0.5))
  }
  
  list(dat_imputed = dat_imputed, impute_model = impute_model, impute_plot = imp_median_weibull, train_dat = train_complete)
}


# Function to summarize data
summarize_data = function(dat_edit) {
  delay_summary <- dat_edit %>%
    group_by(rep_week_local) %>%
    summarise(n = n(),
              n_delay_available = sum(!is.na(delay)),
              prop_delay_available = scales::percent(mean(!is.na(delay))),
              mean = mean(delay, na.rm = TRUE),
              median = median(delay, na.rm = TRUE),
              q25 = quantile(delay, prob=0.25, na.rm = TRUE),
              q75 = quantile(delay, prob=0.75, na.rm = TRUE))
  delay_summary
}

#' Estimate R(t)
estimate_Rt = function(nc_res, data_date) {
  gt_ni = c(mean = 4.70, sd = 2.90)
  gt = R0::generation.time("lognormal", gt_ni)
  set.seed(131015)
  source("general/est.R0.TD.R") # customized version of R0::est.R0.TD
  source("general/functions-R0.R")
  now = nc_res$add_inf$now - nc_res$add_inf$safePredictLag
  time_2ndcase = nc_res$nc_smry$ntInf$date[which.min(cumsum(nc_res$nc_smry$ntInf$observed)>=2)]
  
  most_secondary_transmissions_occurred <- gt$time[min(which(cumsum(gt$GT)
                                                             >= 0.95))]
  # End of R(t) estimation
  end_Rt_estim <- now - most_secondary_transmissions_occurred
  
  # Estimate R0 using the Wallinga & Teunis method.
  # sample time-series from posterior of now-cast object
  # each row is one smaple from posterior
  n_sample = 300
  ntInf_mcmc = rstan::extract(nc_res$nc_obj, "ntInf")$ntInf
  mcmc_samples = sample(1:nrow(ntInf_mcmc), size = n_sample, replace = FALSE)
  posterior_samples = pmax(cbind(do.call(rbind, lapply(1:n_sample, function(x) nc_res$nc_smry$ntInf$med[1:(length(nc_res$nc_smry$ntInf$med)-nc_res$add_inf$D+nc_res$add_inf$safePredictLag)])),
                            ntInf_mcmc[mcmc_samples,]),0)
  
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
        t       = nc_res$nc_smry$ntInf$date,
        correct = TRUE,
        begin   = time_2ndcase,
        end     = end_Rt_estim,
        nsim    = 500)
    }
  )
  
  Rt_df = do.call(rbind, lapply(
    Rt_est_list,
    function(x) {
      as_tibble(x$R.simu) %>% mutate(Date = seq(min(nc_res$nc_smry$ntInf$date), now, 1)) %>%
        filter(Date<=end_Rt_estim) %>% pivot_longer(cols = starts_with("V"), values_to = "Rt")
    }))
  
  Rt_df_smry <- Rt_df %>%
    group_by(Date) %>%
    summarize(
      Rt_lower = quantile(Rt, .025),
      Rt_upper = quantile(Rt, .975),
      Rt = mean(Rt))
  Rt_df_smry
}  

#' create_plots

create_plots = function(est_ntInf, Rt, imputed_data, now, safePredictLag) {
  plot_dat = left_join(est_ntInf, imputed_data %>% group_by(date=disease_start_imp) %>% summarise(obs_imp_onset=n()))
  plot_dat = left_join(plot_dat, imputed_data %>% group_by(date=disease_start) %>% summarise(obs_onset=n()))
  plot_dat = left_join(plot_dat, imputed_data %>% group_by(date=rep_date_reg) %>% summarise(reported_new_cases=n())) %>%
    mutate(obs_onset=replace_na(obs_onset, 0),
           reported_new_cases = replace_na(reported_new_cases, 0))
  
  
  fill_colours = c("med" = "#018571",
                   "obs_onset" = "lightblue",
                   "obs_imp_onset" = "grey")
  fill_labels = c("med" = "Nowcast",
                  "obs_onset" = "Erkrankungsbeginn erhoben",
                  "obs_imp_onset" = "Erkrankungsbeginn imputiert")
  
  plot_nc = ggplot(plot_dat %>% filter(date<=now - safePredictLag)) +
    scale_x_date(breaks = seq(max(plot_dat$date), min(plot_dat$date), by = "-1 week"), date_labels = "%d.%m.") +  
    xlab("Datum") +
    ylab("Anzahl Fälle") +
    # geom_vline(aes(xintercept=ymd("2020-03-22")), lty = 2, col = "lightgrey") +
    coord_cartesian(xlim = c(ymd("2020-02-24"), now)) +
    geom_col(aes(date, obs_imp_onset, fill = "obs_imp_onset")) +
    geom_col(aes(date, obs_onset, fill = "obs_onset")) +
    geom_line(aes(date, med), col = fill_colours["med"]) + 
    geom_ribbon(aes(date, ymin = q025, ymax = q975, fill = "med", col = "med"), 
                alpha = .25, lwd = .25) +
    scale_color_manual(values = fill_colours,
                       labels = fill_labels) +
    scale_fill_manual(values = fill_colours,
                      labels = fill_labels) +
    
    theme_bw() +
    theme(legend.title = element_blank(), legend.position = "bottom") +
    geom_point(aes(x=now, y = 0), color = "gold", pch = 14,
               size = 1.5, show.legend = FALSE) +
    guides(linetype=FALSE, color = FALSE) +
    labs(caption = "Statistisches Beratungslabor StaBLab, LMU München;  Department of Mathematics, Stockholm University \nDaten: Bayerisches Landesamt für Gesundheit und Lebensmittelsicherheit LGL")
  
  
  colorvec = c("Krankheitsbeginn\nerhoben+imputiert" = "black",
               "Krankheitsbeginn\nNowcast" = "#018571",
               "Registrierte Fälle" = "red")
  plot_nc_rep = ggplot(plot_dat %>% filter(date<=now - safePredictLag)) +
    scale_x_date(breaks = seq(max(plot_dat$date), min(plot_dat$date), by = "-1 week"), date_labels = "%d.%m.") +  
    xlab("Datum") +
    ylab("Anzahl Fälle") +
    # geom_vline(aes(xintercept=ymd("2020-03-22")), lty = 2, col = "lightgrey") +
    coord_cartesian(xlim = c(ymd("2020-02-24"), now)) +
    geom_line(aes(date, reported_new_cases, col = "Registrierte Fälle")) +
    geom_line(aes(date, obs_imp_onset, col = "Krankheitsbeginn\nerhoben+imputiert")) +
    geom_line(aes(date, med, col = "Krankheitsbeginn\nNowcast")) +
    scale_color_manual(values = colorvec) +
    geom_ribbon(aes(date, ymin = q025, ymax = q975), fill = "#018571", col = "#018571", lwd = .25, alpha = .25) +
    theme_bw() +
    geom_point(aes(x=now, y = 0), color = "gold", pch = 14,
               size = 1.5, show.legend = FALSE) +
    theme(legend.title = element_blank(), legend.position = "bottom") +
    labs(caption = "Statistisches Beratungslabor StaBLab, LMU München;  Department of Mathematics, Stockholm University \nDaten: Bayerisches Landesamt für Gesundheit und Lebensmittelsicherheit LGL")
  
  plot_Rt = ggplot(Rt, aes(x = Date, y = Rt)) +
    geom_ribbon(aes(ymin = Rt_lower, ymax = Rt_upper), fill = "lightgray", alpha = .9) +
    geom_line() +
    xlab("Zeit") + ylab(expression(R(t))) +
    geom_hline(yintercept=1, lty=2) +
    scale_y_continuous(breaks=seq(0,6, by=1)) +
    scale_x_date(breaks = seq(max(Rt$Date), min(Rt$Date), by = "-1 week"), 
                 labels = paste0(strftime(seq(max(Rt$Date), min(Rt$Date), by = "-1 week"), 
                                          format = "%d.%m."), 
                                 " -\n", 
                                 strftime(seq(max(Rt$Date), min(Rt$Date), by = "-1 week")+10, 
                                          format = "%d.%m."),
                                 "  ")) +  
    coord_cartesian(ylim = c(0,5)) +
    theme_bw() +
    labs(caption = "Statistisches Beratungslabor StaBLab, LMU München; Department of Mathematics, Stockholm University \nDaten: Bayerisches Landesamt für Gesundheit und Lebensmittelsicherheit LGL")
  gg_list = list(plot_nc, plot_nc_rep, plot_Rt)
  gg_list
}
