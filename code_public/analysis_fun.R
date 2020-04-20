#' Helper functions

#' Format results of nowcast call (stsNC object) to tibble
#' @param nc stsNC object containing nowcast results
#' @param nowcastDates the dates for which a nowcast was performed (as supplied to the original nowcast call)
#' @return tibble containing basic description of nowcast call and per fay of nowcastDates the date, observed number
#' of newly infected cases, posterior median of nowcast call (predicted number of newly infected cases) and
#' associated 95% prediction interval
tidy_nc = function(nc, nowcastDates) {
  nc %>% tidy() %>% as_tibble() %>%
    mutate(pi_lower = nc@pi[,,1],  pi_upper=  nc@pi[,,2],
           type = paste0(nc@control$call$method, "_D=",nc@control$call$D,"_m=",nc@control$call$m)) %>%
    filter(epoch %in% nowcastDates) %>%
    dplyr::select(type, date, observed, predicted = upperbound, predicted_lower= pi_lower, predicted_upper = pi_upper)
}

#' Plot estimated distribution of lambda_t
#' Plot the estimated number of newly infected cases per day of a nowcast model and show the available data
#' of newly infected cases per day used in estimation
#' @param nc stsNC object containing nowcast results
#' @param safePredictLag difference between `now` and the last desired prediction of new cases
plot_lambda_t = function(nc, safePredictLag) {plot(nc, legend.opts=NULL,
                              xaxis.tickFreq=list("%d"=atChange,"%m"=atChange),
                              xaxis.labelFreq=list("%d"=at2ndChange),
                              xaxis.labelFormat="%d-%b",
                              xlab="Time (days)",
                              lty=c(1,1,1,1),
                              lwd=c(1,1,2),
                              las=2,
                              cex.axis=0.7,
                              ylab="No. symptom onsets",
                              ylim=c(0,max(observed(nc),
                                           upperbound(nc),
                                           predint(nc),
                                           na.rm=TRUE)))
  if (nc@control$method == "bayes.trunc.ddcp") {
    lambda <- attr(nc@delayCDF[["bayes.trunc.ddcp"]],"model")$lambda
    showIdx <- seq(ncol(lambda) - safePredictLag)
    matlines(
      showIdx,t(lambda)[showIdx,], col="gray", lwd=c(1,2,1), lty=c(2,1,2))
  }
}


#' Estimate different nowcast models
#' This function takes imputed data of reporting date and dates of disease onset and
#' estimates nocasts based on three different models:
#' 1) Lawless with assumed maximum delay of D=21 estimated based on the most current m=21 days.
#' 2) Fully bayesian approach with D=21 and tps prior for estimating lambda_t. The delay distribution is
#' modeled based on discrete breakpoints every two weeks (going backwards from `now`) and dummy effects of
#' weekdays of the reporting date.
#' 3) As 2) but with iidLogGa prior for estimating lambda_t
#' `now` is set to the maximum of the reporting date (rep_date column in imputed_data) - 1 day.
#' @param imputed_data imputed data based on function perform_imputation
#' @param data_date the date when the data was received
#' @param safePredictLag difference between `now` and the last desired prediction of newly infected cases, default to 1
#' @param save_results logical; should the results and plots of the nowcast be safed in a (new) folder at filepath
#' `../results/nowcast_estimates/data_date,add_folder_path/`, `../results/nowcast_plots/data_date,add_folder_path/`
#' @param add_folder_path see save_results
#' @return list including list of estimated nowcasts, nowcast results as tibble used for plotting, two plots of the nowcast results
#' and list with additional information on the nowcast (`now`, `safePredictLag`, `nowcastDates`, `so_range`)
estimate_nowcasts = function(imputed_data,
                             data_date,
                             safePredictLag = 2,
                             save_results = FALSE,
                             add_folder_path = "") {

  ## Prepare folder to save results
  if (save_results & !dir.exists(paste0("../results/nowcast_estimates/",data_date, add_folder_path,"/"))) {
    dir.create(paste0("../results/nowcast_estimates/",data_date, add_folder_path,"/"))
  }

  if (save_results & !dir.exists(paste0("../results/nowcast_plots/",data_date,add_folder_path,"/"))) {
    dir.create(paste0("../results/nowcast_plots/",data_date,add_folder_path,"/"))
  }

  # Select relevant columns for nowcast
  dat_mod = imputed_data %>%
    dplyr::select(rep_date, rep_date_weekday, disease_start_imp)
  # Now is day before most current due to incomplete reporting on actual day
  now <- data_date - 1

  dat_mod = dat_mod %>% filter(rep_date <= now)
  so_range <- range(ymd("2020-02-24"), now)
  nowcastDates <- seq(from = so_range[1], to = now - safePredictLag, by = 1)

  sts <- linelist2sts(
    as.data.frame(dat_mod),
    dateCol="disease_start_imp",
    aggregate.by = "1 day",
    dRange = so_range)

  nc_control = list(
    N.tInf.max = 4e3,
    N.tInf.prior = structure("poisgamma",
                             mean.lambda = mean(observed(sts)),
                             var.lambda = 5*var(observed(sts))),
    gd.prior.kappa = .1,
    predPMF = TRUE,
    dRange = so_range)

  print("Perform lawless nowcast")
  nc_lawless = surveillance::nowcast(now = now,
                                     when = nowcastDates, data = as.data.frame(dat_mod),
                                     dEventCol = "disease_start_imp",
                                     dReportCol = "rep_date",
                                     aggregate.by = "1 day",
                                     D = 21,
                                     method = "lawless",
                                     m = 21,
                                     control = nc_control)
  if (save_results) {
    png(paste0("../results/nowcast_plots/", data_date, add_folder_path,"/nc_lawless_delay.png"), width = 800, height = 350)
    plot(nc_lawless, type="delay", dates=seq(so_range[1], now, by="1 day"), w=1)
    dev.off()
    png(paste0("../results/nowcast_plots/", data_date, add_folder_path, "/nc_lawless_lambda.png"), width = 800, height = 350)
    plot_lambda_t(nc_lawless, safePredictLag = safePredictLag)
    dev.off()
  }

  # Bayesian model iidLogGa prior
  # Two breakpoints in delay distribution every last 2 weeks
  cp_vec = sort(seq(now, ymd("2020-03-01"), by = "-2 week")[-1])
  source("./general/nowcast_w.R")
  source("./general/plot_stsNC.R")
  # Bayesian model weekday effects
  print("Perform bayesian tps + weekday nowcast")
  ## Additional part of the discrete survival model design matrix W = [ W_cp W_extra]
  D <- 21
  nc_dates <- seq(nc_control$dRange[1], now,  by="1 day")
  # Make designmatrix for weekdays
  # Make a factor for each day of the week
  wdays <- levels(lubridate::wday(nc_dates, label=TRUE))[-1]
  # Create the extension of the W matrix
  Wextra <- array(NA, dim=c(length(nc_dates), length(wdays), D+1),  
                  dimnames=list(as.character(nc_dates), wdays, paste("delay",0:D,sep="")))
  
  # Replace easter holidays with sunday indicator
  adjust_day = function(nc_dates_t, D) {
    if_else((nc_dates[t] + 0:D) %in% c(ymd("2020-04-10"), ymd("2020-04-13")), ymd("2020-04-12"), nc_dates[t] + 0:D) 
  }
  
  # Loop over all times and lags
  for (t in seq_len(length(nc_dates))) {
    for (w in seq_len(length(wdays))) {
      #
      Wextra[t,w, ] <- as.numeric( lubridate::wday(adjust_day(nc_dates[t], D), label=TRUE) == wdays[w])
    }
  }
  # TPS
  nc_bayes_tps_wd <- nowcast(now = now,
                             when = nowcastDates,
                             data = as.data.frame(dat_mod),
                             dEventCol = "disease_start_imp",
                             dReportCol = "rep_date",
                             aggregate.by = "1 day",
                             D = 21,
                             ##  Use the discrete time survival model with change-points
                             method = "bayes.trunc.ddcp",
                             control = modifyList(nc_control, list(
                               ddcp = list(
                                 ddChangepoint = cp_vec,
                                 Wextra = Wextra,
                                 logLambda = "tps",
                                 eta.mu=rep(0,length(cp_vec) + if (is.null(Wextra)) 0 else ncol(Wextra)),
                                 eta.prec=diag(rep(1,length(cp_vec) + if (is.null(Wextra)) 0 else ncol(Wextra))),
                                 mcmc=c(burnin=5000,sample=20000,thin=1, adapt=3000)
                               )
                             ))
  )
  if (save_results) {
    png(paste0("../results/nowcast_plots/", data_date, add_folder_path, "/nc_bayes_tps_wd_delay.png"), width = 800, height = 350)
    stsNC_plotDelay(nc_bayes_tps_wd,  dates=seq(so_range[1], now, by="1 day"), w=1)
    dev.off()
    png(paste0("../results/nowcast_plots/", data_date, add_folder_path, "/nc_bayes_tps_wd_lambda.png"), width = 800, height = 350)
    plot_lambda_t(nc_bayes_tps_wd, safePredictLag = safePredictLag)
    dev.off()
  }
  # iidLogGa
  print("Perform bayesian iidLogGa + weekday nowcast")
  nc_bayes_iidlogga_wd <- nowcast(now = now, when = nowcastDates,
                             data = as.data.frame(dat_mod),
                             dEventCol = "disease_start_imp",
                             dReportCol = "rep_date",
                             aggregate.by = "1 day",
                             D = 21,
                             ##  Use the discrete time survival model with change-points
                             method = "bayes.trunc.ddcp",
                             control = modifyList(nc_control, list(
                               ddcp = list(
                                 ddChangepoint = cp_vec,
                                 Wextra = Wextra,
                                 logLambda = "iidLogGa",
                                 eta.mu=rep(0,length(cp_vec) + if (is.null(Wextra)) 0 else ncol(Wextra)),
                                 eta.prec=diag(rep(1,length(cp_vec) + if (is.null(Wextra)) 0 else ncol(Wextra))),
                                 mcmc=c(burnin=5000,sample=20000,thin=1, adapt=3000)
                               )
                             ))
  )
  if (save_results) {
    png(paste0("../results/nowcast_plots/", data_date, add_folder_path, "/nc_bayes_iidlogga_wd_delay.png"), width = 800, height = 350)
    stsNC_plotDelay(nc_bayes_iidlogga_wd,  dates=seq(so_range[1], now, by="1 day"), w=1)
    dev.off()
    png(paste0("../results/nowcast_plots/", data_date, add_folder_path, "/nc_bayes_iidlogga_wd_lambda.png"), width = 800, height = 350)
    plot_lambda_t(nc_bayes_iidlogga_wd, safePredictLag)
    dev.off()
  }


  ## Combine nowcast results and plot
  nc_lawless_td = tidy_nc(nc_lawless, nowcastDates = nowcastDates) %>%
    mutate(type = "lawless")
  nc_bayes_tps_wd_td = tidy_nc(nc_bayes_tps_wd, nowcastDates = nowcastDates) %>%
    mutate(type ="bayes_wd_tps")
  nc_bayes_iidlogga_wd_td = tidy_nc(nc_bayes_iidlogga_wd, nowcastDates = nowcastDates) %>%
    mutate(type ="bayes_wd_iidLogGa")

  reported = imputed_data %>% group_by(rep_date) %>% summarise(sum_rep = n()) %>%
    filter(rep_date %in% nc_lawless_td$date)
  act_observed = imputed_data %>% group_by(disease_start) %>% summarise(sum_rep = n()) %>%
    filter(disease_start %in% nc_lawless_td$date)

  plot_dat = rbind(nc_lawless_td, nc_bayes_tps_wd_td, nc_bayes_iidlogga_wd_td)

  plot_dat_pred = plot_dat %>% dplyr::select(type, date, predicted, predicted_lower, predicted_upper) %>%
    pivot_longer(cols = c("predicted", "predicted_upper", "predicted_lower"), names_to = "pred_type")
  plot_dat = rbind(reported %>% rename(date=rep_date, value = sum_rep) %>%
                     mutate(type = "reported_new_cases", pred_type = "predicted") %>%
                     dplyr::select(type, date, pred_type, value),
                   act_observed %>% rename(date = disease_start, value = sum_rep) %>%
                     mutate(type = "obs_onset", pred_type = "predicted") %>%
                     dplyr::select(type, date, pred_type, value),
                   plot_dat %>% dplyr::select(date, observed) %>% group_by(date) %>% summarise(value = mean(observed)) %>%
                     mutate(type = "obs_imp_onset", pred_type = "predicted") %>%
                     dplyr::select(type, date, pred_type, value),
                   plot_dat_pred)


  colours = c("reported_new_cases" = "red",
              "obs_imp_onset" = "black",
              "lawless" = "#A6611A",
              "bayes_wd_iidLogGa" = "#DFC27D",
              "bayes_wd_tps" = "#018571",
              "obs_onset" = "lightblue")
  plot_nc = ggplot(aes(date, value, col = type, lty = pred_type),
                   data = plot_dat %>% filter(pred_type %in% c("predicted"),
                                              type %in% c("reported_new_cases",
                                                          "obs_imp_onset",
                                                          "bayes_wd_iidLogGa",
                                                          "bayes_wd_tps"))) +
    geom_line() +
    scale_color_manual(values = c(colours), 
                       labels = c("reported_new_cases",
                                  "obs_imp_onset",
                                  "bayes_wd_iidLogGa" = "onset_bayes_wd_iidLogGa",
                                  "bayes_wd_tps" = "onset_bayes_wd_tps")) +
    scale_linetype_manual(values = c("predicted" = 1,
                                     "predicted_lower"=2,
                                     "predicted_upper"=2)) +
    guides(linetype=FALSE) +
    scale_x_date(date_breaks = "1 week") +
    xlab("Date") +
    ylab("Cases") +
    coord_cartesian(xlim = c(ymd("2020-02-22"), now)) +
    theme_light() +
    guides(linetype=FALSE) +
    geom_point(aes(x=now, y = 0), color = "gold", pch = 14,
               size = 1.5, show.legend = FALSE)

  plot_nc_ci = ggplot(aes(date, value, col = type, lty = pred_type),
                      data = plot_dat %>%
                        filter(type %in% c("reported_new_cases",
                                           "obs_imp_onset",
                                           "bayes_wd_iidLogGa",
                                           "bayes_wd_tps"))) +
    geom_line() +
    scale_color_manual(values = c(colours), 
                                    labels = c("reported_new_cases",
                                               "obs_imp_onset",
                                               "bayes_wd_iidLogGa" = "onset_bayes_wd_iidLogGa",
                                               "bayes_wd_tps" = "onset_bayes_wd_tps")) +
    guides(linetype=FALSE) +
    scale_linetype_manual(values = c("predicted" = 1,
                                     "predicted_lower"=2,
                                     "predicted_upper"=2)) +
    scale_x_date(date_breaks = "1 week") +
    xlab("Date") +
    ylab("Cases") +
    coord_cartesian(xlim = c(ymd("2020-02-22"), now)) +
    theme_light() +
    guides(linetype=FALSE) +
    geom_point(aes(x=now, y = 0), color = "gold", pch = 14,
               size = 1.5, show.legend = FALSE)

  line_colours = c("bayes_wd_tps" = "#018571",
                   "obs_onset" = "lightgrey",
                   "obs_imp_onset" = "lightgrey")
  fill_colours = c("bayes_wd_tps" = "#018571",
                   "obs_onset" = "lightblue",
                   "obs_imp_onset" = "grey")
  fill_labels = c("bayes_wd_tps" = "Nowcast",
                  "obs_onset" = "Erkrankungsbeginn Erhoben",
                  "obs_imp_onset" = "Erkrankungsbeginn Imputiert")
  plot_dat2 = pivot_wider(plot_dat, id_cols = c(date, type), names_from = pred_type)
  
  plot_tps = ggplot(aes(date, predicted, col = type, fill = type),
                    data = plot_dat2 %>%
                      filter(type %in% c("bayes_wd_tps"))) +
    scale_color_manual(values = line_colours) +
    scale_fill_manual(values = fill_colours,
                      labels = fill_labels) +
    scale_x_date(date_breaks = "1 week") +
    xlab("Datum") +
    ylab("Fälle") +
    coord_cartesian(xlim = c(ymd("2020-02-22"), now)) +
    geom_col(aes(date, predicted), data = plot_dat2 %>% filter(type == "obs_imp_onset")) +
    geom_col(aes(date, predicted), data = plot_dat2 %>% filter(type == "obs_onset")) +
    geom_line(size=.8) + 
    geom_ribbon(aes(ymin=predicted_lower, ymax = predicted_upper), alpha = .5, linetype = 1, lwd = 0.4) + 
    theme_light() +
    theme(legend.title = element_blank()) +
    guides(linetype=FALSE, color = FALSE) +
    geom_point(aes(x=now, y = 0), color = "gold", pch = 14,
               size = 1.5, show.legend = FALSE)


  if (save_results) {
    ggsave(plot_nc,
           file = paste0("../results/nowcast_plots/", data_date, add_folder_path,
                         "/Nowcast_Covid19_Bavaria_",now+1,".png"),
           width = 8, height = 3, dpi = 500)
    ggsave(plot_nc_ci,
           file = paste0("../results/nowcast_plots/", data_date, add_folder_path,
                         "/Nowcast_Covid19_Bavaria_InclPI_", now+1,".png"),
           width = 8, height = 3, dpi = 500)
    ggsave(plot_tps,
           file = paste0("../results/nowcast_plots/", data_date, add_folder_path,
                         "/Nowcast_Covid19_Bavaria_obs_imp", now+1,".png"),
           width = 8, height = 3, dpi = 500)

  }

  res_list = list(nc_obj = list(nc_lawless = nc_lawless,
                                nc_bayes_tps_wd = nc_bayes_tps_wd,
                                nc_bayes_iidlogga_wd = nc_bayes_iidlogga_wd),
                  plot_dat = plot_dat,
                  plot_nc = plot_nc,
                  plot_nc_ci = plot_nc_ci,
                  plot_tps = plot_tps,
       add_inf = list(now = now, safePredictLag = safePredictLag, nowcastDates = nowcastDates, so_range = so_range)
  )
  if (save_results) {
    save(res_list, file = paste0("../results/nowcast_estimates/", data_date, add_folder_path, "/nc_res_list.RData"))
  }
  res_list
}


#' Perform impuation of disease onset
#' This function takes edited individual-specific data with information on disease reporting date (column name 'rep_date'),
#' reporting week ('rep_week'), weekday of reporting date ('rep_date_weekday'), an individuals age ('age') and the disease onset
#' date (i.e. symptom onset date) ('disease_start') and delay (time between symptom onset and reporting date) if available and 
#' estimates a Weibull GAMLSS for the expected delay. There exist two different options for the linear predictor of the GAMLSS model.
#' Afterwards missing disease onset days are imputed based on the fitted GAMLSS model.
#' @param dat_edit dataframe with columns rep_date, rep_date_weekay, age, disease_start and delay (NA if not available)
#' @param type string specifiying the GAMLSS imputation model, either  "week_weekday_age" for smooth effects of the reporting week and individiuals age
#' and a dummy effect for the reporting weekday, or "week" for only a smooth effect of the reporting week
#' @return list of imputation results, first entry: imputed data; second entry: GAMLSS model for imputation; third entry: ggplot illustrating
#' GAMLSS results; fourth entry: training data used for fitting the GAMLSS.
#' 
perform_imputation = function(dat_edit, type = "week_weekday_age") {
  # select training data with available delay
  # Select training data
  train <- dat_edit %>%
  dplyr::select(rep_date, rep_week, disease_start, delay, rep_date_weekday, age) %>%
    filter(delay >= 0 ) %>%
    mutate(delay_prep_log = if_else(delay == 0, 1e-2, delay),
           event_observed = 1)
  train_complete <- train %>% dplyr::select(delay_prep_log, rep_week, rep_date, rep_date_weekday, age) %>%
    na.omit()
  ## Estimate imputation model
  if (type == "week_weekday_age") {
    impute_model <- gamlss::gamlss(
      delay_prep_log  ~ gamlss::cs(rep_week) + gamlss::cs(age) + rep_date_weekday,
      sigma.formula = ~ gamlss::cs(rep_week) + gamlss::cs(age) + rep_date_weekday,
      family = gamlss.dist::WEI2,
      data = train_complete,
      method = mixed(80,80)
    )} else if(type == "week") {
      impute_model <- gamlss::gamlss(
        delay_prep_log  ~ gamlss::cs(rep_week), 
        sigma.formula = ~ gamlss::cs(rep_week),
        family = gamlss.dist::WEI2,
        data = train_complete,
        method = mixed(80,80)
      )} else {
        stop("No valid imputation model specified")
      }
  ## Plot imputation results
  if (type == "week_weekday_age") {
    plot_imp_mod_dat <- data.frame(expand.grid(rep_week = unique(train$rep_week), age=seq(0,100,5),
                                                rep_date_weekday = unique(train$rep_date_weekday)))
    plot_imp = plot_imp_mod_dat %>% mutate(mu =exp(predict(impute_model, newdata = plot_imp_mod_dat, 
                                                           what = "mu",
                                                           data = train_complete)),
                                           sigma = exp(predict(impute_model, newdata = plot_imp_mod_dat, 
                                                               what="sigma",
                                                               data = train_complete)),
                                           med = gamlss.dist:::qWEI2(0.5, mu=mu, sigma=sigma),
                                           q25 = gamlss.dist:::qWEI2(0.25, mu=mu, sigma=sigma),
                                           q75 = gamlss.dist:::qWEI2(0.75, mu=mu, sigma=sigma),
                                           rep_week = factor(rep_week, levels = sort(unique(rep_week)), 
                                                             labels = sort(unique(rep_week))),
                                           rep_date_weekday = factor(rep_date_weekday,
                                                                     levels = c("Monday", "Tuesday", "Wednesday",
                                                                                "Thursday", "Friday", "Saturday",
                                                                                "Sunday")))
    
    imp_median_weibull = ggplot(aes(age, med, col = rep_week), data = plot_imp) +
      geom_line() + facet_grid(.~rep_date_weekday) +
      scale_color_discrete(name="Week") +
      theme_bw() + theme(legend.position = "bottom") +
      guides(color = guide_legend(nrow = 1)) +
      xlab(label = "Age") + ylab("Delay") +
      coord_cartesian(xlim=c(20,90))
  } else {
    plot_imp_mod_dat <- data.frame(expand.grid(rep_week = unique(train$rep_week)))
    plot_imp = plot_imp_mod_dat %>% mutate(mu =exp(predict(impute_model, newdata = plot_imp_mod_dat, 
                                                           what = "mu", data = train_complete)),
                                           sigma = exp(predict(impute_model, newdata = plot_imp_mod_dat, 
                                                               what="sigma", data = train_complete)),
                                           med = gamlss.dist:::qWEI2(0.5, mu=mu, sigma=sigma))
    imp_median_weibull = NULL
  }
  ## Impute dat_edit
  if (type == "week_weekday_age") {
    dat_imputed <- dat_edit %>%
      mutate(
        # Get estimated person-specific location and scale
        mu = exp(predict(impute_model, newdata = dat_edit %>% dplyr::select(rep_week, age, rep_date_weekday),
                         what = "mu", data = train_complete)),
        sigma = exp(predict(impute_model, newdata = dat_edit %>% dplyr::select(rep_week, age, rep_date_weekday), 
                            what = "sigma", data = train_complete)),
        # Randomly sample from person-specific delay distribution
        delay_imputed = gamlss.dist::rWEI2(nrow(dat_edit), mu = mu, sigma = sigma),
        # Select imputed or reported delay
        delay_mod = ifelse(is.na(delay), delay_imputed, delay),
        # Get day of symptom onset assuming all cases reported at noon (i.e. delay <0.5 -> reporting at same 
        # day as symptom onset)
        disease_start_imp = ymd(rep_date - delay_mod + 0.5))
  } else {
    dat_imputed <- dat_edit %>%
      mutate(
        mu = exp(predict(impute_model, newdata = dat_edit %>% dplyr::select(rep_week), what = "mu", 
                         data = train_complete)),
        sigma = exp(predict(impute_model, newdata = dat_edit %>% dplyr::select(rep_week), what = "sigma", 
                            data = train_complete)),
        delay_imputed = gamlss.dist::rWEI2(nrow(dat_edit), mu = mu, sigma = sigma),
        delay_mod = ifelse(is.na(delay), delay_imputed, delay),
        disease_start_imp = ymd(rep_date - delay_mod + 0.5))
  }

  list(dat_imputed = dat_imputed, impute_model = impute_model, impute_plot = imp_median_weibull, train_dat = train_complete)
}

## Function to summarise data
summarize_data = function(dat_edit) {
  delay_summary <- dat_edit %>%
    group_by(rep_week) %>%
    summarise(n = n(),
              n_delay_available = sum(!is.na(delay)),
              prop_delay_available = scales::percent(mean(!is.na(delay))),
              mean = mean(delay, na.rm = TRUE),
              median = median(delay, na.rm = TRUE),
              q25 = quantile(delay, prob=0.25, na.rm = TRUE),
              q75 = quantile(delay, prob=0.75, na.rm = TRUE))

 delay_summary
}

