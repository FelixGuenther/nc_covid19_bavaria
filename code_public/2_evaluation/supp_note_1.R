#' ---
#' title: "Nowcasting the COVID-19 pandemic in Bavaria"
#' subtitle: "Supplemental Note: Evaluation of the nowcasting approach"
#' author: "Felix Günther, Andreas Bender, Katharina Katz, Helmut Küchenhoff, Michael Höhle"
#' date: ""
#' ---

#' ## Introduction
#' We perform an evaluation of the Bayesian hierarchical nowcasting based on synthetic data mimicking the reported 
#' Bavarian COVID--19 data and retrospectively on the official data from the LGL that was reported until July, 31. 
#' 
#' ## Synthetic data
#' 
#' #### Data generating process
#' 
#' For simulation of the synthetic data, we utilized a smoothed version of the observed number of reported disease onsets per day 
#' and specified a reporting delay model similar to the model described in equation 3 of the manuscript. In this 
#' discrete time hazard model for the reporting delay, we specified a linear time effect on the log-hazard 
#' with five change points over time.
#' The aggregated number of newly diseased cases (green, solid) and number of reported cases per day (red, dotted) 
#' are:

#+ message=FALSE, warning=FALSE, echo=FALSE, fig.height = 3.5, fig.width = 9, fig.cap="Daily numbers of newly diseased cases and reported cases in synthetic data" 
library(xtable)
library(knitr)
library(lubridate)
library(tidyverse)
load("../../data_public/evaluation/dat_mod_synthetic.RData")
dat_mod_synth = dat_mod
rm(dat_mod)
ggplot() + geom_line(aes(date, n_dis), dat_mod_synth %>% group_by(date=disease_start) %>%
                           summarise(n_dis=n()) %>%
                          filter(date <= ymd("2020-07-01")), col = "green") +
              geom_line(aes(date, n_rep), dat_mod_synth %>% group_by(date=rep_date) %>%
                          summarise(n_rep=n()) %>%
                          filter(date <= ymd("2020-07-01")), col = "red", lty = 2) +
              theme_bw() + ylab("No. cases")

#' The changepoints in the linear time effect of the delay distribution lead to (smooth) changes
#' in the delay distribution. We did not add any further effects of e.g., weekdays to the delay distribution model, and
#' simulated reporting dates for each case with disease onset on day $t$ directly from the delay distribution without
#' adding any further variability.
#'  
#' The utilized delay distribution can be illustrated by plotting the empirical median and 25%- and 75%-quantile
#' of the sampled reportig delay over time:

#+ message=FALSE, warning=FALSE, echo=FALSE, fig.height = 3.5, fig.width = 9, fig.cap="Empirical reporting delay of the synthetic data over time"
dat_mod_synth %>% mutate(delay=rep_date-disease_start) %>%
  group_by(date = disease_start) %>%
  summarise(med= as.numeric(median(delay)),
            q25 = as.numeric(quantile(delay, 0.25)),
            q75 = as.numeric(quantile(delay, 0.75))) %>%
  ggplot() + geom_line(aes(date, med), col = "green") +
  geom_ribbon(aes(date, ymin = q25, ymax = q75), alpha=.25, fill="green") +
  coord_cartesian(xlim=c(ymd("2020-02-28"), ymd("2020-06-30"))) + 
  theme_bw() + ylab("Reporting delay")

#' ### Results
#' 
#' We esimated nowcasts for all dates $t$ from March, 17 (22 days after disease onset of first case) until June, 30
#' by restricting the data to all cases reported until the respective date and compared the nowcast predicitions for all days
#' $t-6,\ldots, t-2$ to the actual true numbers of newly diseased cases per day.
#' 
#' Nowcasts are performed based on six different models (see description in the manuscript), and the performance of the models
#' is compared via proper scoring rules, the root mean squared error and coverage frequencies of 95% prediction intervals.
#' 
#' In the following table, we present results based on the nowcasts for all dates and restricted to the time period between
#' March, 17 and April, 30, where most of the dynamic in the epidemic curve happened:
 
#+ message=FALSE, warning=FALSE, echo=FALSE, fig.height = 3.5, fig.width = 9
synth_files = c(
  "../../results_public/2_evaluation/synthetic_data/negBinom_rW_cp_2w.RData",
  "../../results_public/2_evaluation/synthetic_data/negBinom_rW_cp_true.RData",
  "../../results_public/2_evaluation/synthetic_data/negBinom_rW_rW.RData",
  "../../results_public/2_evaluation/synthetic_data/poiss_rW_const.RData",
  "../../results_public/2_evaluation/synthetic_data/poiss_rW_cp_2w.RData",
  "../../results_public/2_evaluation/synthetic_data/poiss_rW_cp_true.RData"
)

ntInf_synth = do.call(rbind, lapply(synth_files,
function(x) {
  load(file = x)
  do.call(rbind, lapply(res_list, function(x) x$ntInf))
}))

# Mutate data
ntInf_synth = ntInf_synth %>% mutate(lag = as.numeric(now-date),
                                     model = factor(model, levels = c("poisson_rW_const",
                                                                      "poisson_rW_cp2W",
                                                                      "negBinom_rW_cp2W",
                                                                      "negBinom_rW_rW",
                                                                      "poisson_rW_cptrue",
                                                                      "negBinom_rW_cptrue")))
# Compute coverage of 95% PIs for all nowcast predictions 2-6 days before 'now'
cov_full = left_join(ntInf_synth, dat_mod_synth %>% 
                       group_by(date = disease_start) %>% 
                       summarise(n_dis = n())) %>%
  mutate(coverage = q025<=n_dis & n_dis <= q975) %>%
  filter(lag %in% 2:6) %>%
  group_by(model) %>% summarise(frac_cov = mean(coverage))

cov_may = left_join(ntInf_synth, dat_mod_synth %>% group_by(date = disease_start) %>% summarise(n_dis = n())) %>%
  mutate(coverage = q025<=n_dis & n_dis <= q975) %>%
  filter(lag %in% 2:6, date <ymd("2020-05-01")) %>%
  group_by(model) %>% summarise(frac_cov = mean(coverage))

# Summarise Scoring Rules
crps_full = ntInf_synth %>% filter(lag %in% c(2:6)) %>% group_by(model) %>% summarise(crps = mean(crps_score))
crps_may = ntInf_synth %>% filter(lag %in% c(2:6),
                            date < ymd("2020-05-01")) %>% group_by(model) %>% summarise(crps = mean(crps_score))
logsc_full = ntInf_synth %>% filter(lag %in% c(2:6)) %>% group_by(model) %>% summarise(logS = mean(log_score))
logsc_may = ntInf_synth %>% filter(lag %in% c(2:6),
                             date < ymd("2020-05-01")) %>% group_by(model) %>% summarise(logS = mean(log_score))
rmse_full = ntInf_synth %>% filter(lag %in% c(2:6)) %>%
  left_join(dat_mod_synth %>% group_by(date = disease_start) %>%
              summarise(truth = n())) %>%
  mutate(diff = med - truth) %>%
  group_by(model) %>% summarise(rmse = sqrt(mean(diff^2)))
rmse_may = ntInf_synth %>% filter(lag %in% c(2:6),
                            date < ymd("2020-05-01")) %>%
  left_join(dat_mod_synth %>% group_by(date = disease_start) %>%
              summarise(truth = n())) %>%
  mutate(diff = med - truth) %>%
  group_by(model) %>% summarise(rmse = sqrt(mean(diff^2)))

smry_synth_paper = tibble(#model = crps_full$model,
  DistrNtd = c("Poisson", "Poisson", rep("Neg. Binomial", 2), "Poisson", "Neg. Binomial"),
  DelayMod = c("No changes", "Lin. effect of time Changepoints 2 weeks", "Lin. effect of time Changepoints 2 weeks",
               "Daily changes (first order RW)",
               "Lin. effect of time true changepoints", "Lin. effect of time true changepoints"),
  CRPS = crps_full$crps,
  logS = logsc_full$logS,
  RMSE = rmse_full$rmse,
  CovNtInf = cov_full$frac_cov,
  CovRt = rep("-", 6))


kable(tibble(model = crps_full$model, crps_f = crps_full$crps, crps_m = crps_may$crps,
             logS_f = logsc_full$logS, logS_m = logsc_may$logS,
             rmse_f = rmse_full$rmse, rmse_m = rmse_may$rmse,
             cov_f = cov_full$frac_cov, cov_m = cov_may$frac_cov), digits = 2, caption="Quantification of the performance of six different nowcast models on synthetic data. Shown are the average metrics over all nowcast dates (*_f) and restricted to the period until May, 1 (*_m). Reported scores are the continuous ranked probability score, logarithmic score, root mean squared error of posterior median, and coverage frequencies of 95% prediction intervals. Estimated models are Poisson and Negative Binomial models with 1) (assumed) constant delay distribution, 2) linear time-effects with changepoints every two weeks before now (*_cp2W), 3) daily changes in delay distribution based on a first-order random walk prior (*_rW) and 4) true changepoints of the data generating process (*_cptrue).")

#' In the following plot we show nowcast predictions and 95%-PIs 3 days before the current date over time for each model, and compare them
#' with the true number of newly diseased cases:

#+ message=FALSE, warning=FALSE, echo=FALSE, fig.height = 8, fig.width = 6, fig.cap="Nowcast predictions and 95% PIs for nowcasts 3 days before 'now' compared to true number of disease onsets (solid) and number of reported cases (dotted)."
ntInf_synth %>%
  dplyr::filter(lag==3) %>% ggplot() + geom_line(aes(date, med, col = model)) +
  geom_ribbon(aes(date, ymin=q025, ymax=q975, fill = model), alpha = .3) +
  geom_line(aes(date, n_dis), dat_mod_synth %>% group_by(date=disease_start) %>%
              summarise(n_dis=n()) %>%
              filter(date <= ymd("2020-07-01"))) +
  geom_line(aes(date, n_rep), dat_mod_synth %>% group_by(date=rep_date) %>%
              summarise(n_rep=n()) %>%
              filter(date <= ymd("2020-07-01")), lty = 2) +
  theme_bw() +
  facet_grid(model~.) +
  coord_cartesian(xlim = c(ymd("2020-02-24"), ymd("2020-05-01")), ylim = c(0,1200)) +
  theme(legend.title = element_blank(),legend.position = "bottom") + ylab("Number cases")

#' \pagebreak
#' 


#' Note that in the quantitative evaluation based on scores and coverage frequencies, we focus on nowcasts 2-6 days before now and not only
#' 3 days before.
#' 
#' We can also compare the estimated delay distribution (for the most current date, at each nowcast day $t$ the delay
#' distribution is estimated for the complete past) with the empirical delay distribution from the sampled data:

#+ message=FALSE, warning=FALSE, echo = FALSE, fig.cap="Comparison of estimated delay distribution (black) and empirical delay distribution (blue) per day for each model. Shown are the median and 25% and 75% quantiles."
p_synth = do.call(rbind, lapply(synth_files,
                          function(x) {
                            load(x)
                            do.call(rbind, lapply(res_list, function(x) tryCatch(x[[2]] %>% filter(date==max(date)-1),
                                                                                 error = function(e) rep(NA, 24))))
                          }))
p_synth = p_synth %>% filter(!is.na(date)) %>% arrange(date)
p_synth_long = p_synth %>% pivot_longer(cols = starts_with("D", ignore.case = F))
delay_est_synth = p_synth_long  %>% group_by(date, model) %>%
  summarise(med_delay = min(which(cumsum(value)>=.5))-1,
            q10_delay = min(which(cumsum(value)>=.10))-1,
            q25_delay = min(which(cumsum(value)>=.25))-1,
            q40_delay = min(which(cumsum(value)>=.40))-1,
            q60_delay = min(which(cumsum(value)>=.60))-1,
            q75_delay = min(which(cumsum(value)>=.75))-1,
            q90_delay = min(which(cumsum(value)>=.90))-1
  )
delay_est_synth = delay_est_synth %>% mutate(model = factor(model,
                                                levels = c("poisson_rW_const",
                                                           "poisson_rW_cp2W",
                                                           "negBinom_rW_cp2W",
                                                           "negBinom_rW_rW",
                                                           "poisson_rW_cptrue",
                                                           "negBinom_rW_cptrue")))
delay_true_synth = dat_mod_synth %>% mutate(delay = pmin(21, as.numeric(rep_date - disease_start))) %>%
  group_by(date = disease_start) %>% summarise(med_delay = median(delay),
                                               q10_delay = quantile(delay, .1),
                                               q25_delay = quantile(delay, .25),
                                               q40_delay = quantile(delay, .4),
                                               q60_delay = quantile(delay, .6),
                                               q75_delay = quantile(delay, .75),
                                               q90_delay = quantile(delay, .90)) %>%
  filter(date %in% delay_est_synth$date)
ggplot() +
  geom_ribbon(aes(date, ymin=q25_delay, ymax=q75_delay), delay_est_synth, alpha = .1) +
  geom_line(aes(date, med_delay), delay_est_synth) +
  facet_wrap(~model, ncol = 2) + theme_bw() +
  geom_ribbon(aes(date, ymin=q25_delay, ymax=q75_delay), delay_true_synth, alpha = .1, fill = "green") +
  geom_line(aes(date, med_delay), delay_true_synth, col = "green", lwd = 1.5)
 
#' \pagebreak

#'We summarize the results the following way: when we supply the true, in reality unknown, changepoints of the delay distribution to model fitting, the nowcast performs 
#'best with respect to our evaluation metrics. 
#'It shows the lowest log and CRPS score, lowest root mean squared error and shows the desired coverage frequencies 
#'for the 95%-prediction intervals. With the models assuming changepoints in the linear time effect on reporting delay every two weeks before $T$, 
#'we obtain similar, but slightly worse performance. The approach appears to be able to capture the moderate changes in the delay distribution successfully. 
#'Modeling the changes on a daily basis shows a slightly worse performance with respect to the CRPS score and PI coverage frequencies. The prediction intervals are
#'wider and there exists some evidence for convergence problems at the beginning of the nowcast period where to median delay was strongly overestimated leading
#'to an upward bias in the predicted number of newly diseased cases per day. This might be related to the overly complex model for the reporting delay
#'and very few observations to estimate the reporting delay up to this time. Assuming a constant reporting delay distribution over time and ignoring the changes leads to the worst performance with biggest scores 
#'and low coverage frequencies of the prediction intervals. The number of estimated newly diseased cases is overestimated during the whole nowcasting
#'period starting at end of March. When specifying an adequate model for the delay distribution, the distributional assumptions regarding $N_{t,d}$ play a minor role for
#'performance on the synthetic data.
#'
#' \pagebreak
#' 
#' ## Retrospective evaluation based on Bavarian data

#+ message=FALSE, warning=FALSE, echo=FALSE, fig.height = 3.5, fig.width = 9
bav_files = c("../../results_public/2_evaluation/bavaria/poiss_rW_const.RData",
              "../../results_public/2_evaluation/bavaria/poiss_rW_cp_2w.RData",
              "../../results_public/2_evaluation/bavaria/negBinom_rW_rW.RData",
              "../../results_public/2_evaluation/bavaria/negBinom_rW_rW_wd.RData",
              "../../results_public/2_evaluation/bavaria/negBinom_rW_cp_2w.RData",
              "../../results_public/2_evaluation/bavaria/negBinom_rW_cp_2w_wd.RData")
# True data summary
dat_bav_smry = read_tsv("../../data_public/evaluation/dat_bavaria_truth_epcurve_reported.csv")

ntInf = do.call(rbind, lapply(bav_files,
                              function(x) {
  load(file = x)
  do.call(rbind, lapply(res_list, function(x) x$ntInf))
  }))
ntInf = ntInf %>% mutate(lag = as.numeric(now-date),
                         model = factor(model, levels = c("poisson_rW_const",
                                                          "poisson_rW_cp2W",
                                                          "negBinom_rW_cp2W",
                                                          "negBinom_rW_rW",
                                                          "negBinom_rW_cp2W_wd",
                                                          "negBinom_rW_rW_wd")))


#' For the retrospective evaluation of nowcasting we utilize all official data that was reported until July, 31 and restrict
#' the evaluation period until June, 30, assuming that for all days before, the true number of new cases with disease onset are reported
#' based on the available data at end of July. Furthermore, we focus on all reported cases with available disease onset information.
#' The aggregated case counts in this period are:

#+ message=FALSE, warning=FALSE, echo=FALSE, fig.height = 3.5, fig.width = 9, fig.cap="Daily numbers of newly diseased cases and reported cases in Bavarian data. Numbers of disease onsets are derived retrospectively based on data available on July, 31."
ggplot(dat_bav_smry) + geom_line(aes(date, n_dis), col = "green") +
  geom_line(aes(date, n_rep), col = "red", lty = 2) +
  theme_bw() + ylab("No. cases")


#' The empirical reporting delay distribution can be illustrated as in case of the synthetic data above:
 
#+ message=FALSE, warning=FALSE, echo=FALSE, fig.height = 3.5, fig.width = 9, fig.cap="Empirical reporting delay distribution (between disease onset and reporting at LGL) for the Bavarian COVID-19 data."
delay_true_smry = read_tsv("../../data_public/evaluation/dat_bavaria_truth_delay_smry.csv")
delay_true_smry %>%
  ggplot() + geom_line(aes(date, med_delay), col = "green") +
  geom_ribbon(aes(date, ymin = q25_delay, ymax = q75_delay), alpha=.25, fill="green") +
  coord_cartesian(xlim=c(ymd("2020-03-15"), ymd("2020-07-01"))) + 
  theme_bw() + ylab("Reporting delay")

#' Note that the illustrated delay is, in contrast to the delay reported in Table 2, the time between disease onset and reporting at
#' LGL (regional health authority). In Table 2, we described the delay between disease onset and reporting at the local health authority that is relevant for the
#' disease onset imputation model.

#' ### Results
#' 
#' We esimated nowcasts for all dates $t$ from March, 17 (22 days after disease onset of first case) until June, 30
#' by restricting the data to all cases reported until the respective date and compared the nowcast predicitions for all days
#' $t-5,\ldots, t-2$ to the numbers of newly diseased cases per day reported until July, 31.
#' Nowcasts are performed based on six different models (see description in the Paper), and the performance of the models
#' is compared as above. In addition we compute the coverage of the 95%-credibility intervals of the estimated time-varying 
#' reproduction number with the estimate obtained from utilizing all available data unti July, 31.

#+ message=FALSE, warning=FALSE, echo=FALSE, fig.height = 3.5, fig.width = 9
cov_full = left_join(ntInf, dat_bav_smry) %>%
  mutate(coverage = q025<=n_dis & n_dis <= q975) %>%
  filter(lag %in% 2:6) %>%
  group_by(model) %>% summarise(frac_cov = mean(coverage))

cov_may = left_join(ntInf, dat_bav_smry) %>%
  mutate(coverage = q025<=n_dis & n_dis <= q975) %>%
  filter(lag %in% 2:6, date <ymd("2020-05-01")) %>%
  group_by(model) %>% summarise(frac_cov = mean(coverage))

# Summarise Scoring Rules
crps_full = ntInf %>% filter(lag %in% c(2:6)) %>% group_by(model) %>% summarise(crps = mean(crps_score))
crps_may = ntInf %>% filter(lag %in% c(2:6),
                 date < ymd("2020-05-01")) %>% group_by(model) %>% summarise(crps = mean(crps_score))
logsc_full = ntInf %>% filter(lag %in% c(2:6)) %>% group_by(model) %>% summarise(logS = mean(log_score))
logsc_may = ntInf %>% filter(lag %in% c(2:6),
                             date < ymd("2020-05-01")) %>% group_by(model) %>% summarise(logS = mean(log_score))
rmse_full = ntInf %>% filter(lag %in% c(2:6)) %>%
  left_join(dat_bav_smry) %>%
  mutate(diff = med - n_dis) %>%
  group_by(model) %>% summarise(rmse = sqrt(mean(diff^2)))
rmse_may = ntInf %>% filter(lag %in% c(2:6),
                            date < ymd("2020-05-01")) %>%
  left_join(dat_bav_smry) %>%
  mutate(diff = med - n_dis) %>%
  group_by(model) %>% summarise(rmse = sqrt(mean(diff^2)))

load("../../data_public/evaluation/dat_bavaria_true_Rt.RData")
Rt = do.call(rbind, lapply(bav_files,
                           function(x) {
                             load(file = x)
                             do.call(rbind, lapply(res_list, function(x) x$Rt))
                           }))
Rt = Rt %>% mutate(lag=now-Date,
                   model = factor(model, levels = c("poisson_rW_const",
                                                    "poisson_rW_cp2W",
                                                    "negBinom_rW_cp2W",
                                                    "negBinom_rW_rW",
                                                    "negBinom_rW_cp2W_wd",
                                                    "negBinom_rW_rW_wd"))) %>%
  left_join(true_Rt %>% rename(Date=date,
                             Rt_true=Rt))

cov_Rt = Rt %>% filter(lag==11) %>% mutate(cov_95 = Rt_lower<=Rt_true & Rt_upper>=Rt_true) %>%
  group_by(model) %>% summarise(frac_cov_95 = mean(cov_95))

smry_bav_paper = tibble(
  DistrNtd = c("Poisson", "Poisson", rep("Neg. Binomial", 4)),
  DelayMod = c("No changes", "Lin. effect of time changepoints 2 weeks", "Lin. effect of time changepoints 2 weeks",
               "Daily changes (first order RW)",
               "Lin. effect of time changepoints 2 weeks + weekday effect",
               "Daily changes (first order RW) + weekday effect"),
  CRPS = crps_full$crps,
  logS = logsc_full$logS,
  RMSE = rmse_full$rmse,
  CovNtInf = cov_full$frac_cov,
  CovRt = round(cov_Rt$frac_cov_95,2))

kable(tibble(model = crps_full$model, crps_f = crps_full$crps, crps_m = crps_may$crps,
       logS_f = logsc_full$logS, logS_m = logsc_may$logS,
       rmse_f = rmse_full$rmse, rmse_m = rmse_may$rmse,
       cov_f = cov_full$frac_cov, cov_m = cov_may$frac_cov, 
       cov_Rt = cov_Rt$frac_cov_95), digits = 2, caption="Retrospective quantification of the performance of six different nowcast models on Bavarian COVID-19 data. Shown are the average metrics over all nowcast dates (*_f) and restricted to the period until May, 1 (*_m). Reported scores are the continuous ranked probability score, logarithmic score, root mean squared error of posterior median, and coverage frequencies of 95% prediction intervals, as well as coverage frequencies of the estimated R(t) at the most current date. Estimated models are Poisson and Negative Binomial models with 1) (assumed) constant delay distribution, 2) linear time-effects with changepoints every two weeks before now (*_cp2W), 3) daily changes in delay distribution based on a first-order random walk prior (*_rW) and 2) and 3) with additional effects of the weekday of case reporting (*_wd).")

#' The comparison of the daily nowcast predictions 3 day before now and the retrospective truth, as well 
#' as the estimated delay distribution are shown in the following figures:

#+ message=FALSE, warning=FALSE, echo=FALSE, fig.height = 7, fig.width = 6, fig.cap="Nowcast predictions and 95% PIs for nowcasts 3 days before 'now' compared to retrospectively true number of disease onsets (solid) and number of reported cases (dotted) in Bavarian data."
ntInf %>%
  dplyr::filter(lag==3) %>% ggplot() + geom_line(aes(date, med, col = model)) +
  geom_ribbon(aes(date, ymin=q025, ymax=q975, fill = model), alpha = .3) +
  geom_line(aes(date, n_dis), dat_bav_smry) +
  geom_line(aes(date, n_rep), dat_bav_smry, lty = 2) +
  theme_bw() +
  facet_grid(model~.) +
  coord_cartesian(xlim = c(ymd("2020-02-24"), ymd("2020-05-01")), ylim = c(0,2500)) +
  theme(legend.title = element_blank(), legend.pos = "bottom") + ylab("Number cases")
  
#'

#+ message=FALSE, warning=FALSE, echo=FALSE, fig.cap="Comparison of estimated delay distribution (black) and retrospective empirical delay distribution (blue) per day for each model on Bavarian data. Shown are the median and 25% and 75% quantiles."
p = do.call(rbind, lapply(bav_files,
                          function(x) {
                            load(x)
                            do.call(rbind, lapply(res_list, function(x) tryCatch(x[[2]] %>% filter(date==max(date)),
                                                                                 error = function(e) rep(NA, 24))))
                          }))
p = p %>% filter(!is.na(date)) %>% arrange(date)
p_long = p %>% pivot_longer(cols = starts_with("D", ignore.case = F))
delay_est = p_long  %>% group_by(date, model) %>%
  summarise(med_delay = min(which(cumsum(value)>=.5))-1,
            q10_delay = min(which(cumsum(value)>=.10))-1,
            q25_delay = min(which(cumsum(value)>=.25))-1,
            q40_delay = min(which(cumsum(value)>=.40))-1,
            q60_delay = min(which(cumsum(value)>=.60))-1,
            q75_delay = min(which(cumsum(value)>=.75))-1,
            q90_delay = min(which(cumsum(value)>=.90))-1
  )
delay_est = delay_est %>% mutate(model = factor(model,
                                                levels = c("poisson_rW_const",
                                                           "poisson_rW_cp2W",
                                                           "negBinom_rW_cp2W",
                                                           "negBinom_rW_rW",
                                                           "negBinom_rW_cp2W_wd",
                                                           "negBinom_rW_rW_wd")))
delay_true = read_tsv("../../data_public/evaluation/dat_bavaria_truth_delay_smry.csv")

ggplot() +
  geom_ribbon(aes(date, ymin=q25_delay, ymax=q75_delay), delay_est, alpha = .2) +
  geom_line(aes(date, med_delay), delay_est) +
  facet_wrap(~model, ncol = 2) + theme_bw() +
  geom_ribbon(aes(date, ymin=q25_delay, ymax=q75_delay), delay_true, alpha = .2, fill = "green") +
  geom_line(aes(date, med_delay), delay_true, col = "green", lwd = 1.5) + ylab("Reporting delay")


#' We find that the Poisson model assuming no changes in the reporting delay distribution performs bad. Daily case numbers are strongly
#' overestimated.
#' This is in line with the apparent changes in the reporting delay between disease onset and reporting at LGL over time. 
#' Comparing the Poisson model with two-week changepoints with a similar model using a Negative Binomial distribution for 
#' $N_{t,d}$ we find the latter to perform better with respect to the evaluation metrics. Furthermore, it is apparent that
#' prediction intervals in the poisson model are too narrow. This improves when using a Negative Binomial model with overdispersion.
#' Adding weekday effects to the delay distribution improves the performance of the models as well. Comparing the Negative Binomial model with daily changes in the 
#' delay distribution with the two week changepoint model, we found better coverage frequencies for the former (mainly because of wider prediction intervals) 
#' but lower CRPS score and RMSE for the latter. The reporting delay estimation allowing for daily changes based on a first-order random walk 
#' appeared to be somewhat instable at the beginning of the nowcasting period, as also seen in the evaluation based on synthetic data.
#' 
#' Looking visually at the predictions of the nowcast, we find that the predictions close to now (e.g. 3 days before now), as illustrated
#' in the Fig. 7, overestimate daily case counts in the crucial timeperiod between March 15 and April 1. This is also true for our preferred model,
#' with 2 week changepoints and weekday effects in the delay distribution. In this period the case counts of newly
#' diseased individuals stabilized, but daily reported cases were still increasing steadily. The situation of a stabilizing number of new disease onsets and simultaneously decreasing average
#' reporting delay between disease onset and case registration (that becomes now - in retrospect- apparent) is not easily identifiable 
#' based on incomplete data in real time surveillance. 
#' It is, however, apparent that the nowcasting approach is valuable in understanding the
#' dynamics of the pandemic better compared to the daily counts of newly reported cases. It helps to illustrate uncertainty 
#' with respect to the current state of the pandemic and the predictions of our preferred models that account for changes in reporting delay are not too far off. Furthermore coverage frequencies of 95% prediction intervals 
#' close to 90% appear acceptable. Looking more closely at the results of e.g., the Negative Binomial model with 2-week changepoints and weekday
#' effects in the delay distribution, we find that based on the nowcast, the pandemic situation seems to stabilize from around March, 20 on, and the predictions
#' close to now are already starting to decreas at around April, 1, where the daily number of newly reported cases were still at their peak.
#' 
#' 
#' Comparing the estimated $R(t)$ at most current $t$ based on the different nowcast models with the retrospective \textit{truth} based on all reported data, we find coverage probabilities of the $95\%-$ credibility intervals bigger 
#' than $90\%$ for all models that consider changes in the delay distribution over time. 
#' The estimation of $R(t)$ is, however, biased when it is based on a biased nowcasting approach, e.g., due to 
#' ignoring changes in the delay distribution. In the Poisson model assuming no changes in delay distribution $R(t)$ biased upwards
#' during the whole timeperiod. 

#+ message=FALSE, warning=FALSE, echo=FALSE, fig.height = 3.5, fig.width = 9, fig.cap="Estimated R(t) over time based on all nowcast models at most current date and associated 95% CI. Comparison with R(t) based on all reported disease onsets until July, 31."
Rt %>% filter(lag == 11) %>% ggplot() +
  geom_line(aes(Date, Rt, col = model)) +
  geom_ribbon(aes(Date, ymin = Rt_lower, ymax = Rt_upper, fill = model), alpha = .1) +
  geom_line(aes(Date, Rt_true), lty = 2, lwd = 1.2) +
  theme_bw() + theme(legend.title = element_blank(), legend.position = "bottom")


