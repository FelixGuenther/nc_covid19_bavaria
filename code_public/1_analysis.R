library(tidyverse)
library(surveillance)
library(lubridate)
library(broom)
source("./analysis_fun.R")

# Read data
dat = read_tsv("../data_public/data_public.csv")

# Edit data set
dat_edit = dat %>% mutate(rep_date_weekday = weekdays(rep_date),
                          rep_week = as.numeric(strftime(rep_date, format = "%V")),
                          delay = as.numeric(rep_date - disease_start))

# Summarise available delay data
dat_sum = summarize_data(dat_edit)

# Impute data based on Weibull GAMLSS
imputation_res = perform_imputation(dat_edit = dat_edit,
                                    type = "week_weekday_age")
# Plot imputation_results
# Figure 1
imputation_res[[3]]
# Get imputated data
imputed_data = imputation_res[[1]]

# Estimate Nowcasts
# Define current date
data_date = ymd("2020-04-09")
# Compute nowcasts (requires working rjags installation on device),
# Results and diagnostic plots can automatically be safed (see documentation in './analysis_fun.R')
nc_res = estimate_nowcasts(imputed_data = imputed_data,
                           data_date = data_date,
                           safePredictLag = 2, 
                           save_results = FALSE)
# Plot Nowcasts results
# Fig 2
nc_res[[4]]
# Fig 3
nc_res[[5]]

## Estimation R0
library(R0)
library(EpiEstim)

nc_info <- nc_res$add_inf
# Extract nowcast object
nc <- nc_res$nc_obj$nc_bayes_tps_wd

# Generation time from Nishiura et al. (2020)
gt_ni <- c(mean = 4.70, sd = 2.90)
gt <- R0::generation.time("lognormal", gt_ni)


set.seed(131015)
source("general/est.R0.TD.R") # customized version of R0::est.R0.TD
source("./general/functions-R0.R")
# extract nowcast dates
now <- nc_info$now - nc_info$safePredictLag
#Time the 2nd case - could remove imported cases more explictly
time_2ndcase <- epoch(nc)[min(which(cumsum(observed(nc)) >= 2))]

#Number of days of the serial interval to ensure
most_secondary_transmissions_occurred <- gt$time[min(which(cumsum(gt$GT)
                                                           >= 0.95))]
# End of R(t) estimation
end_Rt_estim <- now - most_secondary_transmissions_occurred + 1

# Estimate R0 using the Wallinga & Teunis method.
# sample time-series from posterior of now-cast object
# each row is one sample from posterior
posterior_samples <- sample_from_nc_posterior(nc, n_sample = 1000)

library(future)
plan("multicore")
# Apply R0 estimatio to to each sampled timeseries of new daily cases
options(mc.cores = 6)
Rt_est_list <- furrr::future_map(
  .x = seq_len(nrow(posterior_samples)),
  .f = ~{
    ts <- posterior_samples[.x,]
    est.R0.TD( # use custom function that returns R.simu
      epid    = ts,
      GT      = gt,
      t       = epoch(nc),
      correct = TRUE,
      begin   = time_2ndcase,
      end     = end_Rt_estim,
      nsim    = 1000)
  }
)

# Combine R0 estimates into one data frame
Rt_df <- map_dfr(
  Rt_est_list,
  ~ tidy_Rsimu(.x, nc, time_2ndcase, end_Rt_estim), .id="posterior_sample") %>%
  gather(iter, Rt, -posterior_sample, - Date)
# Aggregate over samples
Rt_df_smry <- Rt_df %>%
  group_by(Date) %>%
  summarize(
    Rt_lower = quantile(Rt, .025),
    Rt_upper = quantile(Rt, .975),
    Rt = mean(Rt))

# get R(t) at most current date
last_Rt <- Rt_df_smry %>% filter(Date == max(Date))
# Plot R(t)
gg_Rt <- ggplot(Rt_df_smry, aes(x = Date, y = Rt)) +
  geom_ribbon(aes(ymin = Rt_lower, ymax = Rt_upper), fill = "lightgray", alpha = .9) +
  geom_line() +
  xlab("Date") + ylab(expression(R(t))) +
  geom_hline(yintercept=1, lty=2) +
  scale_y_continuous(breaks=seq(0,6, by=1)) +
  scale_x_date(date_breaks = "2 week") +
  coord_cartesian(ylim = c(0,5)) + theme_bw()
gg_Rt
