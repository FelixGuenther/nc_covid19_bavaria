library(tidyverse)
library(data.table)
library(R0)
library(nleqslv)
library(lubridate)
library(broom)
library(future)
library(gamlss)
library(rstan)
library(surveillance)

source("./analysis_fun.R")

# Read data
dat = read_tsv("../data_public/data_public.csv")

# Edit data set
set.seed(241)
dat_edit = dat %>% rename(rep_date_reg = rep_date) %>%
  mutate(rep_date_local = pmax(rep_date_reg-sample(0:4, size = n(), replace = T, prob = c(.1,.4,.3,.1,.1)),
                               disease_start, na.rm = T),
         rep_date_reg_weekday = weekdays(rep_date_reg),
         rep_date_local_weekday  = weekdays(rep_date_local),
         rep_week_reg = week(rep_date_reg),
         rep_week_local = week(rep_date_local),
         delay = as.numeric(rep_date_local - disease_start))

# Summarise available delay data
dat_smry = summarize_data(dat_edit)
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
# Results can be automatically be safed (see documentation in './analysis_fun.R')
nc_res = estimate_nowcasts(imputed_data = imputed_data,
                           data_date = data_date,
                           safePredictLag = 2, 
                           save_results = FALSE)
# Estimate Rt
Rt = estimate_Rt(nc_res, data_date)

# Plot Nowcast results
res_plots = create_plots(nc_res$nc_smry$ntInf, Rt = Rt, imputed_data = imputed_data, 
                         now = nc_res$add_inf$now, 
                         safePredictLag = nc_res$add_inf$safePredictLag)
# Plot nowcast results
res_plots[[1]]
res_plots[[2]]
res_plots[[3]]

# Save nowcast results
nc_res_dat = nc_res$nc_smry$ntInf
write.table(nc_res_dat, file = paste0("../results_public/nowcasting_results_", data_date,".csv"), 
            dec = ".", sep = "\t", quote = FALSE, row.names = FALSE)

