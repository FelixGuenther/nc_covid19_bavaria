# Load packages
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
library(ggpubr)


source("./analysis_fun.R")

# Read data
dat = read_tsv("../../data_public/main_analysis/synth_dat.csv")

# Edit data set
set.seed(241)
dat_edit = dat %>%
  mutate(rep_date_local = pmax(rep_date_reg-sample(0:4, size = n(), replace = T, prob = c(.1,.4,.3,.1,.1)),
                               disease_start, na.rm = T),
         rep_date_reg_weekday = weekdays(rep_date_reg),
         rep_date_local_weekday  = weekdays(rep_date_local),
         rep_week_reg = week(rep_date_reg),
         rep_week_local = week(rep_date_local),
         delay = as.numeric(rep_date_local - disease_start))

# Summarise available delay data
# Table 2
dat_smry = summarize_data(dat_edit)
write.table(dat_smry, file = "../../results_public/1_main_analysis/tables/table_2_delay_smry_synth.csv", sep="\t", row.names = FALSE, quote = FALSE)
# Impute data based on Weibull GAMLSS
imputation_res = perform_imputation(dat_edit = dat_edit,
                                    type = "week_weekday_age")
# Plot imputation results
# Figure 1
fig_1 = imputation_res[[3]] + theme(axis.text=element_text(size=12),
                                axis.title=element_text(size=12), 
                                legend.text = element_text(size=12), legend.title = element_text(size=12))
ggsave(fig_1, filename = "../../results_public/1_main_analysis/figures/fig_1_imp_synth.jpg", 
       width = 10, height = 3.5, dpi = 250)


# Get imputated data
imputed_data = imputation_res[[1]]

# Estimate Nowcasts
# Define current date
data_date = ymd("2020-04-09")
# Compute nowcasts (requires working rjags installation on device),
# Results can be automatically saved (see documentation in './analysis_fun.R')
nc_res = estimate_nowcasts(imputed_data = imputed_data,
                           data_date = data_date,
                           safePredictLag = 2, 
                           save_results = FALSE)

# Figure 2
ntInf_true = read_tsv("../../data_public/main_analysis/onsets_true_retrospec_07-31.csv")
plot_dat = left_join(nc_res$nc_smry$ntInf, 
                     imputed_data %>% 
                       filter(rep_date_reg<=data_date-1) %>%
                       group_by(date=disease_start_imp) %>% 
                       summarise(obs_imp_onset=n())) %>%
  left_join(imputed_data %>% 
                       filter(rep_date_reg<=data_date-1) %>% 
                       group_by(date=disease_start) %>% 
                       summarise(obs_onset=n())) %>%
  left_join(imputed_data %>% 
              filter(rep_date_reg<=data_date-1) %>% 
              group_by(date=rep_date_reg) %>% summarise(reported_new_cases=n())) %>%
  mutate(obs_onset=replace_na(obs_onset, 0),
         reported_new_cases = replace_na(reported_new_cases, 0))

fill_colours = c("med" = "#018571",
                 "obs_onset" = "lightblue",
                 "obs_imp_onset" = "grey",
                 "onset_retr" = "black",
                 "reported" = "darkred")
fill_labels = c("med" = "Nowcast",
                "obs_onset" = "Disease start \nobserved",
                "obs_imp_onset" = "Disease start \nimputed",
                "onset_retr" = "Disease start\nretrospective",
                "reported" = "Reported")

fig_2 = ggplot(plot_dat %>% filter(date<=ymd("2020-04-06"))) +
  scale_x_date(breaks = seq(max(plot_dat$date), min(plot_dat$date), by = "-1 week"), date_labels = "%d.%m.") +  
  xlab("Date") +
  ylab("Number cases") +
  coord_cartesian(xlim = c(ymd("2020-02-24"), ymd("2020-04-08"))) +
  geom_col(aes(date, obs_imp_onset, fill = "obs_imp_onset")) +
  geom_col(aes(date, obs_onset, fill = "obs_onset")) +
  geom_line(aes(date, reported_new_cases), col = fill_colours["reported"]) +
  geom_line(aes(date, ntInf_true), data = ntInf_true, 
            col = fill_colours["onset_retr"], lty = 2) +
  geom_ribbon(aes(date, ymin = reported_new_cases, ymax = reported_new_cases, fill = "reported", 
                  col = "reported"), 
              alpha = .25, lwd = .25) +
  geom_ribbon(aes(date, ymin = ntInf_true, ymax = ntInf_true, fill = "onset_retr", col = "onset_retr"), 
              alpha = 0, lwd = 0, data = ntInf_true) +
  geom_line(aes(date, med), col = fill_colours["med"]) +
  geom_ribbon(aes(date, ymin = q025, ymax = q975, fill = "med", col = "med"), 
              alpha = .25, lwd = .25) +
  scale_color_manual(values = fill_colours,
                     labels = fill_labels) +
  scale_fill_manual(values = fill_colours,
                    labels = fill_labels) +
  theme_bw() +
  theme(legend.title = element_blank(), legend.position = "bottom") +
  geom_point(aes(x=ymd("2020-04-08"), y = 0), color = "gold", pch = 14,
             size = 1.5, show.legend = FALSE) +
  guides(linetype=FALSE, color = FALSE) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12))

ggsave(fig_2, filename = "../../results_public/1_main_analysis/figures/fig_2_nc_synth.jpg", dpi = 250, width = 12, height = 5)

# Estimate Rt
Rt = estimate_Rt(nc_res, data_date)

# Figure 3
fig_3 = ggplot(Rt, aes(x = Date, y = Rt)) +
  geom_ribbon(aes(ymin = Rt_lower, ymax = Rt_upper), fill = "lightgray", alpha = .9) +
  geom_line() +
  xlab("Time") + ylab(expression(R(t))) +
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
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12))
ggsave(fig_3, filename = "../../results_public/1_main_analysis/figures/fig_3_Rt_synth.jpg", width = 9, height = 4, dpi = 500)

# Sensitivity analysis
# Reduce dataset to 1) all persons with explicit information on symptoms; 2) exclusion of all persons explicitly without symptoms
dat_sens_1 = dat_edit %>% filter(sympt_expl)
dat_sens_2 = dat_edit %>% filter(!no_sympt_expl)
# Perform imputation
set.seed(127744)
imputation_sens1 = perform_imputation(dat_sens_1, type = "week_weekday_age")
set.seed(18457)
imputation_sens2 = perform_imputation(dat_sens_2, type = "week_weekday_age")

imputed_data_sens1 = imputation_sens1[[1]]
imputed_data_sens2 = imputation_sens2[[1]]

# Estimate nowcasts
set.seed(1835)
nc_res_sens1 = estimate_nowcasts(imputed_data = imputed_data_sens1, data_date = data_date, 
                                 save_results = FALSE, safePredictLag = 2)

set.seed(22374)
nc_res_sens2 = estimate_nowcasts(imputed_data = imputed_data_sens2, data_date = data_date, 
                                 save_results = FALSE, safePredictLag = 2)

Rt_sens1 = estimate_Rt(nc_res_sens1, data_date = data_date)
Rt_sens2 = estimate_Rt(nc_res_sens2, data_date = data_date)

# Plot results (Figure 4)
plot_dat = rbind(nc_res$nc_smry$ntInf %>% dplyr::select(date, med, q025, q975) %>% mutate(type = "Main"),
                 nc_res_sens1$nc_smry$ntInf %>% dplyr::select(date, med, q025, q975) %>% mutate(type = "Sensitivity 1"),
                 nc_res_sens2$nc_smry$ntInf %>% dplyr::select(date, med, q025, q975) %>% mutate(type = "Sensitivity 2"))

plot_dat_rt = rbind(Rt %>% dplyr::select(date=Date, Rt, Rt_upper, Rt_lower) %>% mutate(type = "Main"),
                    Rt_sens1 %>% dplyr::select(date=Date, Rt, Rt_upper, Rt_lower) %>% mutate(type = "Sensitivity 1"),
                    Rt_sens2 %>% dplyr::select(date=Date, Rt, Rt_upper, Rt_lower) %>% mutate(type = "Sensitivity 2"))
# Nowcasts
sens_nt = ggplot(plot_dat) + geom_line(aes(date, med, col = type)) + 
  geom_ribbon(aes(date, ymin = q025, ymax = q975, fill = type), alpha=.3) +
  xlab("Date") +
  ylab("Number cases") +
  coord_cartesian(xlim = c(ymd("2020-02-24"), ymd("2020-04-08"))) +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12)) + theme(legend.title = element_blank(), legend.position = "bottom") +
  guides(linetype=FALSE, color = FALSE)

# R(t)
sens_r = ggplot(plot_dat_rt) + geom_line(aes(date, Rt, col = type)) +
  geom_ribbon(aes(date, ymin=Rt_lower, ymax = Rt_upper, fill = type), alpha=.3) + 
  xlab("Time") + ylab(expression(R(t))) +
  geom_hline(yintercept=1, lty=2) +
  scale_y_continuous(breaks=seq(0,6, by=1)) +
  scale_x_date(breaks = seq(max(Rt$Date), min(Rt$Date), by = "-1 week"), 
               labels = paste0(strftime(seq(max(Rt$Date), min(Rt$Date), by = "-1 week"), 
                                        format = "%d.%m."), 
                               "-\n", 
                               strftime(seq(max(Rt$Date), min(Rt$Date), by = "-1 week")+10, 
                                        format = "%d.%m."),
                               "  ")) +  
  coord_cartesian(ylim = c(0,5)) +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12)) + theme(legend.title = element_blank(), legend.position = "bottom") +
  guides(linetype=FALSE, color = FALSE)
# Combine to Figure 4
fig_4 = ggarrange(sens_nt + theme(legend.position = "none"), 
                  sens_r + theme(legend.position = "right"), ncol = 2, labels = "AUTO")
ggsave(fig_4, filename = "../../results_public/1_main_analysis/figures/fig_4_sens_synth.jpg", dpi=250, width = 10, height = 3.5)
