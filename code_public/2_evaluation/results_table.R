library(lubridate)
library(tidyverse)
# Eva
# Load 'True' synthetic data
load("../../data_public/evaluation/dat_mod_synthetic.RData")
dat_mod_synth = dat_mod
rm(dat_mod)

# Load resuls of models applied to synthetic data
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


ntInf_synth = ntInf_synth %>% mutate(lag = as.numeric(now-date),
                                     model = factor(model, levels = c("poisson_rW_const",
                                                                      "poisson_rW_cp2W",
                                                                      "negBinom_rW_cp2W",
                                                                      "negBinom_rW_rW",
                                                                      "poisson_rW_cptrue",
                                                                      "negBinom_rW_cptrue")))

# Quantitative Summary
cov_full = left_join(ntInf_synth, dat_mod_synth %>% group_by(date = disease_start) %>% summarise(n_dis = n())) %>%
  mutate(coverage = q025<=n_dis & n_dis <= q975) %>% 
  filter(lag %in% 2:6) %>%
  group_by(model) %>% summarise(frac_cov = mean(coverage))



# Summarise Scoring Rules
crps_full = ntInf_synth %>% filter(lag %in% c(2:6)) %>% group_by(model) %>% summarise(crps = mean(crps_score))
logsc_full = ntInf_synth %>% filter(lag %in% c(2:6)) %>% group_by(model) %>% summarise(logS = mean(log_score))
rmse_full = ntInf_synth %>% filter(lag %in% c(2:6)) %>% 
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


smry_synth_paper

# Bavarian data
# True data summary
dat_bav_smry = read_tsv("../../data_public/evaluation/dat_bavaria_truth_epcurve_reported.csv")

# Load model results
bav_files = c("../../results_public/2_evaluation/bavaria/poiss_rW_const.RData",
              "../../results_public/2_evaluation/bavaria/poiss_rW_cp_2w.RData",
              "../../results_public/2_evaluation/bavaria/negBinom_rW_rW.RData",
              "../../results_public/2_evaluation/bavaria/negBinom_rW_rW_wd.RData",
              "../../results_public/2_evaluation/bavaria/negBinom_rW_cp_2w.RData",
              "../../results_public/2_evaluation/bavaria/negBinom_rW_cp_2w_wd.RData")

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
# Compute coverage frequencies of 95%-PIs
cov_full = left_join(ntInf, dat_bav_smry) %>%
  mutate(coverage = q025<=n_dis & n_dis <= q975) %>%
  filter(lag %in% 2:6) %>%
  group_by(model) %>% summarise(frac_cov = mean(coverage))

# Summarise Scoring Rules
crps_full = ntInf %>% filter(lag %in% c(2:6)) %>% group_by(model) %>% summarise(crps = mean(crps_score))
logsc_full = ntInf %>% filter(lag %in% c(2:6)) %>% group_by(model) %>% summarise(logS = mean(log_score))
rmse_full = ntInf %>% filter(lag %in% c(2:6)) %>%
  left_join(dat_bav_smry) %>%
  mutate(diff = med - n_dis) %>%
  group_by(model) %>% summarise(rmse = sqrt(mean(diff^2)))

# Coverage frequencies of R(t) estimation
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

# Full table
table_3_full = rbind(smry_synth_paper, smry_bav_paper)
write.table(table_3_full, file = "../../results_public/2_evaluation/tables/table_3_eval.csv", dec = ".", sep= "\t", row.names = FALSE, quote = FALSE)
