# Load required packages and helper functions
library(tidyverse)
library(segmented)
library(future.apply)
source("./breakpoint_fun.R")

# Load data
dat = read_tsv("../results_public/nowcasting_results_bavaria_2020-05-05.csv")
# Adjust data to contain only columns 'date' and 'nc'
dat = dat %>% mutate(date = date, nc = bayes_wd_tps_predicted) %>%
  dplyr::select(date, nc) %>%
  filter(date<=lubridate::ymd("2020-05-01"))

# Estimate breakpoint model based on discrete optimization
plan("multisession")
options(mc.cores = 7)
res_disc4 = estimate_bp_disc_optim(data = dat, bp = 4)
# Save Results of discrete optimization
save(res_disc4, file = "../results_public/res_disc_rki_2020-05-05_4bp.RData")
load("../results_public/res_disc_rki_2020-05-05_4bp.RData")
# Estimate breakpoint model based on segmented 'function' with starting values
# from discrete optimization n.boot = 50, it.max = 500, seed = 1495
res_segmented_4 = estimate_bp_segmented(data = dat, bp = 4, start_bp = res_disc4$bps,
                                      plot_main = "Tägliche Neuerkrankungen Bayern",
                                      n_boot = 50, segmented_seed = 1495)
# Plot results breakpoint analysis
res_segmented_4$plot
