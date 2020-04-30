# Load required packages and helper functions
library(tidyverse)
library(segmented)
library(future.apply)
source("./breakpoint_fun.R")
# Load data
dat = read_tsv("../results_public/nowcasting_results_2020-04-09.csv")
# Adjust data to contain only columns 'date' and 'nc'
dat = dat %>% mutate(nc = bayes_wd_tps_predicted) %>%
  dplyr::select(date, nc)

# Estimate breakpoint model (with three breakpoints) based on discrete optimization
plan("multicore")
options(mc.cores = 6)
res_disc = estimate_bp_disc_optim(data = dat, bp = 3)

# breakpoints with lowest deviance
res_disc$bps
res_disc$deviance
res_disc$overdispersion

# Estimate breakpoint model based on 'segmented' function with starting values
# from discrete optimization
res_segmented = estimate_bp_segmented(data = dat, bp = 3, start_bp = res_disc$bps,
                                      plot_main = "Tägliche Neuerkrankungen Bayern",
                                      n_boot = 50, segmented_seed = 1495)
res_segmented$plot
res_segmented$coef
res_segmented$breakpoints


