library(Rcpp)
library(surveillance)
library(rstan)
library(lubridate)
library(tidyverse)
library(nleqslv)
library(scoringRules)
library(R0)
library(future)
library(future.apply)

# Load data
load("../../../data_public/evaluation/dat_bavaria_persspec_synth.RData")
# Load model
mod = readRDS("../stan_models/negBinom_rWalk_cp.rds")
# Load functions to estimate nowcast
source("../eval_fun.R")

# Define dates for which to estimate nowcast
eval_dates = seq(ymd("2020-02-24") + 22, ymd("2020-06-30"), by = 1)
# Restrict to fewer dates for faster computation
eval_dates = eval_dates[1]

# Setup cluster to perform estimation
# Library Path
libs <- .libPaths()[1]
libs_cmd <- sprintf(".libPaths(c(%s))", paste(sprintf('"%s"', libs), collapse = ", "))

# Define number of workers
# Note that each call to estimate_nowcast starts four parallel chains of MCMC sampling 
# (Available number of cores has to be 4 times the number of workers in cluster)
cl <- future::makeClusterPSOCK(1, rscript_args = c("-e", shQuote(libs_cmd)))
plan(cluster, workers = cl)
# Apply over dates to perform nowcasting based on data available until 'now'
res_list = future_lapply(sample(eval_dates), FUN = function(x) {
  print(paste0("Start nowcast: ",x))
  tryCatch(estimate_nowcast(now=x, 
                            data = dat_mod, 
                            begin_date = ymd("2020-02-24"),
                            D=21, 
                            predLag = 2, 
                            model = mod, 
                            mod_name = "negBinom_rW_cp2W"),
           error = function(e) e)
})

save(res_list, file = "../../../results_public/2_evaluation/bavaria/negBinom_rW_cp_2w_synthetic.RData")
