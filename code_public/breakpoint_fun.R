#' Function to estimate a piecewise-linear model quasipoisson model fort the analysis
#' of new COVID-19 cases per day. Breakpoints are found based on discrete optimization
#' over all potential combinations of breakpoints with at least 4 days in between 2 breakpoints
#' @param data dataset of subsequent days including a column 'date' and a column 'nc' which contains
#' nowcast results for the estimated number of new cases per day
#' @param bp number of breakpoints (larger numbers than 4 require long computation time and available memory)
#' @return list with entries: 'mod_min_dev': the fitted quasi-poisson glm model with minimum deviance; 'bps': the breakpoints 
#' of the model with minimum deviance (in days), 'deviance': the deviance of the model with minimum deviance; 
#' 'overdispersion': overdispersion of model with minimum deviance; 'res_disc_opt': tibble with deviance and overdispersion
#' for each combination of breakpoints
estimate_bp_disc_optim = function(data, bp=3) {
  dat = data %>% mutate(t=1:n())
  day_vec = lapply(1:bp, function(x) 4:(nrow(data)-4))
  bp_scen = tibble(bp_1 = day_vec[[1]])
  for(i in 2:bp) {
    bp_scen = expand_grid(bp_scen, !!paste0("bp_", i) := day_vec[[i]]) %>%
      dplyr::filter(!!as.symbol(paste0("bp_", i)) > !!as.symbol(paste0("bp_", i-1)),
                    (!!as.symbol(paste0("bp_", i)) - !!as.symbol(paste0("bp_", i-1))) > 2)
  }
  create_data_est_model = function(bp_vec, data = dat) {
    # Create data matrix, linear trends between changepoints
    dat_mod = cbind(data %>% dplyr::select(nc, t),
                    do.call(cbind, lapply(1:length(bp_vec), 
                                          function(x) pmax(0, data$t - bp_vec[x] + 1))))
    colnames(dat_mod) = c("nc", "t", paste0("t", 1:length(bp_vec)))
                            
    # Estimate quasipoisson model
    mod = tryCatch(glm(nc~., family = quasipoisson, data = dat_mod), error = NULL)
    # Out: changepoints and corresponding deviance
    c(bp_vec, deviance = ifelse(is.null(mod), NA, mod$deviance), 
      dispersion = ifelse(is.null(mod), NA, summary(mod)$dispersion))
  }
  res = future_apply(bp_scen, 1, create_data_est_model)
  res = as_tibble(t(res))
  # Get Breakpoints with lowest deviance
  bp_ts = res %>% dplyr::filter(rank(deviance) == 1)
  dat_mod = cbind(dat %>% dplyr::select(nc, t),
                  do.call(cbind, lapply(1:bp, 
                                        function(x) pmax(0, dat$t - unlist(bp_ts[1,1:bp])[x] + 1))))
  colnames(dat_mod) = c("nc", "t", paste0("t", 1:bp))
  # Re-estimate model
  mod = glm(nc~., family = quasipoisson, data = dat_mod)
  
  list(mod_min_dev = mod, bps = unlist(bp_ts[1:bp]),
       deviance = mod$deviance, overdispersion = summary(mod)$dispersion, res_disc_opt = res)
}

#' Function to estimate a piecewise-linear model quasipoisson model fort the analysis
#' of new COVID-19 cases per day using the 'segmented' function of the 'segmented' package.
#' @param data dataset of subsequent days including a column 'date' and a column 'nc' which contains
#' nowcast results for the estimated number of new cases per day
#' @param bp number of breakpoints (larger numbers than 4 require long computation time and available memory)
#' @param start_bp starting values for the estimated breakpoints (might be derived from discrete optimization), length
#' of the vector has to correspond with the 'bp' parameter
#' @param plot_main desired title of the computed ggplot object
#' @param n_boot number of bootstrap samples in the bootstrap restarting algorithm of the 'segmented' function. Default to 50.
#' Bootstrapping can help to escape local minima of the objective function. Depending on the data situation, optimization can
#' fail to converge for few or several of the bootstrap samples. Please check resulting warning messages and results 
#' carefully and potentially adjust the number of breakpoints. Changing starting values and seed can help to check stability 
#' of results.
#' @param segmented_seed seed used in the 'segmented' function for reproducibility of (bootstrap) results.
#' @return list with entries: 'segmented_model': the result object of the 'segmented' funtion, 
#' 'coef': tibble with information on multiplicative effects on daily new cases in the segments of the
#' estimated epidemic curve, 'breakpoints': tibble of the estimated breakpoints, 'plot': ggplot
#' of the estimated epidemic curve


estimate_bp_segmented = function(data, bp = 3, start_bp = c(18,23,54), 
                                 plot_main = "Tägliche Neuerkrankungen",
                                 n_boot = 50,
                                 segmented_seed = 1135235) {
  stopifnot(length(start_bp)==bp)
  dat = data %>% mutate(t=1:n())
  mod1 = glm(nc~t, family=quasipoisson, data=dat)
  segmented = segmented(mod1, seg.Z=~t,
                        npsi = 3, psi = start_bp,
                        control = seg.control(n.boot = n_boot, it.max = 1000, seed = segmented_seed))
  
  # Results
  # get regression coefficients
  gamma = summary(segmented)$coefficients[2:(2+bp),1]
  # linear effect on log scale per interval
  st = cumsum(gamma)
  # get varcov of gammas
  m = vcov(segmented)[2:(2+bp),2:(2+bp)]
  # derive variance of linear combinations / linear effect per interval
  lin_comb = lower.tri(matrix(rep(0,length(st)^2), ncol=length(st)), diag = T)
  var = diag(lin_comb %*% m %*% t(lin_comb))
  # multiplicative effects on cases per day and 95%-CI
  mod_res = tibble(factor = exp(st), factor_CI_lwr = exp(st - 2*sqrt(var)),
                   Factor_CI_upr = exp(st+2*sqrt(var)))
  # breakpoints
  break_points <- as_tibble(confint(segmented)) %>% 
    rename(bp=Est., bp_lwr=`CI(95%).low`, bp_upr = `CI(95%).up`) %>%
    mutate(bp_date = min(dat$date) -1 + round(bp), 
           bp_lwr_date = min(dat$date)-1 + floor(bp_lwr), 
           bp_upr_date = min(dat$date)-1 + ceiling(bp_upr))
  
  bp_dates_infect = break_points %>% mutate(bp_inf = bp_date - 5,
                                          bp_lwr_inf = bp_lwr_date-5,
                                          bp_upr_inf = bp_upr_date-5)
  bp_dates_res = bp_dates_infect %>%
    mutate(BP =
             paste0(as.character(round(bp,1)), " (", as.character(bp_date), ")"),
           BP_CI_lwr =
             paste0(as.character(round(bp_lwr,1)), " (", as.character(bp_lwr_date), ")"),
           BP_CI_upr =
             paste0(as.character(round(bp_upr,1)), " (", as.character(bp_upr_date), ")")) %>%
    dplyr::select(BP, BP_CI_lwr, BP_CI_upr)
  
  dat = dat %>% mutate(pred_seg = segmented$fitted.values)
  bp_plot = ggplot() +
    geom_col(data = dat, aes(x=t, y=nc),
             col = "grey", fill = "lightgrey") +
    geom_rect(aes(xmin=break_points$bp_lwr,
                  xmax=break_points$bp_upr, ymin=0, ymax=Inf), fill = "steelblue", alpha = 0.25)+
    geom_line(aes(x=seq(0,max(dat$t), 0.1), 
                  y=predict(segmented, newdata = data.frame(t=seq(0,max(dat$t), 0.1)),
                            type = "response")), 
              col = "black", lwd = 1.1) +
    geom_vline(xintercept=break_points$bp, lty=2, col ="steelblue", size=0.8) +
    theme_bw()+
    scale_x_continuous(breaks = seq(max(dat$t), min(dat$t), by = -14), 
                       labels = strftime(seq(max(dat$date), min(dat$date), by = "-2 weeks"), 
                                         format = "%d.%m.")) +
    scale_y_continuous(expand=expansion(mult = c(0, .1))) +
    geom_text(mapping=aes(x=break_points$bp,  y = 20, label = strftime(break_points$bp_date, format = "%d.%m.")), 
              size = 5, angle = 90, vjust =-.5, hjust = -0.1) +
    labs(x = "Zeit", y = "Anzahl Fälle (geschätzt)", title = plot_main,
         caption = "Statistisches Beratungslabor StaBLab, LMU München;  Department of Mathematics, Stockholm University \nDaten: Bayerisches Landesamt für Gesundheit und Lebensmittelsicherheit LGL") +
    theme(
      axis.text=element_text(size = rel(1.1)),
      axis.title=element_text(size = rel(1.3)))
  list(segmented_model = segmented, coef=mod_res, breakpoints = bp_dates_res, plot = bp_plot)
}


