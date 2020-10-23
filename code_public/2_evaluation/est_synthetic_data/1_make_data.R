library(tidyverse)
# Load epidemic curve (number of individuals with disease onset per day)
ep_curve = read_tsv("../../data_public/evaluation/epidemic_curve_synthetic_data.csv")

# Create person-specific dataset
dat = tibble(disease_start = rep(ep_curve$date,  times = ep_curve$n_dis))

# Sample reporting date
# Create day-specific reporting delay distribution based on discrete time hazard model
T = nrow(ep_curve)
now = max(ep_curve$date) 
ddChangepoint = c(lubridate::ymd("2020-03-01", 
                                 "2020-03-15",
                                 "2020-03-25",
                                 "2020-04-15", 
                                 "2020-04-30"))
D = 21
t02s = seq(now-T+1, to = now, by = 1)
W <- array(NA, dim = c(T, length(ddChangepoint), D + 1), 
           dimnames = list(as.character(t02s),
                           c(as.character(ddChangepoint)), 
                           paste("delay", 0:D, sep = "")))
for (t in 0:T-1) {
  for (i in 1:length(ddChangepoint)) {
    W[t + 1, i, ] <- pmax((t02s[t + 1] + 0:D) - ddChangepoint[i] + 1, 0)
  }
}  

m = NULL
maxDelay = D
p = matrix(NA, nrow = T, D+1)
haz = matrix(NA, nrow = T, D+1)

gamma = c(-5.5,-4.5,-4.5,-4,-3.5,-1,-1.5,-1.5,-1.5,-1.5,-1,-1.5,-1,-1,-0.5,-0.5,-1.5,-1.5,-1.5,-.5,-.5,0)

eta = c(-.02, .045, 0.02, -.02, -.045)

for (t in max(1,T-m):T) {
  #Time dependent delay distribution
  p[t,1] <- plogis(gamma[1] + eta %*% W[t,,1])
  for (d in 1:(maxDelay-1)) {
    haz[t,d+1] <- plogis(gamma[d+1] + eta %*% W[t,,d+1])
    p[t,d+1] <- (1-sum(p[t,1:d]))*haz[t,d+1]
  }
  p[t,maxDelay+1] <- (1-sum(p[t,1:maxDelay]))*1 #since haz[maxDelay+1]=1
}


p_tib = as.data.frame(p) %>% as_tibble
colnames(p_tib) = paste0("D",0:21)
p_tib = p_tib %>% mutate(disease_start = t02s) 
# Join person-specific data and defined delay distribution
samp_rep_dat = full_join(dat, p_tib, by=c("disease_start"))

# Sample reporting delay per person
set.seed(122)
samp_delay = apply(samp_rep_dat, MARGIN = 1, function(x) sample(0:21, size = 1, prob = x[paste0("D",0:21)]))
# Derive reporting date per person
dat = dat %>% mutate(rep_delay = samp_delay,
                     rep_date = disease_start + rep_delay)
# Plot number of disease onsets and number of reported cases per day
dat %>% group_by(date = disease_start) %>% 
  summarise(n_dis = n()) %>% 
  full_join(dat %>% group_by(date=rep_date) %>% summarise(n_rep = n())) %>%
  mutate(n_dis = replace_na(n_dis, 0),
         n_rep = replace_na(n_rep, 0)) %>%
  ggplot() + geom_line(aes(date, n_dis), col = "green") +
  geom_line(aes(date, n_rep), col = "orange") 
# Create final dataset
dat_mod = dat %>% mutate(rep_date_weekday = weekdays(rep_date)) %>%
  dplyr::select(rep_date, rep_date_weekday, disease_start)
# Save dataset
save(dat_mod, file = "../../data_public/evaluation/dat_mod_synthetic.RData")  
