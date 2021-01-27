#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(rstan)
# library(loo)
source("analysis/fit_model/stan_utility.R")

# stan settings
rstan_options(auto_write = TRUE)
# options(mc.cores = parallel::detectCores()-2)

cl <- parallel::makeCluster(2, setup_strategy = "sequential")

# read data
type <- "benthic"
prep_data <- read_csv(paste0("analysis/fit_model/input/",type,"/prep_data.csv"))
data_list <- read_rdump(paste0("analysis/fit_model/input/",type,"/data_list.R"))

# set reference temperature
data_list$temp_ref <- 12






#==========
#========== Fit model
#==========

# process data
# data_list$event_site <- rep(1, length(data_list$event_site))
data_list$event_site <- unique(data_list$event)
data_list$n_sites <- length(data_list$event_site)

# specify model list
model <- "pi_curve.stan"

# model path
model_path <- paste0("analysis/fit_model/stan/",model)

# MCMC specifications
chains <- 4
iter <- 10
adapt_delta <- 0.95
max_treedepth <- 15

# fit model
fit <- stan(file = model_path, data = data_list, seed=1, chains = chains, iter = iter, 
            control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth))

# summary of fit
fit_summary <- summary(fit, probs=c(0.16, 0.5, 0.84))$summary %>% 
{as_tibble(.) %>%
    mutate(var = rownames(summary(fit)$summary))}

# check Rhat & n_eff
fit_summary %>% filter(Rhat > 1.01) %>% select(Rhat, n_eff, var) %>% arrange(-Rhat)
fit_summary %>% filter(n_eff < 0.5*(chains*iter/2)) %>% select(Rhat, n_eff, var) %>% arrange(n_eff) %>%
  mutate(eff_frac = n_eff/(chains*iter/2))

# additional diagnostics
check_div(fit)
check_treedepth(fit,max_treedepth)
check_energy(fit)




#==========
#========== Examine Chains
#==========

fixed_pars_fn <- function(x){
  if(x == "pi_curve.stan") return(c("g_b","g_r","sig_b0","sig_a","sig_d","sig_r","sig_obs","lp__"))
  if(x == "pi_curve_tanh.stan") return(c("g_b","g_r","sig_b0","sig_a","sig_r","sig_obs","lp__"))
  if(x == "pi_presp.stan") return(c("g_b","g_r","sig_b0","sig_a","sig_d","sig_r","sig_w","sig_obs","lp__"))
  }

# fixed parameters by step
fixed_pars <- rstan::extract(fit, pars=fixed_pars_fn(model)) %>%
  lapply(as_tibble) %>%
  bind_cols() %>%
  set_names(fixed_pars_fn(model)) %>%
  mutate(chain = rep(1:chains, each = iter/2), step = rep(c(1:(iter/2)), chains))

# examine chains for parameters
fixed_pars %>%
  gather(par, value, -chain, -step) %>%
  filter(par != "lp__") %>%
  ggplot(aes(step, value, color=factor(chain)))+
  facet_wrap(~par, scales="free_y")+
  geom_line(alpha=0.5)+
  theme_bw()

# posterior densities
fixed_pars %>%
  gather(par, value, -chain, -step) %>%
  filter(par != "lp__") %>%
  ggplot(aes(value))+
  facet_wrap(~par, scales="free")+
  stat_density(alpha=0.5, geom = "line")+
  theme_bw()





#==========
#========== Export
#==========

# event_names
event_names <- unique(prep_data$event_name)

# clean variable names in summary
fit_clean <- fit_summary %>%
  rename(lower16 = `16%`, middle = `50%`, upper84 = `84%`)  %>%
  mutate(name = strsplit(var, "\\[|\\]|,") %>% map_chr(~.x[1]),
         index = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         event_name = ifelse(name %in% c("log_beta0","log_alpha","log_delta","log_rho","log_omega",
                                  "beta0","alpha","delta","rho","omega","b0","a","d","r","w"), 
                             event_names[index], NA),
         site = strsplit(event_name, " ") %>% map_chr(~.x[1]),
         sampledate = lubridate::ymd(strsplit(event_name, " ") %>% map_chr(~.x[2]))) %>%
  select(name, index, event_name, site, sampledate, middle, lower16, upper84)








