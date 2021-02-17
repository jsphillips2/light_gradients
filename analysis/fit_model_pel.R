#=========================================================================================
#========== Preliminaries
#=========================================================================================

# load packages
library(tidyverse)
library(rstan)

# stan settings
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-2)

# import data
dd <- read_csv("data/clean_data/pelagic_clean.csv")  %>%
  filter(site != "e5",
         sampledate != "2018-08-15") %>%
  mutate(date = factor(sampledate),
         event = factor(paste0(date, site))) %>%
  select(event, date, site,  par, temp, do_flux) %>%
  na.omit()

# data frame for summary
dd_sum <- dd %>%
  split(.$event) %>%
  lapply(function(x) {
    x %>% 
      tidyr::expand(par = seq(min(par), max(par), length.out = 500)) %>%
      mutate(date = unique(x$date),
             site = unique(x$site),
             event = unique(x$event))
  }) %>%
  bind_rows()

#=========================================================================================





#=========================================================================================
#========== Package data for fitting
#=========================================================================================

# define variables
n_obs <- nrow(dd)
x <- cbind(as.numeric(dd$event),
           as.numeric(dd$event),
           as.numeric(dd$event))
par <- dd$par
do_flux <- dd$do_flux
n_sum <- nrow(dd_sum)
x_sum <- cbind(as.numeric(dd_sum$event),
               as.numeric(dd_sum$event),
               as.numeric(dd_sum$event))
par_sum <- dd_sum$par

# package variables
data_list <- list(
  n_obs = n_obs,
  x = x,
  par = par,
  do_flux = do_flux,
  n_sum = n_sum,
  x_sum = x_sum,
  par_sum = par_sum
)

#=========================================================================================





#=========================================================================================
#========== Fit model
#=========================================================================================

# MCMC specifications
chains <- 4
iter <- 2000
adapt_delta <- 0.95
max_treedepth <- 10

# # fit model
# fit <- stan(file = "analysis/stan/pi_curve.stan", data = data_list, seed=2e3,
#             chains = chains, iter = iter,
#             control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth))
# 
# 
# # summary of fit
# fit_summary <- summary(fit, probs=c(0.16, 0.5, 0.84))$summary %>%
#   {as_tibble(.) %>%
#       mutate(var = rownames(summary(fit)$summary))}
#
# # export
# write_rds(list(dd = dd,
#                dd_sum = dd_sum,
#                data_list = data_list,
#                fit = fit,
#                fit_summary = fit_summary),
#            paste0("analysis/output/pelagic.rds"))

#=========================================================================================