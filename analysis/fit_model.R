#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(rstan)
# source("analysis/model_fn.R")

# stan settings
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-2)

# import data
dd <- read_csv("data/clean_data/benthic_clean.csv") %>%
  mutate(date = factor(sampledate),
         event = factor(paste0(date, site))) %>% 
  select(date, site, event, par, temp, do_flux) %>%
  na.omit()

# data frame for summary
dd_sum <- dd %>%
  split(.$event) %>%
  lapply(function(x) {
    x %>% 
      tidyr::expand(par = seq(min(par), max(par), length.out = 50)) %>%
      mutate(date = unique(x$date),
             site = unique(x$site),
             event = unique(x$event))
  }) %>%
  bind_rows()




#==========
#========== Package data for fitting
#==========

# define variables
temp_ref <- mean(dd$temp)
kb <- 8.62e-5
n_obs <- nrow(dd)
x <- cbind(as.numeric(dd$event),
           as.numeric(dd$event),
           as.numeric(dd$event))
tscale <- 1 / (kb * (dd$temp + 273.15)) - 1 / (kb * (temp_ref + 273.15)) 
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
  tscale = tscale,
  par = par,
  do_flux = do_flux,
  n_sum = n_sum,
  x_sum = x_sum,
  par_sum = par_sum
)





#==========
#========== Fit model
#==========

# MCMC specifications
chains <- 4
iter <- 2000
adapt_delta <- 0.9
max_treedepth <- 10

# # fit model
# fit <- stan(file = "analysis/pi_curve.stan", data = data_list, seed=2e3,
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
# write_rds(list(data_list = data_list,
#                fit = fit,
#                fit_summary = fit_summary),
#            paste0("analysis/output/benthic.rds"))



benthic_rds <- read_rds("analysis/output/benthic.rds")

fit <- benthic_rds$fit
fit_summary <- benthic_rds$fit_summary




nep_sum <- fit_summary %>%
  filter(str_detect(.$var, "nep_sum")) %>%
  mutate(name = strsplit(var, "\\[|\\]|,") %>% map_chr(~.x[1]),
         id = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2]))) %>%
  filter(name %in% c("nep_sum")) %>%
  rename(lo = `16%`,
         mi = `50%`,
         hi = `84%`) %>%
  select(id, lo, mi, hi) %>%
  unique() %>%
  left_join(dd_sum %>%
              mutate(id = row_number())) %>%
  select(-id) %>%
  unique()


dd %>%
  ggplot(aes(par, do_flux, color = date))+
  facet_wrap(~site, nrow = 2)+
  geom_ribbon(data = nep_sum,
              aes(x = par,
                  ymin = lo,
                  ymax = hi,
                  fill = date),
              linetype = 0,
              alpha = 0.2,
              inherit.aes = F)+
  geom_point()+
  geom_line(data = nep_sum,
            aes(y = mi))+
  theme_bw()+
  scale_color_viridis_d()+
  scale_fill_viridis_d()


pars <-  fit_summary %>%
  filter(var %in% c(paste0("b[",1:11,"]"),
                    paste0("a[",1:11,"]"),
                    paste0("r[",1:11,"]"))) %>%
  mutate(name = strsplit(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
         id = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2]))) %>%
  rename(lo = `16%`,
         mi = `50%`,
         hi = `84%`) %>%
  select(name, id, lo, mi, hi) %>%
  left_join(dd %>%
              select(date, site, event) %>%
              unique() %>%
               mutate(id = row_number())) %>%
  select(-id) %>%
  unique() 


pars %>%
  mutate(dummya = as.numeric(as.factor(site))) %>%
  group_by(site) %>%
  mutate(dummy = as.numeric(date) - mean(as.numeric(date))) %>%
  ungroup() %>%
  mutate(dummy = dummy / abs(max(dummy)) / 3,
         dummyb = dummya + dummy) %>%
  ggplot(aes(dummyb, mi, color = site))+
  facet_wrap(~name, scales = "free_y", nrow = 3)+
  geom_point(size = 2)+
  geom_errorbar(aes(ymin = lo,
                  ymax = hi),
                width = 0)+
  theme_bw()+
  scale_color_manual(values = c("dodgerblue",
                                "firebrick",
                                "black"))
  


 


fit_summary %>%
  filter(str_detect(.$var, "e_b")) 


fit_summary %>%
  filter(str_detect(.$var, "e_r")) 
