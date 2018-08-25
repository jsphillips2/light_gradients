#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(rstan)

# stan settings
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-1)

# read data
benthic_clean = read_csv("data/clean_data/benthic_clean.csv") 





#==========
#========== Prpare Data
#==========

# create event x site data frame
events = benthic_clean %>%
  expand(nesting(sampledate, site)) %>%
  arrange(site, sampledate) %>%
  mutate(event = 1:length(site))

# create list of main data
main_data = benthic_clean %>%
  arrange(site, sampledate, id) %>%
  # define events
  left_join(events) %>%
  # site as numeric
  mutate(site = ifelse(site == "st33", 2, 1)) %>%
  select(site, event, do_flux, par, temp) %>%
  as.list() 
# add temp_ref, based on sonde data
main_data$temp_ref = 12

# create list of indices
index_data = main_data %>%
{list(N = length(.$do_flux),
      Nsites = length(unique(.$site)),
      Nevents = length(unique(.$event)),
      event_sites = {events %>% mutate(site = ifelse(site == "st33", 2, 1))
        }$site)}

# priors
priors = list(prior_beta_site = c(0.35, 0.35),
              prior_alpha_site = c(1, 1),
              # prior_opt_site = c(1000, 10),
              prior_rho_site = c(0.12, 0.12),
              prior_sig_beta = c(0.01, 0.03),
              prior_sig_alpha = c(0.1, 0.3),
              # prior_sig_opt = c(1, 1),
              prior_sig_rho = c(0.01, 0.03),
              prior_gamma_1 = c(1.9, 0.005),
              prior_gamma_2 = c(1.12, 0.005),
              prior_sig_obs = c(0.02, 0.02))

# combine into full data list 
data_list = index_data %>% 
  append(main_data) %>%
  append(priors)
data_list$site = NULL

# function for initial values
init_fn = function(x){x %>% 
    {list(
      beta_site = 
        runif(.$Nsites, 
              min = .$prior_beta_site[1]*0.5, 
              max = .$prior_beta_site[1]*1.5),
      alpha_site = 
        runif(.$Nsites, 
              min = .$prior_alpha_site[1]*0.5, 
              max = .$prior_alpha_site[1]*1.5),
      # opt_site = 
      #   runif(.$Nsites, 
      #         min =.$prior_opt_site[1]*0.5, 
      #         max = .$prior_opt_site[1]*1.5),
      rho_site = 
        runif(.$Nsites,
              min = .$prior_rho_site[1]*0.5, 
              max = .$prior_rho_site[1]*1.5),
      sig_beta = 
        runif(1, 
              min = .$prior_sig_beta[1]*0.5, 
              max = .$prior_sig_beta[1]*1.5),
      sig_alpha = 
        runif(1, 
              min = .$prior_sig_alpha[1]*0.5, 
              max = .$prior_sig_alpha[1]*1.5),
      # sig_opt = 
      #   runif(1, 
      #         min = .$prior_sig_opt[1]*0.5, 
      #         max = .$prior_sig_opt[1]*1.5),
      sig_rho = 
        runif(1, 
              min = .$prior_sig_rho[1]*0.5, 
              max = .$prior_sig_rho[1]*1.5),
      gamma_1 = 
        runif(1, 
              min = max(1, .$prior_gamma_1[1]*0.9), 
              max = .$prior_gamma_1[1]*1.1),
      gamma_2 = 
        runif(1, 
              min = max(1, .$prior_gamma_2[1]*0.9), 
              max = .$prior_gamma_2[1]*1.5),
      sig_obs = 
        runif(1, 
              min = .$prior_sig_obs[1]*0.5, 
              max = .$prior_sig_obs[1]*1.5)
    )}}
inits = function(){init_fn(data_list)}

#==========
#========== Fit model
#==========

# initial specifications
model = "photo_irrad"
model_path = paste0("analysis/fit_model/",model,".stan")
chains = 1
iter = 3000

# fit model
fit = stan(file = model_path, data = data_list, seed=1, chains = chains,
           init = inits, iter = iter, control = list(adapt_delta = 0.8))

# summary of fit
fit_summary = summary(fit, probs=c(0.16, 0.5, 0.84))$summary %>% 
{as_data_frame(.) %>%
    mutate(var = rownames(summary(fit)$summary))}

# check Rhat
fit_summary %>% filter(Rhat > 1.05) %>% select(Rhat, n_eff, var)




