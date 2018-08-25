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

# create data list
data_list = benthic_clean %>%
  arrange(site, sampledate, id) %>%
  # define events
  left_join(benthic_clean %>%
              expand(nesting(sampledate, site)) %>%
              arrange(site, sampledate) %>%
              mutate(event = 1:length(site))) %>%
  # site as numeric
  mutate(site = ifelse(site == "st33", 2, 1)) %>%
  # create list
  {list(N = nrow(.),
      Nsites = length(unique(.$site)),
      Nevents = length(unique(.$event)),
      site = .$site,
      do_flux = .$do_flux,
      par = .$par,
      temp = .$temp,
      # temp ref from sonde data
      temp_ref = 12)}

# create priors list
# weakly informative priors informed by analysis of sonde data
list(
  prior_beta_site = 0.35,
  prior_alpha_site = 0.005,
  prior_opt_site = 500,
  prior_rho_site = 0.12,
)


