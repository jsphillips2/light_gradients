#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(rstan)

# import data
benthic_clean <- read_csv("data/clean_data/benthic_clean.csv") 
pelagic_clean <- read_csv("data/clean_data/pelagic_clean.csv") %>%
  filter(site != "e5") %>%
  filter(!is.na(par), sampledate != "2018-08-15")





#==========
#========== Package for Stan
#==========

# function for preparing data
prep_fn <- function(data_, path_){
  
  # prepare data
  data = data_ %>%
    mutate(event_name = as.factor(paste(site, sampledate)),
           site_name = as.factor(site),
           event = as.numeric(event_name),
           site = as.numeric(site_name)) %>%
    select(site_name, site, event_name, event, do_flux, par, temp) %>%
    na.omit() %>%
    arrange(event)
  
  # export
  return(data)
  
}

# define function for packaging
pckg_fn <- function(data_, path_){
  
  # extract variables
  do_flux = data_$do_flux
  temp = data_$temp
  par = data_$par
  event = data_$event
  event_site = {data_ %>% group_by(event) %>% summarize(site = unique(site))}$site
  n_obs = length(do_flux)
  n_events = length(event_site)
  n_sites = length(unique(event_site))
  
  # export
  stan_rdump(c("n_obs","n_events","n_sites","do_flux","temp","par","event","event_site"), file=path_)
}

# export benthic
# write_csv(prep_fn(benthic_clean), "analysis/fit_model/input/benthic/prep_data.csv")
# pckg_fn(prep_fn(benthic_clean), "analysis/fit_model/input/benthic/data_list.R")

# export pelagic
# write_csv(prep_fn(pelagic_clean), "analysis/fit_model/input/pelagic/prep_data.csv")
# pckg_fn(prep_fn(pelagic_clean), "analysis/fit_model/input/pelagic/data_list.R")


