#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(lubridate)
library(nlme)

# import data
# define relative light
benthic_clean = read_csv("data/clean_data/benthic_clean.csv") 
pelagic_clean = read_csv("data/clean_data/pelagic_clean.csv") 

# define base theme
theme_base = theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.margin = margin(0,0,0,0),
        text = element_text(size=12),
        strip.text = element_text(size=10),
        legend.text = element_text(size=10),
        axis.text=element_text(size=10, color="black"),
        axis.title.y=element_text(margin=margin(0,15,0,0)),
        axis.title.x=element_text(margin=margin(15,0,0,0)))





#==========
#========== Analysis: Benthic
#==========

# plot
benthic_clean %>%
  mutate(month = month(sampledate, label=T)) %>%
  ggplot(aes(par, do_flux, color = month))+
  facet_wrap(~site, nrow = 2)+
  geom_hline(yintercept = 0, size = 0.5, alpha = 0.5)+
  geom_point(size=3.5, alpha = 0.7)+
  scale_y_continuous("Net Ecosystem Production")+
  scale_color_manual("",values=c("dodgerblue3","firebrick3","black"))+
  theme_base

# define gradient for P-I curve
mod = deriv(~beta*(1.08^2)*tanh((alpha/1000)*light/(beta*(1.08^2))) - rho*(1.11^2),
            c("beta", "alpha", "rho"), 
            function(beta, alpha, rho, light){})

# fit model
# allow paraemters to differ between sites and by date
# code dates as contrasts (jul and aug vs. june)
# cannot include interaction between site and august becuase there is only data for Reyk
m_full = nlme(
  model = do_flux ~ mod(beta, alpha, rho, par),
  fixed = c(beta ~ site + jul + site:jul + aug, 
            alpha ~ site, 
            rho ~ site + jul + site:jul + aug),
  random = rho ~ 1|dummy, 
  data = benthic_clean %>%
    mutate(dummy = 1,
           site = ifelse(site=="st33", 1, 0),
           month = month(sampledate, label=T),
           jul = ifelse(month == "Jul", 1, 0),
           aug = ifelse(month == "Aug", 1, 0)),
  start = c(0.1, 0, 0, 0, 0, 1, 0, 0.1, 0, 0, 0, 0)
)

# remove time effects
m_slim_1 = nlme(
  model = do_flux ~ mod(beta, alpha, rho, par),
  fixed = c(beta ~ site, 
            alpha ~ site, 
            rho ~ site),
  random = rho ~ 1|dummy, 
  data = benthic_clean %>%
    mutate(dummy = 1,
           site = ifelse(site=="st33", 1, 0),
           month = month(sampledate, label=T),
           jul = ifelse(month == "Jul", 1, 0),
           aug = ifelse(month == "Aug", 1, 0)),
  start = c(0.1, 0, 1, 0, 0.1, 0)
)

# likelihood ratio test
anova(m_full, m_slim_1)

# data frame of modeled values
benthic_pred = benthic_clean %>% 
  split(.$sampledate) %>% 
  lapply(function(x){data_frame(site = unique(x$site),
                                sampledate = unique(x$sampledate),
                                par = seq(0, max(x$par), 1))}) %>%
  bind_rows() %>%
  mutate(month = month(sampledate, label=T),
         jul = ifelse(month == "Jul", 1, 0),
         aug = ifelse(month == "Aug", 1, 0))
  
benthic_pred$do_flux = 
  predict(m_full, newdata = 
            benthic_pred %>% 
            mutate(dummy = 1,
                   site = ifelse(site=="st33", 1, 0)  
                   )) 

# plot
benthic_clean %>%
  mutate(month = month(sampledate, label=T)) %>%
  ggplot(aes(par, do_flux, color = month))+
  facet_wrap(~site, nrow = 2)+
  geom_hline(yintercept = 0, size = 0.5, alpha = 0.5)+
  geom_point(size=3.5, alpha = 0.7)+
  geom_line(data = benthic_pred, size = 1)+
  scale_y_continuous("Net Ecosystem Production")+
  scale_color_manual("",values=c("dodgerblue3","firebrick3","gray40"))+
  theme_base





#==========
#========== Analysis: Pelagic
#==========

# plot
pelagic_clean %>%
  mutate(month = month(sampledate, label=T)) %>%
  ggplot(aes(par, do_flux, color = month))+
  facet_wrap(~site, nrow = 2)+
  geom_hline(yintercept = 0, size = 0.5, alpha = 0.5)+
  geom_point(size=3.5, alpha = 0.7)+
  scale_y_continuous("Net Ecosystem Production")+
  scale_color_manual("",values=c("dodgerblue3","firebrick3","gray40"))+
  theme_base

# define gradient for P-I curve
mod = deriv(~beta*light/ 
              ((beta/(alpha*opt^2))*light^2 +  
                 (1 - 2*(beta/(alpha*opt)))*light +  
                 beta/alpha) - rho,
            c("beta", "alpha","opt" , "rho"), 
            function(beta, alpha, opt, rho, light){})

m_full = nlme(
  model = do_flux ~ mod(beta, alpha, opt, rho, par),
  fixed = c(beta ~ site, 
            alpha ~ site, 
            opt ~ site,
            rho ~ site),
  random = rho ~ 1|dummy, 
  data = pelagic_clean %>%
    filter(is.na(do_flux) == F,
           !(site=="st33" & par > 400 & do_flux > 0.02),
           !(par == 0 & do_flux > 0)) %>%
    mutate(dummy = 1,
           site = ifelse(site=="st33", 1, 0),
           month = month(sampledate, label=T),
           jul = ifelse(month == "Jul", 1, 0),
           aug = ifelse(month == "Aug", 1, 0)),
  start = c(0.07, 0, 0.0005, 0, 120, 0, 0.04, 0)
)

anova(m_full)

# data frame of modeled values
pelagic_pred = pelagic_clean %>% 
  split(.$sampledate) %>% 
  lapply(function(x){data_frame(site = unique(x$site),
                                sampledate = unique(x$sampledate),
                                par = seq(0, max(x$par), 1))}) %>%
  bind_rows() %>%
  mutate(month = month(sampledate, label=T),
         jul = ifelse(month == "Jul", 1, 0),
         aug = ifelse(month == "Aug", 1, 0))

pelagic_pred$do_flux = 
  predict(m_full, newdata = 
            pelagic_pred %>% 
            mutate(dummy = 1,
                   site = ifelse(site=="st33", 1, 0)  
            )) 

# plot
pelagic_clean %>%
  mutate(month = month(sampledate, label=T)) %>%
  ggplot(aes(par, do_flux, color = month))+
  facet_wrap(~site, nrow = 2)+
  geom_hline(yintercept = 0, size = 0.5, alpha = 0.5)+
  geom_point(size=3.5, alpha = 0.7)+
  geom_line(data = pelagic_pred, size = 1, color = "black")+
  scale_y_continuous("Net Ecosystem Production")+
  scale_color_manual("",values=c("dodgerblue3","firebrick3","gray40"))+
  theme_base



