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

benthic_clean %>%
  ggplot(aes(par, do_flux, color=interaction(site, sampledate)))+
  geom_hline(yintercept = 0, size = 0.5, alpha = 0.5)+
  geom_point(size=3.5)+
  scale_color_manual("", values=c("dodgerblue","firebrick","gray","orange","magenta"))+
  scale_y_continuous("Net Ecosystem Production")+
  theme_base

# define gradient for P-I curve
mod = deriv(~beta*(1.08^2)*tanh((alpha/1000)*light/(beta*(1.08^2))) - rho*(1.11^2),
            c("beta", "alpha", "rho"), 
            function(beta, alpha, rho, light){})

# fit model
# allow paraemters to differ between sites
# don't yet have proper light data, so approximate
m = nlme(
  model = do_flux ~ mod(beta, alpha, rho, par),
  fixed = c(beta ~ site*time, alpha ~ site*time, rho ~ site*time),
  random = rho ~ 1|dummy, 
  data = benthic_clean %>%
    mutate(dummy = 1,
           site = ifelse(site=="st33", 1, 0),
           time = as.numeric(sampledate) - min(as.numeric(sampledate))),
  start = c(0.1, 0, 0, 0, 1, 0, 0, 0, 0.1, 0, 0, 0)
)



# examine values
summary(m)
anova(m)

# data frame of modeled values
benthic_pred = benthic_clean %>%
  expand(nesting(sampledate, site), par = seq(min(par), max(par), 0.1)) 
benthic_pred$do_flux = 
  predict(m, newdata = 
            benthic_pred %>% 
            mutate(dummy = 1,
                   site = ifelse(site=="st33", 1, 0),
                   time = as.numeric(sampledate) - min(as.numeric(sampledate))
                   )) 

# plot
benthic_clean %>%
  ggplot(aes(par, do_flux, color=interaction(site, sampledate)))+
  # facet_wrap(site ~ sampledate, scales="free_x")+
  geom_hline(yintercept = 0, size = 0.5, alpha = 0.5)+
  geom_point(size=3.5)+
  geom_line(data = benthic_pred, size = 1)+
  scale_color_manual("", values=c("dodgerblue","firebrick","gray","orange","magenta"))+
  scale_y_continuous("Net Ecosystem Production")+
  theme_base

summary(m)$tTable



#==========
#========== Analysis: Pelagic
#==========

pelagic_clean %>%
  ggplot(aes(par, do_flux, color=interaction(site, sampledate)))+
  facet_wrap(site ~ sampledate, scales="free_x")+
  geom_hline(yintercept = 0, size = 0.5, alpha = 0.5)+
  geom_point(size=3.5)+
  scale_color_manual("", values=c("dodgerblue","firebrick","gray","orange","magenta"))+
  scale_y_continuous("Net Ecosystem Production")+
  theme_base

# define gradient for P-I curve
mod = deriv(~beta*(1.08^2)*tanh((alpha/1000)*light/(beta*(1.08^2))) - rho*(1.11^2),
            c("beta", "alpha", "rho"), 
            function(beta, alpha, rho, light){})

# fit model
# allow paraemters to differ between sites
# don't yet have proper light data, so approximate
m = nlme(
  model = do_flux ~ mod(beta, alpha, rho, par),
  fixed = c(beta ~ site, alpha ~ site, rho ~ site),
  random = rho ~ 1|dummy, 
  data = pelagic_clean %>%
    filter(is.na(do_flux) == F,
           sampledate > "2018-06-28") %>%
    mutate(dummy = 1,
           site = ifelse(site=="st33", 1, 0)) ,
  start = c(0.01, 0, 1, 0, 0.01, 0)
)

summary(m)
anova(m)


