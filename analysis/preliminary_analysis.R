#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(nlme)

# import data
# define relative light
benthic = read_csv("data/benthic_grad.csv") 
pelagic = read_csv("data/pelagic_grad.csv") 
profiles = read_csv("data/light_profile.csv")
shading = read_csv("data/shading.csv")

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
#========== Benthic: Plot
#==========

benthic %>%
  mutate(light = ifelse(site=="st33", 200*relative_light, 100*relative_light)) %>%
  ggplot(aes(light, do_flux, color=interaction(site, sampledate)))+
  geom_hline(yintercept = 0, size = 0.5, alpha = 0.5)+
  geom_point(size=3.5)+
  scale_color_manual("", values=c("dodgerblue","firebrick","gray"))+
  scale_y_continuous("Net Ecosystem Production")+
  theme_base




#==========
#========== Benthic: Fit model
#==========

# define gradient for P-I curve
mod = deriv(~beta*(1.08^2)*tanh((alpha/1000)*light/(beta*(1.08^2))) - rho*(1.11^2),
            c("beta", "alpha", "rho"), 
            function(beta, alpha, rho, light){})

# fit model
# allow paraemters to differ between sites
# don't yet have proper light data, so approximate
m = nlme(
  model = do_flux ~ mod(beta, alpha, rho, light),
  fixed = c(beta ~ site, alpha ~ site, rho ~ site),
  random = rho ~ 1|dummy, 
  data = benthic %>%
    mutate(light = ifelse(site=="st33", 
                          200*relative_light, 
                          100*relative_light),
           dummy = 1,
           site = ifelse(site=="st33", 1, 0)) %>%
    filter(
      !(light > 50 & do_flux < 0)
    ),
  start = c(0.1, 0, 1, 0, 0.1, 0)
)

# examine values
summary(m)
anova(m)

# model predictions
nd = crossing(site = c("st33", "reyk"),
              light = 0:200,
              dummy = 1)
nd$do_flux = predict(m, newdata=nd)

# plot
benthic %>%
  mutate(light = ifelse(site=="st33", 200*relative_light, 100*relative_light)) %>%
  ggplot(aes(light, do_flux, color=site))+
  geom_hline(yintercept = 0, size = 0.5, alpha = 0.5)+
  geom_line(data = nd, size=0.75)+
  geom_point(size=3.5)+
  scale_color_manual(values=c("dodgerblue", "firebrick"))+
  scale_y_continuous("Net Ecosystem Production")+
  theme_base






#==========
#========== Plot Data
#==========

pelagic %>% 
  filter(site=="st33") %>%
  ggplot(aes(relative_light, do_flux ))+
  geom_hline(yintercept = 0, size = 0.5, alpha = 0.5)+
  geom_point(size=3.5)+
  scale_y_continuous("Net Ecosystem Production")+
  theme_base

pelagic %>% 
  filter(site=="reyk", sampledate=="2018-06-28") %>%
  ggplot(aes(relative_light, do_flux ))+
  geom_hline(yintercept = 0, size = 0.5, alpha = 0.5)+
  geom_point(size=3.5)+
  scale_y_continuous("Net Ecosystem Production")+
  theme_base

pelagic %>% 
  filter(site=="reyk", sampledate=="2018-07-17") %>%
  ggplot(aes(relative_light, do_flux ))+
  geom_hline(yintercept = 0, size = 0.5, alpha = 0.5)+
  geom_point(size=3.5)+
  scale_y_continuous("Net Ecosystem Production")+
  theme_base

