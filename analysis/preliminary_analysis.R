#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(lubridate)
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
#========== Process data
#==========

# calculate effect of shading
shading_trt  = shading %>%
  mutate(light_frac = ifelse(par_in > par_out, 1, par_in/par_out)) %>%
  group_by(light_trt) %>%
  summarize(light_frac = mean(light_frac))

# benthic light
benthic_light = benthic %>%
  group_by(sampledate, site) %>%
  summarize(time_start = min(time_start),
            time_end = max(time_end),
            depth = unique(inc_depth)) %>%
  left_join(profiles %>%
              filter(sampledate != "2018-06-27", site != "grim") %>%
              mutate(site = ifelse(site=="btl", "reyk", site))) %>%
  mutate(within = time > time_start - 60*30 & time < time_end + 60*30) %>%
  group_by(site, sampledate) %>%
  summarize(par = mean(par))

# pelagic light
pelagic_light = pelagic %>%
  filter(is.na(time_start)==F, is.na(time_end)==F, site != "grim") %>%
  group_by(sampledate, site) %>%
  summarize(time_start = min(time_start),
            time_end = max(time_end),
            depth = unique(inc_depth)) %>%
  left_join(profiles %>%
              filter(sampledate != "2018-06-27", site != "grim") %>%
              mutate(site = ifelse(site=="btl", "reyk", site))) %>%
  mutate(within = time > time_start - 60*30 & time < time_end + 60*30) %>%
  group_by(site, sampledate) %>%
  summarize(par = mean(par))

# add light data to benthic
benthic_full = benthic %>%
  left_join(benthic_light) %>%
  left_join(shading_trt) %>%
  mutate(par = light_frac*par)

# add light data to benthic
pelagic_full = pelagic %>%
  filter(site != "grim") %>%
  left_join(pelagic_light) %>%
  left_join(shading_trt) %>%
  mutate(par = light_frac*par) 





#==========
#========== Analysis: Benthic
#==========

benthic_full %>%
  ggplot(aes(par, do_flux, color=interaction(site, sampledate)))+
  geom_hline(yintercept = 0, size = 0.5, alpha = 0.5)+
  geom_point(size=3.5)+
  scale_color_manual("", values=c("dodgerblue","firebrick","gray","orange"))+
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
  fixed = c(beta ~ site*sampledate, alpha ~ site*sampledate, rho ~ site*sampledate),
  random = rho ~ 1|dummy, 
  data = benthic_full %>%
    mutate(dummy = 1,
           site = ifelse(site=="st33", 1, 0)) ,
  start = c(0.1, 0, 0, 0, 1, 0, 0, 0, 0.1, 0, 0, 0)
)



# examine values
summary(m)
anova(m)

# data frame of modeled values
benthic_pred = benthic_full %>%
  expand(nesting(sampledate, site), par = seq(min(par), max(par), 0.1)) 
benthic_pred$do_flux = predict(m, newdata = benthic_pred %>%
                                 mutate(dummy = 1,
                                        site = ifelse(site=="st33", 1, 0)))

# plot

benthic_full %>%
  ggplot(aes(par, do_flux, color=interaction(site, sampledate)))+
  # facet_wrap(site ~ sampledate, scales="free_x")+
  geom_hline(yintercept = 0, size = 0.5, alpha = 0.5)+
  geom_point(size=3.5)+
  geom_line(data = benthic_pred, size = 1)+
  scale_color_manual("", values=c("dodgerblue","firebrick","gray","orange"))+
  scale_y_continuous("Net Ecosystem Production")+
  theme_base





#==========
#========== Analysis: Pelagic
#==========

pelagic_full %>%
  ggplot(aes(par, do_flux, color=interaction(site, sampledate)))+
  facet_wrap(site ~ sampledate, scales="free_x")+
  geom_hline(yintercept = 0, size = 0.5, alpha = 0.5)+
  geom_point(size=3.5)+
  scale_color_manual("", values=c("dodgerblue","firebrick","gray","orange"))+
  scale_y_continuous("Net Ecosystem Production")+
  theme_base

