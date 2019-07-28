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

# set theme
theme_set(theme_bw() %+replace% 
            theme(panel.grid = element_blank(),
                  strip.background = element_blank(),
                  legend.margin = margin(0,0,0,0),
                  strip.text = element_text(size=10),
                  legend.text = element_text(size=10),
                  axis.text=element_text(size=10, color="black"),
                  axis.title.y=element_text(angle = 90 ,margin=margin(0,15,0,0)),
                  axis.title.x=element_text(margin=margin(15,0,0,0))))





#==========
#========== Analysis: Benthic
#==========

# plot
benthic_clean %>%
  filter(!is.na(par)) %>%
  mutate(site = factor(site), site = droplevels(site)) %>%
  ggplot(aes(par, do_flux, color = factor(sampledate)))+
  facet_wrap(~site, nrow = 2)+
  geom_hline(yintercept = 0, size = 0.5, alpha = 0.5)+
  geom_point(size=3.5, alpha = 0.6)+
  scale_y_continuous("Net Ecosystem Production")+
  scale_color_manual("",values=c("dodgerblue3","dodgerblue3","firebrick3","firebrick3",
                                 "gray40","gray40","magenta4","forestgreen"))+
  stat_smooth(method = "gam", formula = y~s(x, k =3), se = F)

#==========
#========== Analysis: Pelagic
#==========

# plot
pelagic_clean %>%
  filter(!is.na(par), sampledate != "2018-08-15") %>%
  mutate(site = factor(site), site = droplevels(site)) %>%
  ggplot(aes(par, do_flux, color = factor(sampledate)))+
  facet_wrap(~site, nrow = 2)+
  geom_hline(yintercept = 0, size = 0.5, alpha = 0.5)+
  geom_point(size=3.5, alpha = 0.7)+
  scale_y_continuous("Net Ecosystem Production")+
  scale_color_manual("",values=c("dodgerblue3","dodgerblue3","firebrick3","firebrick3",
                                 "gray40","gray40","orange","magenta4","forestgreen"))+
  stat_smooth(method = "gam", formula = y~s(x, k =3), se = F)

