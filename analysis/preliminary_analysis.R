#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(readxl)

# import data
# define relative light
benthic = read_excel("data/field_incubations_5jul18.xlsx", sheet = "BenthicGradient",
                     na = "NA") %>%
  mutate(light_trt  = ifelse(light_trt == 10, 9, light_trt),
         relative_light = max(light_trt) - light_trt)
pelagic = read_excel("data/field_incubations_5jul18.xlsx", sheet = "Pelagic",
                     na = "NA") %>%
  mutate(light_trt  = ifelse(light_trt == 10, 9, light_trt),
         relative_light = max(light_trt) - light_trt)

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
#========== Plot Data
#==========

benthic %>% 
  filter(site=="st33") %>%
  ggplot(aes(relative_light, do_flux))+
  geom_hline(yintercept = 0, size = 0.5, alpha = 0.5)+
  geom_point(size=3.5)+
  scale_color_manual("Midge Tubes",values=c("gray60","dodgerblue2","firebrick"))+
  scale_y_continuous("Net Ecosystem Production")+
  theme_base

benthic %>% 
  filter(site=="reyk") %>%
  ggplot(aes(relative_light, do_flux))+
  geom_hline(yintercept = 0, size = 0.5, alpha = 0.5)+
  geom_point(size=3.5)+
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
  filter(site=="reyk") %>%
  ggplot(aes(relative_light, do_flux ))+
  geom_hline(yintercept = 0, size = 0.5, alpha = 0.5)+
  geom_point(size=3.5)+
  scale_y_continuous("Net Ecosystem Production")+
  theme_base


