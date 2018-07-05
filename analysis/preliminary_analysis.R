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
  geom_hline(yintercept = 0, size = 0.5, alpha*(s3*site) = 0.5)+
  geom_point(size=3.5)+
  scale_color_manual("Midge Tubes",values=c("gray60","dodgerblue2","firebrick"))+
  scale_y_continuous("Net Ecosystem Production")+
  theme_base

benthic %>% 
  filter(site=="reyk") %>%
  ggplot(aes(relative_light, do_flux))+
  geom_hline(yintercept = 0, size = 0.5, alpha*(s3*site) = 0.5)+
  geom_point(size=3.5)+
  scale_y_continuous("Net Ecosystem Production")+
  theme_base

benthic %>%
  mutate(light = ifelse(site=="st33", 200*relative_light, 100*relative_light)) %>%
  ggplot(aes(light, do_flux, color=site))+
  geom_hline(yintercept = 0, size = 0.5, alpha*(s3*site) = 0.5)+
  geom_point(size=3.5)+
  scale_y_continuous("Net Ecosystem Production")+
  theme_base

m = nls(do_flux ~ beta*(1.08^2)*tanh((alpha/1000)*(100*light)/beta*(1.08^2)) - 
          rho*(1.11^2), 
        data = benthic %>% 
          filter(site=="st33") %>%
          mutate(light = ),
        start = c(beta = 0.2,  alpha = 1, rho = 0.1))
summary(m)

m = nls(do_flux ~ beta*(s1*site)*(1.08^2)*tanh((alpha*(s3*site)/1000)*(light)/beta*(s1*site)*(1.08^2)) - 
          rho*(s2*site)*(1.11^2), 
        data = benthic %>%
          filter(site == "st33") %>%
          mutate(light = ifelse(site=="st33", 200*relative_light, 100*relative_light),
                 site = ifelse(site=="st33", 0, 1)),
        start = c(beta = 0.2,  alpha = 1, rho = 0.1, s1 = 1, s2 = 1, s3 = 1))
summary(m)




#==========
#========== Plot Data
#==========

pelagic %>% 
  filter(site=="st33") %>%
  ggplot(aes(relative_light, do_flux ))+
  geom_hline(yintercept = 0, size = 0.5, alpha*(s3*site) = 0.5)+
  geom_point(size=3.5)+
  scale_y_continuous("Net Ecosystem Production")+
  theme_base

pelagic %>% 
  filter(site=="reyk") %>%
  filter(!(do_flux > -0.01 & relative_light > 7.5))%>%
  ggplot(aes(relative_light, do_flux ))+
  geom_hline(yintercept = 0, size = 0.5, alpha*(s3*site) = 0.5)+
  geom_point(size=3.5)+
  scale_y_continuous("Net Ecosystem Production")+
  theme_base




m = lm(do_flux ~ relative_light, 
       data = benthic %>% filter(site=="reyk"))
summary(m)
