#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(lubridate)
library(nlme)

# import data
# define relative light
benthic = read_csv("data/raw_data/benthic_grad.csv") 
pelagic = read_csv("data/raw_data/pelagic_grad.csv") 
profiles = read_csv("data/raw_data/light_profile.csv") 
shading = read_csv("data/raw_data/shading.csv") 

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
#========== Process data: Light
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
              mutate(site = ifelse(site=="btl", "reyk", site),
                     site = ifelse(site == "kal", "st33", site))) %>%
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
              mutate(site = ifelse(site=="btl", "reyk", site),
                     site = ifelse(site == "kal", "st33", site))) %>%
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
#========== Process data: DO Flux
#==========

# benthic
bethic_clean = benthic_full %>%
  mutate(duration = as.numeric(time_end - time_start),
         # calculate change in do
         # multiply by 1000 to convert to mg m^-3
         # multiply by column depth in m's to convert to mg m^-2
         # divde by 1000 to convert to g m^-2
         # divide by duration to convert to g m^-2 h-1
         do_flux = 1000*(column_depth/100)*(do_end - do_start)/1000/duration,
         # mean temperature
         temp = (temp_start + temp_end)/2) %>%
  select(id, sampledate, site, rack, par, temp,midge_tubes, do_flux)

# pelagic
pelagic_clean = pelagic_full %>%
  mutate(duration = as.numeric(time_end - time_start),
         # calculate change in do
         # multiply by 1000 to convert to mg m^-3
         # multiply by column depth in m's to convert to mg m^-2
         # divde by 1000 to convert to g m^-2
         # divide by duration to convert to g m^-2 h-1
         do_flux = 1000*(column_depth/100)*(do_end - do_start)/1000/duration,
         # mean temperature
         temp = (temp_start + temp_end)/2) %>%
  select(id, sampledate, site, rack, par, temp, do_flux)





#==========
#========== Export
#==========

# write_csv(bethic_clean, "data/clean_data/benthic_clean.csv")
# write_csv(pelagic_clean, "data/clean_data/pelagic_clean.csv")





