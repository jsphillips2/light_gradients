#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(lubridate)

# import data
benthic <- read_csv("data/raw_data/extracted/benthic_grad.csv") 
pelagic <- read_csv("data/raw_data/extracted/pelagic_grad.csv") 
profiles <- read_csv("data/raw_data/extracted/light_profile.csv") 
shading <- read_csv("data/raw_data/extracted/shading.csv") 





#==========
#========== Process data: Light (from Li-Cor meter)
#==========

# calculate effect of shading
shading_trt <- shading %>%
  mutate(light_frac = ifelse(par_in > par_out, 1, par_in/par_out)) %>%
  group_by(light_trt) %>%
  summarize(light_frac = mean(light_frac))

# benthic light
benthic_light <- benthic %>%
  group_by(sampledate, site) %>%
  summarize(time_start = min(time_start),
            time_end = max(time_end),
            depth = unique(inc_depth)) %>%
  left_join(profiles %>%
              filter(sampledate != "2018-06-27", site != "grim") %>%
              mutate(site = ifelse(site=="btl", "reyk", site),
                     site = ifelse(site == "kal", "st33", site))) %>%
  mutate(within = time > time_start - 60*30 & time < time_end + 60*30) %>%
  filter(within == T) %>%
  group_by(site, sampledate) %>%
  summarize(par = mean(par))

# pelagic light
pelagic_light <- pelagic %>%
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
  filter(within == T) %>%
  group_by(site, sampledate) %>%
  summarize(par = mean(par))

# add light data to benthic
benthic_full <- benthic %>%
  left_join(benthic_light) %>%
  left_join(shading_trt) %>%
  mutate(par = light_frac*par)

# add light data to benthic
pelagic_full <- pelagic %>%
  filter(site != "grim") %>%
  left_join(pelagic_light) %>%
  left_join(shading_trt) %>%
  mutate(par = light_frac*par) 





#==========
#========== Process data: Light (HOBO)
#==========

# import all files
hobo_names <- list.files("data/raw_data/hobo")
hobo <- hobo_names %>%
  lapply(function(x){read_csv(paste0("data/raw_data/hobo/", x), skip = 1) %>%
      mutate(site = str_split(x,"_")[[1]][1],
             loc = str_split(x,"_")[[1]][2])})
hobo[[1]]

# rename and select relevant columns
hobo_slim <- hobo %>%
  lapply(function(x){
    x %>%
      rename(date_time = grep("Date",names(hobo[[1]])),
             lux = grep("Intensity",names(hobo[[1]]))) %>%
      select(site, loc, date_time, lux) %>%
      return()
  }) %>%
  bind_rows()

# convert lux to par using standard correction for sunlight (divide by 54)
# approximate external par as par at 0.5cm (multiply by 0.44; based on analysis in sonde_oxygen)
hobo_clean = hobo_slim %>%
  mutate(par = lux/54,
         par = ifelse(loc == "out", 0.44*par, par))

# e5 benthic
e5_benthic <- benthic_full %>%
  filter(site == "e5") %>%
  group_by(sampledate, site) %>%
  summarize(time_start = min(time_start),
            time_end = max(time_end),
            depth = unique(inc_depth)) %>%
  left_join(hobo_clean %>%
              filter(loc == "out")) %>%
  mutate(within = date_time > time_start - 60*30 & date_time < time_end + 60*30) %>%
  filter(within == T) %>%
  rename(par0 = par) %>%
  select(-time_start, -time_end, -lux, -within, -loc) %>%
  left_join(hobo_clean %>%
              filter(loc == "bottom") %>%
              rename(parX = par) %>%
              select(-lux, - loc)) %>%
  mutate(atten = -log(parX/par0)/2.5,
         par = par0*exp(-mean(atten, na.rm=T)*0.5)) %>%
  group_by(sampledate, site) %>%
  summarize(par = mean(par))

# e5 pelagic
e5_pelagic <- pelagic_full %>%
  filter(site == "e5") %>%
  group_by(sampledate, site) %>%
  summarize(time_start = min(time_start),
            time_end = max(time_end),
            depth = unique(inc_depth)) %>%
  left_join(hobo_clean %>%
              filter(loc == "out")) %>%
  mutate(within = date_time > time_start - 60*30 & date_time < time_end + 60*30) %>%
  filter(within == T) %>%
  rename(par0 = par) %>%
  select(-time_start, -time_end, -lux, -within, -loc) %>%
  left_join(hobo_clean %>%
              filter(loc == "bottom") %>%
              rename(parX = par) %>%
              select(-lux, - loc)) %>%
  mutate(atten = -log(parX/par0)/2.5,
         par = par0*exp(-mean(atten, na.rm=T)*0.5)) %>%
  group_by(sampledate, site) %>%
  summarize(par = mean(par))

# add to full data
benthic_full <- benthic_full %>%
  left_join(e5_benthic %>% rename(par2 = par)) %>%
  left_join(shading_trt) %>%
  mutate(par2 = light_frac*par2,
         par = ifelse(sampledate=="2019-07-22", par2, par)) %>%
  select(-par2)

pelagic_full <- pelagic_full %>%
  left_join(e5_pelagic %>% rename(par2 = par)) %>%
  left_join(shading_trt) %>%
  mutate(par2 = light_frac*par2,
         par = ifelse(sampledate=="2019-07-22", par2, par)) %>%
  select(-par2)



#==========
#========== Process data: DO Flux
#==========

# benthic
benthic_clean <- benthic_full %>%
  mutate(duration = as.numeric(time_end - time_start),
         # calculate change in do
         # multiply by 1000 to convert to mg m^-3
         # multiply by column depth in m's to convert to mg m^-2
         # divde by 1000 to convert to g m^-2
         # divide by duration to convert to g m^-2 h-1
         do_flux = 1000*(column_depth/100)*(do_end - do_start)/1000/duration,
         # mean temperature
         temp = (temp_start + temp_end)/2) %>%
  select(id, sampledate, site, rack, par, temp, midge_tubes, do_flux)

# pelagic
pelagic_clean <- pelagic_full %>%
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

# write_csv(benthic_clean, "data/clean_data/benthic_clean.csv")
# write_csv(pelagic_clean, "data/clean_data/pelagic_clean.csv")





