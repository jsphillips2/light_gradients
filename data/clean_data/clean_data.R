#=========================================================================================
#========== Preliminaries
#=========================================================================================

# load packages
library(tidyverse)
library(lubridate)

# import data
benthic <- read_csv("data/raw_data/extracted/benthic_grad.csv") 
pelagic <- read_csv("data/raw_data/extracted/pelagic_grad.csv") 
profiles <- read_csv("data/raw_data/extracted/light_profile.csv") 
shading <- read_csv("data/raw_data/extracted/shading.csv") 

# correct light meter readings taken out in air (-0.05m) prior to 2019
# in-water calibration multiplier for UW-LTREB Li-Cor meter is -344.88
# in-air calibration multiplier for UW-LTREB Li-Cor meter is -261.28
profiles <- profiles %>%
  mutate(par = ifelse(year(sampledate) < 2019, (-261.28/-344.88)*par, par))

#=========================================================================================





#=========================================================================================
#========== Light attenuation & light at depth (from Li-Cor)
#=========================================================================================

# attenuation
light_atten <- profiles %>%
  mutate(par = ifelse(par == 0, 0.5, par)) %>%
  group_by(time, site) %>%
  filter(depth >= 0, length(na.omit(par) > 1)) %>%
  ungroup() %>%
  split(.$time) %>%
  lapply(function(x){split(x, x$site)}) %>%
  lapply(function(x){xx = x[[1]]
  m = lm(log(par) ~ depth, data = xx)
  xxx = xx %>%
    group_by(sampledate, time, site) %>%
    summarize(par0 = exp(coef(m)[1]),
              atten = coef(m)[2])
  return(xxx)
  }) %>%
  bind_rows() 

# expand by depths
light_depth <-light_atten %>%
  expand(sampledate, time, site, par0, atten, depth = profiles$depth) %>%
  mutate(par = par0*exp(atten*depth))

#=========================================================================================





#=========================================================================================
#========== Add light data to metabolism data
#=========================================================================================

# calculate effect of shading
shading_trt <- shading %>%
  mutate(light_frac = ifelse(par_in > par_out, 1, par_in/par_out)) %>%
  group_by(light_trt) %>%
  summarize(light_frac = mean(light_frac))

# benthic light
benthic_light <- benthic %>%
  group_by(sampledate, site) %>%
  summarize(time_start = min(time_start,na.rm=T),
            time_end = max(time_end,na.rm=T),
            depth = unique(inc_depth)) %>%
  left_join(light_depth %>%
              filter(sampledate != "2018-06-27", site != "grim") %>%
              mutate(site = ifelse(site=="btl", "reyk", site),
                     site = ifelse(site == "kal", "st33", site))) %>%
  mutate(within = time > time_start - 60*30 & time < time_end + 60*30) %>%
  filter(within == T) %>%
  group_by(site, sampledate) %>%
  summarize(par = mean(par,na.rm=T))

# pelagic light
pelagic_light <- pelagic %>%
  filter(is.na(time_start)==F, is.na(time_end)==F, site != "grim") %>%
  group_by(sampledate, site) %>%
  summarize(time_start = min(time_start,na.rm=T),
            time_end = max(time_end,na.rm=T),
            depth = unique(inc_depth)) %>%
  left_join(light_depth %>%
              filter(sampledate != "2018-06-27", site != "grim") %>%
              mutate(site = ifelse(site=="btl", "reyk", site),
                     site = ifelse(site == "kal", "st33", site))) %>%
  mutate(within = time > time_start - 60*30 & time < time_end + 60*30) %>%
  filter(within == T) %>%
  group_by(site, sampledate) %>%
  summarize(par = mean(par,na.rm=T))

# add light data to benthic
benthic_full <- benthic %>%
  left_join(benthic_light)

# add light data to benthic
pelagic_full <- pelagic %>%
  filter(site != "grim") %>%
  left_join(pelagic_light)

#=========================================================================================




#=========================================================================================
#========== HOBO Light (for select dates)
#=========================================================================================

# estimate decline in light immediately beneath surface
air_surf <- profiles %>%
  mutate(depth = ifelse(depth == 0.05, 0, depth)) %>%
  filter(depth %in% c(-0.05, 0)) %>%
  spread(depth, par) %>%
  rename(air = `-0.05`,
         surf = `0`) %>%
  na.omit() %>%
  {lm(surf ~ air - 1, data = .)}
coef(air_surf)

# import all files
hobo_names <- list.files("data/raw_data/hobo")
hobo <- hobo_names %>%
  lapply(function(x){read_csv(paste0("data/raw_data/hobo/", x), skip = 1) %>%
      mutate(site = str_split(x,"_")[[1]][1],
             type = str_split(x,"_")[[1]][2])})
hobo[[1]]

# rename and select relevant columns
hobo_clean <- hobo %>%
  lapply(function(x){
    x %>%
      rename(date_time = grep("Date",names(hobo[[1]])),
             lux = grep("Intensity",names(hobo[[1]]))) %>%
      select(site, type, date_time, lux) %>%
      return()
  }) %>%
  bind_rows() %>%
  mutate(sampledate = as_date(date_time)) %>%
  # combine with log file, only keeping mathces
  inner_join(read_csv("data/raw_data/extracted/hobo_log.csv") %>%
               mutate(sampledate = as_date(sampledate))) %>%
  # convert lux to par using standard correction for sunlight (divide by 54)
  # approximate external par as par at water surface
  mutate(par = lux/54,
         par = ifelse(type == "out", coef(air_surf)*par, par))

#=========================================================================================





#=========================================================================================
#========== Add HOBO light (for select dates)
#=========================================================================================


# Reykjahlid 17 July 2018
# benthic
reyk_benthic18 <- benthic_full %>%
  filter(site == "reyk" & sampledate == "2018-07-17") %>%
  group_by(sampledate, site) %>%
  summarize(time_start = min(time_start),
            time_end = max(time_end),
            depth = unique(inc_depth)) %>%
  left_join(hobo_clean %>%
              select(-depth) %>%
              filter(type == "out")) %>%
  mutate(within = date_time > time_start - 60*30 & date_time < time_end + 60*30) %>%
  filter(within == T) %>%
  group_by(site, sampledate) %>%
  summarize(par0 = mean(par)) %>%
  left_join(light_atten %>% 
              filter(site == "reyk" & sampledate == "2018-07-17") %>%
              group_by(site, sampledate) %>%
              summarize(atten = mean(atten))) %>%
  mutate(par = par0*exp(atten*0.5)) %>%
  select(site, sampledate, par)

# pelagic
reyk_pelagic18 <- pelagic_full %>%
  filter(site == "reyk" & sampledate == "2018-07-17") %>%
  group_by(sampledate, site) %>%
  summarize(time_start = min(time_start),
            time_end = max(time_end),
            depth = unique(inc_depth)) %>%
  left_join(hobo_clean %>%
              select(-depth) %>%
              filter(type == "out")) %>%
  mutate(within = date_time > time_start - 60*30 & date_time < time_end + 60*30) %>%
  filter(within == T) %>%
  group_by(site, sampledate) %>%
  summarize(par0 = mean(par)) %>%
  left_join(light_atten %>% 
              filter(site == "reyk" & sampledate == "2018-07-17") %>%
              group_by(site, sampledate) %>%
              summarize(atten = mean(atten))) %>%
  mutate(par = par0*exp(atten*0.5)) %>%
  select(site, sampledate, par)



# e5 22-July-2019 
# benthic
e5_benthic19 <- benthic_full %>%
  filter(site == "e5" & sampledate == "2019-07-22") %>%
  group_by(sampledate, site) %>%
  summarize(time_start = min(time_start),
            time_end = max(time_end),
            depth = unique(inc_depth)) %>%
  left_join(hobo_clean %>%
              select(-depth) %>%
              filter(type == "out")) %>%
  mutate(within = date_time > time_start - 60*30 & date_time < time_end + 60*30) %>%
  filter(within == T) %>%
  rename(par0 = par) %>%
  select(-time_start, -time_end, -lux, -within, -type) %>%
  left_join(hobo_clean %>%
              select(-depth) %>%
              filter(type == "bottom") %>%
              rename(parX = par) %>%
              select(-lux, - type)) %>%
  mutate(atten = -log(parX/par0)/2.5,
         par = par0*exp(-mean(atten, na.rm=T)*0.5)) %>%
  group_by(sampledate, site) %>%
  summarize(par = mean(par)) %>%
  select(site, sampledate, par)

# pelagic
e5_pelagic19 <- pelagic_full %>%
  filter(site == "e5" & sampledate == "2019-07-22") %>%
  filter(site == "e5") %>%
  group_by(sampledate, site) %>%
  summarize(time_start = min(time_start),
            time_end = max(time_end),
            depth = unique(inc_depth)) %>%
  left_join(hobo_clean %>%
              select(-depth) %>%
              filter(type == "out")) %>%
  mutate(within = date_time > time_start - 60*30 & date_time < time_end + 60*30) %>%
  filter(within == T) %>%
  rename(par0 = par) %>%
  select(-time_start, -time_end, -lux, -within, -type) %>%
  left_join(hobo_clean %>%
              select(-depth) %>%
              filter(type == "bottom") %>%
              rename(parX = par) %>%
              select(-lux, - type)) %>%
  mutate(atten = -log(parX/par0)/2.5,
         par = par0*exp(-mean(atten, na.rm=T)*0.5)) %>%
  group_by(sampledate, site) %>%
  summarize(par = mean(par)) %>%
  select(site, sampledate, par)



# add to full data
benthic_full <- benthic_full %>%
  left_join(e5_benthic19 %>% rename(par2 = par)) %>%
  left_join(reyk_benthic18 %>% rename(par3 = par)) %>%
  left_join(shading_trt) %>%
  mutate(par2 = light_frac*par2,
         par3 = light_frac*par3,
         par = ifelse(is.na(par2)==F, par2, ifelse(is.na(par3)==F, par3, par)),
         par = light_frac*par) %>%
  select(-par2, -par3)

pelagic_full <- pelagic_full %>%
  left_join(e5_pelagic19 %>% rename(par2 = par)) %>%
  left_join(reyk_pelagic18 %>% rename(par3 = par)) %>%
  left_join(shading_trt) %>%
  mutate(par2 = light_frac*par2,
         par3 = light_frac*par3,
         par = ifelse(is.na(par2)==F, par2, ifelse(is.na(par3)==F, par3, par)),
         par = light_frac*par) %>%
  select(-par2, -par3)

#=========================================================================================





#=========================================================================================
#========== Process data: DO Flux
#=========================================================================================

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

#=========================================================================================




#=========================================================================================
#========== Export
#=========================================================================================

# write_csv(benthic_clean, "data/clean_data/benthic_clean.csv")
# write_csv(pelagic_clean, "data/clean_data/pelagic_clean.csv")





