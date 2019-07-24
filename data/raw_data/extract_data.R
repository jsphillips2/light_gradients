#==========
#========== Preliminaries
#==========

# load packages
library(readxl)
library(readr)

# set path
path <- "data/raw_data/field_incubations_24jul19.xlsx"

# examine sheets
sheets <- excel_sheets(path)





#==========
#========== Import
#==========

# read sheets
benthic_grad <- read_excel(path, sheet = "BenthicGradient", na=c("","NA"))
pelagic_grad <- read_excel(path, sheet = "Pelagic", na=c("","NA"))
light_profile <- read_excel(path, sheet = "LightProfiles", na=c("","NA"))
shading <- read_excel(path, sheet = "Shading", na=c("","NA"))
hobo_log <- read_excel(path, sheet = "HOBOLog", na=c("","NA")) 

# correct light meter readings taken out in air (-0.05m) prior to 2019
# in-water calibration multiplier for UW-LTREB Li-Cor meter is -344.88
# in-air calibration multiplier for UW-LTREB Li-Cor meter is -261.28
light_profile <- light_profile %>%
  mutate(par = ifelse(year(sampledate) < 2019, (-261.28/-344.88)*par, par))





#==========
#========== Export
#==========

# export sheets
# write_csv(benthic_grad, "data/raw_data/extracted/benthic_grad.csv")
# write_csv(pelagic_grad, "data/raw_data/extracted/pelagic_grad.csv")
# write_csv(light_profile, "data/raw_data/extracted/light_profile.csv")
# write_csv(shading, "data/raw_data/extracted/shading.csv")
# write_csv(hobo_log, "data/raw_data/extracted/hobo_log.csv")





