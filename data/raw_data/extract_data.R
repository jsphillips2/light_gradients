#==========
#========== Preliminaries
#==========

# load packages
library(readxl)
library(readr)

# examine sheets
sheets = excel_sheets("data/raw_data/field_incubations_11aug18.xlsx")





#==========
#========== Import
#==========

# read sheets
benthic_grad = read_excel("data/raw_data/field_incubations_11aug18.xlsx",
                          sheet = "BenthicGradient")
pelagic_grad = read_excel("data/raw_data/field_incubations_11aug18.xlsx",
                          sheet = "Pelagic")
light_profile = read_excel("data/raw_data/field_incubations_11aug18.xlsx",
                          sheet = "LightProfiles")
shading = read_excel("data/raw_data/field_incubations_11aug18.xlsx",
                          sheet = "Shading")





#==========
#========== Export
#==========

# export sheets
write_csv(benthic_grad, "data/raw_data/benthic_grad.csv")
write_csv(pelagic_grad, "data/raw_data/pelagic_grad.csv")
write_csv(light_profile, "data/raw_data/light_profile.csv")
write_csv(shading, "data/raw_data/shading.csv")





