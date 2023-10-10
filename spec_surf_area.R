#### Code for calculating specific surface area from WP4C

# S = specific surface area (m^2/kg)
# p = density of water (1000 kg/m^2)
# k = hamaker constant (-6 * 10^-20 J)
# y = water potential (J/kg) 
  # note - water potential measured as MPa on WP4C
  # 1 J/kg = 1 kPa at p = 1000 kg/m^2
# = water content (g/g)

library(tidyverse)
library(readxl)
library(dplyr)

pnnl.user = 'guil098'

input.path = paste0("C:/Users/",pnnl.user,"/OneDrive - PNNL/Data Generation and Files/ECA/SSA/01_RawData/2023_Data_Raw_SSA_ECA_EC.xlsx")

ssa <- read_excel(input.path) %>% 
  mutate(Parent_ID = str_extract(Sample_ID, ".{5}(?=_)"))


#setting constants
#p = density of water (1000 kg/m^2)
p = 1000 
k = as.numeric(-6*(10^-20)) 

ssa_calc <- ssa %>% 
  mutate(gwc = (tray_soil_wt_g - tray_soil_od_wt_g)/(tray_soil_od_wt_g - tare_wt_g),
         Water_Potential_kPa = (Water_Potential_MPa*1000),
         S = (gwc/p) / ((k/(6 * pi* p * Water_Potential_kPa))^(1/3)),
         ssa_m2_g = S/1000) 


#calculating Cv 
pF_by_site = #water potential log
  ssa_calc %>% 
  group_by(Parent_ID) %>% 
  summarise(mean = mean(Water_Potential_log),
            sd = round(sd(Water_Potential_log),3),
            cv = round(sd/mean,3))

MPa_by_site =
  ssa_calc %>% 
  group_by(Parent_ID) %>% 
  summarise(mean = mean(Water_Potential_MPa),
            sd = round(sd(Water_Potential_MPa),3),
            cv = round(sd/mean,3))

S_by_site =
  ssa_calc %>% 
  group_by(Parent_ID) %>% 
  summarise(mean = mean(S),
            sd = round(sd(S),3),
            cv = round(sd/mean,3))

ssa_m2_g_by_site =
  ssa_calc %>% 
  group_by(Parent_ID) %>% 
  summarise(mean = mean(ssa_m2_g),
            sd = round(sd(ssa_m2_g),3),
            cv = round(sd/mean,3))



