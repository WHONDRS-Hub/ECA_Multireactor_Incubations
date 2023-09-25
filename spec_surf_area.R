#### Code for calculating specific surface area from WP4C

# S = specific surface area (m^2/kg)
# p = density of water (1000 kg/m^2)
# k = hamaker constant (-6 * 10^-20 J)
# y = water potential (J/kg) 
  # note - water potential measured as MPa on WP4C
  # 1 J/kg = 1 kPa at p = 1000 kg/m^2
# = water content (g/g)

library(readxl)
library(dplyr)

pnnl.user = 'laan208'

input.path = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/SSA/01_RawData/2023_Data_Raw_SSA_ECA_EC.xlsx")

ssa <- read_excel(input.path)  



ssa_calc <- ssa %>% 
  mutate(gwc = (tray_soil_wt_g - tray_soil_od_wt_g)/(tray_soil_od_wt_g - tare_wt_g)) %>% 
  mutate(Water_Potential_kPa = `Water_Potential_MPa `*1000) %>% 
  mutate(k = as.numeric(-6*(10^-20))) %>% 
  mutate(p = 1000) %>% 
  mutate(S = (gwc/((k/(6 * pi* (p * Water_Potential_kPa)))^(1/3))*p))


