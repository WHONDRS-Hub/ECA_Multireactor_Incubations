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
  mutate(Parent_ID = str_extract(Sample_ID, ".{5}(?=_)")) %>% 
  filter(!grepl("EC_21_MOI-3|EC_23_MOI-2|EC_24_MOI-3",Sample_ID))


#setting constants
#p = density of water (1000 kg/m^2)
p = 1000 
k = as.numeric(-6*(10^-20)) 

ssa_calc <- ssa %>% 
  mutate(gwc = (tray_soil_wt_g - tray_soil_od_wt_g)/(tray_soil_od_wt_g - tare_wt_g),
         Water_Potential_kPa = (Water_Potential_MPa*1000),
         S = (gwc/p) / ((k/(6 * pi* p * Water_Potential_kPa))^(1/3)),
         ssa_m2_g = S/1000) %>% 
  filter(ssa_m2_g > 0) #eliminate any values below 0 

#calculating Cv 
ssa_m2_g_by_site =
  ssa_calc %>% 
  group_by(Parent_ID) %>% 
  summarise(mean = mean(ssa_m2_g),
            sd = round(sd(ssa_m2_g),3),
            cv = round(sd/mean,3))

ssa_hist <- ggplot(data = ssa_m2_g_by_site, aes(x = mean, y = cv, color = Parent_ID)) +
  geom_jitter()

#histogram to show cv distribution
ssa_cv <- ggplot(data = ssa_m2_g_by_site, aes(x = cv)) +
  geom_histogram(binwidth = 0.05, aes(fill = factor(Parent_ID)), color = "black", alpha = 0.7)

#importing grain size
grain_size <- read.csv("data/v2_CM_SSS_Sediment_Grain_Size.csv", skip = 2, na=c("-9999", "N/A")) %>% 
  select(Sample_Name, Percent_Fine_Sand, Percent_Med_Sand, Percent_Coarse_Sand, Percent_Tot_Sand, Percent_Clay, Percent_Silt) %>% 
  mutate(Percent_Fine_Sand=as.numeric(Percent_Fine_Sand),
         Percent_Med_Sand=as.numeric(Percent_Med_Sand),
         Percent_Coarse_Sand=as.numeric(Percent_Coarse_Sand),
         Percent_Tot_Sand=as.numeric(Percent_Tot_Sand),
         Percent_Clay=as.numeric(Percent_Clay),
         Percent_Silt=as.numeric(Percent_Silt),
         Sample_Name=str_remove(Sample_Name, "_GRN"),
         Parent_ID = str_extract(Sample_Name, "\\d{2}$")) %>% 
  filter(!is.na(Sample_Name),
         !grepl("SSS", Sample_Name))
  
grain_size$Parent_ID <- paste0("EC_", grain_size$Parent_ID)

#joining grain size with ssa
gs_ssa <- grain_size %>% 
  left_join(ssa_m2_g_by_site, by = c('Parent_ID')) %>% 
  filter(!is.na(cv))

  
silt_cv_plot <- ggplot(data = gs_ssa, aes(x = cv, y = Percent_Silt, color = Parent_ID)) +
  geom_jitter()

clay_cv_plot <- ggplot(data = gs_ssa, aes(x = cv, y = Percent_Clay, color = Parent_ID)) +
  geom_jitter()
  
sand_cv_plot <- ggplot(data = gs_ssa, aes(x = cv, y = Percent_Tot_Sand, color = Parent_ID)) +
  geom_jitter() 
  
  
  
  
















