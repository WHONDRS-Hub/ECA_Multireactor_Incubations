##ECA Data Analysis 

#### Read in Data ####

#read in libraries

library(lubridate);library(writexl);library(raster);library(tidyverse);library(devtools)
library(readxl)
library(corrplot)
library(corrr)
library(vegan)
library(FactoMineR)
library(factoextra)
library(ggpmisc)

rm(list=ls());graphics.off()

# Set working directory to data file
#Example:
pnnl.user = 'laan208'

# choose file dates to read in 

respiration.date = '2024-03-05'
respiration.summary = '2024-03-05'
grav.date = '2024-01-29'
grav.summary = '01-29-2024'
fe.date = '03-04-2024'
fe.summary = '03-05-2024'
atp.date = '01-26-2024'
atp.summary = '03-05-2024'

#Read in all data
setwd(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/"))

#All Respiration Rates
all_respiration <- read.csv(paste0("INC/03_ProcessedData/ECA_Sediment_Incubations_Respiration_Rates_ReadyForBoye_",respiration.date,".csv"))

mean_respiration <- read.csv(paste0("INC/03_ProcessedData/ECA_Sediment_Incubations_Respiration_Rates_Summary_ReadyForBoye_",respiration.summary,".csv"))

#ECA Iron
all_iron <- read_csv(paste0("Fe/03_ProcessedData/EC_ReadyForBoye_",fe.date,".csv"))

mean_iron <- read_csv(paste0("Fe/03_ProcessedData/EC_SFE_Summary_ReadyForBoye_",fe.summary,".csv"))

#ICON Grain Size
grain <- read.csv("C:/Github/ECA_Multireactor_Incubations/Data/v3_CM_SSS_Sediment_Grain_Size.csv", skip = 2, header = TRUE)

grain_all <- grain %>% 
  filter(!row_number() %in% c(1:11)) %>% 
  dplyr::select(-c(Field_Name, Material)) %>% 
  filter(!grepl("SSS", Sample_Name)) %>% 
  filter(Sample_Name != "") %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "CM", "EC")) %>% 
  mutate(Sample_Name = str_remove(Sample_Name, "_GRN")) %>% 
  mutate_at(c("Percent_Fine_Sand", "Percent_Med_Sand", "Percent_Coarse_Sand", "Percent_Tot_Sand", "Percent_Silt", "Percent_Clay"), as.numeric)

ssa <- read.csv("C:/GitHub/ECA_Multireactor_Incubations/Data/CM_SSS_Sediment_Specific_Surface_Area.csv", skip = 2, header = TRUE)

ssa_clean = ssa %>% 
  filter(!row_number() %in% c(1:11)) %>% 
  dplyr::select(-c(Field_Name, Material, IGSN)) %>% 
  filter(!grepl("SSS", Sample_Name)) %>% 
  filter(Sample_Name != "") %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "CM", "EC")) %>% 
  separate(Sample_Name, c("Sample_Name", "Rep"), sep = -6)

mean_ssa <- ssa_clean %>% 
  filter(Specific_Surface_Area_m2_per_g != -9999) %>% 
  filter(!grepl("Negative", Specific_Surface_Area_m2_per_g)) %>%
  mutate(Specific_Surface_Area_m2_per_g = as.numeric(Specific_Surface_Area_m2_per_g)) %>%
  group_by(Sample_Name) %>%
  summarise(mean_ssa = mean(Specific_Surface_Area_m2_per_g, na.rm = TRUE))

#Gravimetric Moisture

grav_inc <- read.csv(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/INC/03_ProcessedData/ECA_Drying_Masses_Summary_merged_by_laan208_on_",grav.date,".csv"))

all_data <- left_join(all_respiration, iron_samples, by = "Sample_Name") %>% 
  separate(Sample_Name, c("EC", "kit", "INC"), sep = "_", remove = FALSE) %>%
  unite(Sample_ID, c("EC", "kit")) %>% 
  left_join(grain_all, by = "Sample_ID") %>% 
  left_join(mean_ssa, by = "Sample_ID") %>% 
  left_join(chem_all, by = "Sample_Name") %>% 
  left_join(grav_inc, by = "Sample_Name") %>% 
  mutate(Fe_mg_per_kg = Fe_mg_per_L * (mass_water/Dry_Sediment_Mass_g)) %>% 
  relocate(Fe_mg_per_kg, .after = Fe_mg_per_L) %>% 
  mutate(Initial_Gravimetric_Water = Initial_Water_mass_g/Dry_Sediment_Mass_g) %>% 
  mutate(Final_Gravimetric_Water = Final_Water_mass_g/Dry_Sediment_Mass_g) %>% 
  mutate(Lost_Gravimetric_Water = Initial_Gravimetric_Water - Final_Gravimetric_Water) %>% 
dplyr::select(-c(Methods_Deviation, ECA, kit, Analysis, INC, Treat, mass_water)) 
 
write.csv(all_data,"C:/GitHub/ECA_Multireactor_Incubations/Data/Cleaned Data/All_ECA_Data.csv")  

summary_data <- all_data %>% 
  separate(Sample_Name, c("Sample_Name", "Replicate"), sep = "-") %>% 
  separate(Replicate, c("Treat", "Replicate"), sep = -1) %>% 
  unite(Sample_Name, c("Sample_Name", "Treat"), sep = "-") %>% 
  dplyr::select(-c(Replicate)) %>% 
  drop_na() %>% 
  filter(Respiration_Rate_mg_DO_per_L_per_H > -9999) %>% 
  group_by(Sample_Name) %>% 
  summarise_if(is.numeric, mean)
  
write.csv(summary_data,"C:/GitHub/ECA_Multireactor_Incubations/Data/Cleaned Data/Summary_ECA_Data.csv") 
  

effect_data <- summary_data %>% 
  separate(Sample_Name, c("Sample_Name", "Treat"), sep = "-") %>% 
  filter(Sample_Name != "EC_011_INC") %>% 
  filter(Sample_Name != "EC_012_INC") %>% 
  filter(Sample_Name != "EC_021_INC") %>% 
  group_by(Sample_Name) %>% 
  mutate(Effect_Size_mg_per_L = (abs(Respiration_Rate_mg_DO_per_L_per_H[Treat == "W"]) - abs(Respiration_Rate_mg_DO_per_L_per_H[Treat == "D"]))) %>% 
  mutate(Effect_Size_mg_per_kg = (abs(Respiration_Rate_mg_DO_per_kg_per_H[Treat == "W"]) - abs(Respiration_Rate_mg_DO_per_kg_per_H[Treat == "D"]))) %>% 
  mutate(Fe_Difference_mg_per_L = Fe_mg_per_L[Treat == "W"] - Fe_mg_per_L[Treat == "D"]) %>% 
  mutate(Fe_Difference_mg_per_kg = Fe_mg_per_kg[Treat == "W"] - Fe_mg_per_kg[Treat == "D"]) %>% 
  mutate(SpC_Difference = SpC[Treat == "W"] - SpC[Treat == "D"]) %>% 
  mutate(pH_Difference = pH[Treat == "W"] - pH[Treat == "D"]) %>% 
  mutate(Temp_Difference = Temp[Treat == "W"] - Temp[Treat == "D"]) %>% 
  mutate(Initial_Grav_Water_Difference = Initial_Gravimetric_Water[Treat == "W"] - Initial_Gravimetric_Water[Treat == "D"]) %>% 
  mutate(Final_Grav_Water_Difference = Final_Gravimetric_Water[Treat == "W"] - Final_Gravimetric_Water[Treat == "D"]) %>% 
  mutate(Final_Grav_Water_Difference = Final_Gravimetric_Water[Treat == "W"] - Final_Gravimetric_Water[Treat == "D"]) %>% 
  distinct(Sample_Name, .keep_all = TRUE) %>% 
  dplyr::select(c(Sample_Name, Effect_Size_mg_per_L, Effect_Size_mg_per_kg, Fe_Difference_mg_per_L, Fe_Difference_mg_per_kg, SpC_Difference, pH_Difference, Temp_Difference, Initial_Grav_Water_Difference, Final_Grav_Water_Difference))

write.csv(effect_data,"C:/GitHub/ECA_Multireactor_Incubations/Data/Cleaned Data/Effect_ECA_Data.csv") 
