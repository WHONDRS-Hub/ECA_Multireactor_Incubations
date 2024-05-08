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

respiration.date = '2024-04-05'
#respiration.summary = '2024-03-05'
grav.date = '2024-04-26'
grav.summary = '2024-04-26'
fe.date = '04-12-2024'
#fe.summary = '03-05-2024'
atp.date = '01-26-2024'
#atp.summary = '03-05-2024'

#Read in all data
setwd(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/"))

#All Respiration Rates
# Remove NEON samples - 52, 53, 57
all_respiration <- read.csv(paste0("INC/03_ProcessedData/ECA_Sediment_Incubations_Respiration_Rates_ReadyForBoye_",respiration.date,".csv")) %>% 
  select(c(Sample_Name, SpC, pH, Temp, Respiration_Rate_mg_DO_per_L_per_H, Respiration_Rate_mg_DO_per_kg_per_H, Methods_Deviation)) 

#mean_respiration <- read.csv(paste0("INC/03_ProcessedData/ECA_Sediment_Incubations_Respiration_Rates_Summary_ReadyForBoye_",respiration.summary,".csv"))

#ECA Iron
all_iron <- read_csv(paste0("Fe/03_ProcessedData/EC_ReadyForBoye_",fe.date,".csv")) %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "SFE", "INC")) %>% 
  select(-c(Methods_Deviation))
  

#mean_iron <- read_csv(paste0("Fe/03_ProcessedData/EC_SFE_Summary_ReadyForBoye_",fe.summary,".csv")) %>% 
  #select(-c(Material)) %>% 
  #mutate(Sample_Name = str_replace(Sample_Name, "SFE", "INC"))

#ICON Grain Size
grain <- read.csv("C:/Github/ECA_Multireactor_Incubations/Data/v3_CM_SSS_Sediment_Grain_Size.csv", skip = 2, header = TRUE)

grain_all <- grain %>% 
  filter(!row_number() %in% c(1:11)) %>% 
  dplyr::select(-c(Field_Name, Material)) %>% 
  filter(!grepl("SSS", Sample_Name)) %>% 
  filter(Sample_Name != "") %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "CM", "EC")) %>% 
  mutate(Sample_ID = str_remove(Sample_Name, "_GRN")) %>% 
  mutate_at(c("Percent_Fine_Sand", "Percent_Med_Sand", "Percent_Coarse_Sand", "Percent_Tot_Sand", "Percent_Silt", "Percent_Clay"), as.numeric) %>% 
  select(-c(IGSN, Methods_Deviation, Sample_Name))

ssa <- read.csv("C:/GitHub/ECA_Multireactor_Incubations/Data/CM_SSS_Sediment_Specific_Surface_Area.csv", skip = 2, header = TRUE)

ssa_clean = ssa %>% 
  filter(!row_number() %in% c(1:11)) %>% 
  dplyr::select(-c(Field_Name, Material, IGSN)) %>% 
  filter(!grepl("SSS", Sample_Name)) %>% 
  filter(Sample_Name != "") %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "CM", "EC")) %>% 
  separate(Sample_Name, c("Sample_ID", "Rep"), sep = -6)

mean_ssa <- ssa_clean %>% 
  filter(Specific_Surface_Area_m2_per_g != -9999) %>% 
  filter(!grepl("Negative", Specific_Surface_Area_m2_per_g)) %>%
  mutate(Specific_Surface_Area_m2_per_g = as.numeric(Specific_Surface_Area_m2_per_g)) %>%
  group_by(Sample_ID) %>%
  summarise(mean_ssa = mean(Specific_Surface_Area_m2_per_g, na.rm = TRUE))

#Gravimetric Moisture

grav_inc <- read.csv(paste0("INC/03_ProcessedData/EC_Drying_Masses_Summary_ReadyForBoye_on_",grav.summary,".csv")) %>% 
  select(-c(Methods_Deviation)) 

## ECA ATP ####
atp_all = read.csv(paste0("ATP/03_ProcessedData/EC_ATP_ReadyForBoye_",atp.date,".csv")) %>% 
  select(-c(Material, Methods_Deviation)) %>% 
  mutate(Sample_Name  = str_replace(Sample_Name, "ATP", "INC"))

#atp_summary = read.csv(paste0("ATP/03_ProcessedData/EC_ATP_Summary_ReadyForBoye_",atp.summary,".csv"))


##Start Merging Individual data

all_data <- left_join(all_respiration, all_iron, by = "Sample_Name") %>%
  left_join(grav_inc, by = "Sample_Name") %>% 
  left_join(atp_all, by = "Sample_Name") %>% 
  separate(Sample_Name, c("EC", "kit", "INC"), sep = "_", remove = FALSE) %>%
  unite(Sample_ID, c("EC", "kit")) %>% 
  left_join(grain_all, by = "Sample_ID") %>% 
  left_join(mean_ssa, by = "Sample_ID") %>% 
  mutate(Lost_Gravimetric_Water = Initial_Gravimetric_Moisture - Final_Gravimetric_Moisture)
 
write.csv(all_data,"C:/GitHub/ECA_Multireactor_Incubations/Data/Cleaned Data/All_ECA_Data_05-08-2024.csv", row.names = FALSE)  

# summary_data <- all_data %>% 
#   separate(Sample_Name, c("Sample_Name", "Replicate"), sep = "-") %>% 
#   separate(Replicate, c("Treat", "Replicate"), sep = -1) %>% 
#   unite(Sample_Name, c("Sample_Name", "Treat"), sep = "-") %>% 
#   dplyr::select(-c(Replicate)) %>% 
#   drop_na() %>% 
#   filter(Respiration_Rate_mg_DO_per_L_per_H > -9999) %>% 
#   group_by(Sample_Name) %>% 
#   summarise_if(is.numeric, mean)

medians = all_data %>% 
  filter(!grepl("INC_005|INC_Method_002|INC_008|INC_QA_004|INC_Method_001", Methods_Deviation)) %>% 
  filter(ATP_nanomol_per_L != -9999) %>% 
  filter(Respiration_Rate_mg_DO_per_kg_per_H != -9999) %>% 
  separate(Sample_Name, c("Sample_ID", "Rep"), sep = "-") %>% 
  mutate(Rep = if_else(grepl("D", Rep), "Dry", "Wet")) %>% 
  mutate(Fe_mg_per_L = if_else(grepl("SFE_Below", Fe_mg_per_L), "0.001", Fe_mg_per_L)) %>% 
  mutate(Fe_mg_per_kg = if_else(grepl("SFE_Below", Fe_mg_per_kg), "0.003", Fe_mg_per_kg)) %>%
  mutate(Fe_mg_per_L = if_else(grepl("SFE_Above", Fe_mg_per_L), str_extract(Fe_mg_per_L, "(?<=\\|[^|]{1,100}\\|)\\d+\\.\\d+"), Fe_mg_per_L)) %>% 
  mutate(Fe_mg_per_kg = if_else(grepl("SFE_Above", Fe_mg_per_kg), str_extract(Fe_mg_per_kg, "(?<=\\|[^|]{1,100}\\|)\\d+\\.\\d+"), Fe_mg_per_kg)) %>% 
  mutate(Fe_mg_per_L = as.numeric(Fe_mg_per_L)) %>% 
  mutate(Fe_mg_per_kg = as.numeric(Fe_mg_per_kg)) %>% 
  group_by(Sample_ID, Rep) %>%
  summarise(across(where(is.numeric), median)) %>% 
  rename_with(.cols = c(SpC:Lost_Gravimetric_Water), .fn = ~ paste0("median_", .x))
  
write.csv(medians,"C:/GitHub/ECA_Multireactor_Incubations/Data/Cleaned Data/Medians_ECA_Data.csv") 
  

effect_data <- medians %>% 
  #separate(Sample_Name, c("Sample_Name", "Treat"), sep = "-") %>% 
  filter(Sample_ID != "EC_011_INC") %>% 
  filter(Sample_ID != "EC_012_INC") %>% 
  relocate(median_Lost_Gravimetric_Water, .after = median_Final_Gravimetric_Moisture) %>% 
  #filter(Sample_Name != "EC_021_INC") %>% 
  group_by(Sample_ID) %>% 
  mutate(across(c(median_SpC:median_Lost_Gravimetric_Water), ~. [Rep == "Wet"] - .[Rep  == "Dry"])) %>% 
  rename_with(.cols = c(median_SpC:median_Lost_Gravimetric_Water), .fn = ~ paste0("diff_", .x)) %>% 
  distinct(Sample_ID, .keep_all = TRUE) %>% 
  select(-c(median_Incubation_Water_Mass_g, diff_median_Dry_Sediment_Mass_g, diff_median_Final_Water_mass_g, diff_median_Initial_Water_mass_g))
  

write.csv(effect_data,"C:/GitHub/ECA_Multireactor_Incubations/Data/Cleaned Data/Effect_Median_ECA_Data.csv") 
