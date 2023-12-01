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

effect.date = '2023-11-08'
respiration.date = '2023-11-08'
removed.respiration.date = '2023-11-13'
mg.kg.respiration.date = '2023-12-01'
grav.date = '2023-12-01'

#Read in all data
setwd(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/"))

#effect size - change date to most recent
effect_size <- read_csv(paste0("Optode multi reactor/Optode_multi_reactor_incubation/rates/ReadyForBoye/ECA_Effect_Size_ReadyForBoye_",effect.date,".csv"))


#Respiration rates with removals from dist matrix to calculate effect size 

#respiration <- read.csv(paste0("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/rates/ReadyForBoye/ECA_Sediment_Incubations_Removed_Respiration_Rates_",removed.respiration.date,".csv"))

#all_respiration <- read.csv(paste0("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/rates/ReadyForBoye/ECA_Sediment_Incubations_Respiration_Rates_ReadyForBoye_",respiration.date,".csv"))

all_respiration <- read.csv(paste0("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/rates/ReadyForBoye/ECA_Sediment_Incubations_mg_kg_rates_laan208_on_",mg.kg.respiration.date,".csv"))

all_respiration <- all_respiration %>% 
  dplyr::select(c(Sample_Name, Respiration_Rate_mg_DO_per_L_per_H, Respiration_Rate_mg_DO_per_kg_per_H, mass_water))

#ECA Iron
iron <- read_csv(paste0("Fe/03_ProcessedData/EC_SFE_ReadyForBoye_06-29-2023.csv"))

iron <- iron %>% 
  dplyr::select(c(sample_label, Fe_mg_per_L))

iron_ex <- read_csv(paste0("Fe/03_ProcessedData/20230711_Data_Processed_SFE_ECA_EC_1-153/20230711_Data_Processed_SFE_ECA_EC_1-153.csv"))

iron_ex <- iron_ex %>% 
  rename(Fe_mg_per_L = mg_Fe_per_L) %>%
  filter(is.na(flag)) %>% 
  dplyr::select(c(sample_label, Fe_mg_per_L))

iron <- rbind(iron, iron_ex)

iron <- iron %>% 
separate(sample_label, into = c("EC", "Site", "INC"), sep = "_", remove = FALSE)

#add "0" to start of sample kit names that don't have it

for (i in 1:nrow(iron)){
  
  if (str_count(iron$Site[i], "[0-9]") <= 2){
    
    iron$Site[i] = paste0("0", iron$Site[i])
    
  }
  
  else {
    
    iron$Site[i] = iron$Site[i]
  }
  
}

iron_samples <- iron %>% 
  drop_na(sample_label) %>% 
  separate(INC, c("Replicate", "Analytical"), sep = -1) %>% 
  group_by(Site, Replicate) %>% 
  mutate(Fe_mg_per_L = mean(Fe_mg_per_L)) %>% 
  unite(Sample_Name, c(EC:Replicate), sep = "_") %>% 
  distinct(Sample_Name, .keep_all = TRUE) %>% 
  dplyr::select(-c(Analytical, sample_label)) %>% 
  filter(Fe_mg_per_L > 0.03) %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "SFE", "INC"))
  

#ICON Grain Size
grain <- read.csv("C:/Github/ECA_Multireactor_Incubations/Data/v2_CM_SSS_Sediment_Grain_Size.csv", skip = 2, header = TRUE)

grain <- grain %>% 
  filter(!row_number() %in% c(1:11)) %>% 
  dplyr::select(-c(Field_Name, Material)) 

grain_new <- read.csv(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ICON_ModEx_SSS/09_Grain_Size/03_ProcessedData/20230721_Grain_Size_SBR_RC4_CM_R21/20230721_Data_Processed_Grain_Size_SBR_RC4_CM_R21.csv"))

grain <- grain %>% 
  mutate(Sample_ID = Sample_Name) %>% 
  dplyr::select(-c(Sample_Name))

grain_all <- rbind(grain, grain_new)

grain_all$Sample_ID <- str_replace(grain_all$Sample_ID, "CM", "EC") 

grain_all <- grain_all %>% 
  filter(!grepl("SSS", Sample_ID)) %>% 
  separate(Sample_ID, into = c("EC", "Kit", "GRN"), sep = "_") %>% 
  unite("Sample_ID", EC:Kit, sep = "_") %>% 
  dplyr::select(-c(GRN)) %>% 
  filter(!grepl("NA", Sample_ID))


ssa <- read_csv(paste0("C:/GitHub/ECA_Multireactor_Incubations/Data/eca_ssa_predatapackage.csv"))

ssa <- ssa %>% 
  separate(Parent_ID, c("EC", "Site"), sep = "_")

for (i in 1:nrow(ssa)){
  
  if (str_count(ssa$Site[i], "[0-9]") <= 2){
    
    ssa$Site[i] = paste0("0", ssa$Site[i])
    
  }
  
  else {
    
    ssa$Site[i] = ssa$Site[i]
  }
  
}

ssa <- ssa %>% 
  unite(Sample_ID, c("EC", "Site"), sep = "_")

mean_ssa <- ssa %>% 
  group_by(Sample_ID) %>% 
  mutate(average_ssa = mean(ssa_m2_g)) %>% 
  distinct(Sample_ID, .keep_all = TRUE) %>% 
  dplyr::select(Sample_ID, average_ssa)



#All incubation pH, SpC, temp
chemistry <- read_csv("INC/03_ProcessedData/SpC_pH_Temp.csv")

map_corr = chemistry %>% 
  dplyr::select(c(Sample_Name, SpC, Temp, pH)) %>% 
  filter(!grepl("EV", Sample_Name))

chem_all = map_corr %>% 
  separate(Sample_Name, c("ECA"
                          , "kit", "Analysis"), sep = "_", remove = FALSE) %>% 
  mutate(Treat = case_when(grepl("W",Analysis)~"Wet",
                           grepl("D", Analysis) ~"Dry"))


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
select(-c(Methods_Deviation, ECA, kit, Analysis, INC, Treat, mass_water)) 
 
write.csv(all_data,"C:/GitHub/ECA_Multireactor_Incubations/Data/Cleaned Data/All_ECA_Data.csv")  
