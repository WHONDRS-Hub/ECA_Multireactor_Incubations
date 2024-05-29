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
npoc.tn.date = '2024-03-01'
cn.date = '05-23-2024'

#Read in all data
setwd(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/"))

#All Respiration Rates
# Remove NEON samples - 52, 53, 57
all_respiration <- read.csv(paste0("INC/03_ProcessedData/ECA_Sediment_Incubations_Respiration_Rates_ReadyForBoye_",respiration.date,".csv")) %>% 
  dplyr::select(c(Sample_Name, SpC, pH, Temp, Respiration_Rate_mg_DO_per_L_per_H, Respiration_Rate_mg_DO_per_kg_per_H, Methods_Deviation)) 

#mean_respiration <- read.csv(paste0("INC/03_ProcessedData/ECA_Sediment_Incubations_Respiration_Rates_Summary_ReadyForBoye_",respiration.summary,".csv"))

#ECA Iron
all_iron <- read_csv(paste0("Fe/03_ProcessedData/EC_ReadyForBoye_",fe.date,".csv")) %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "SFE", "INC")) %>% 
  dplyr::select(-c(Methods_Deviation))
  

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
  dplyr::select(-c(IGSN, Methods_Deviation, Sample_Name))

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
  dplyr::select(-c(Methods_Deviation)) 

## ECA ATP ####
atp_all = read.csv(paste0("ATP/03_ProcessedData/EC_ATP_ReadyForBoye_",atp.date,".csv")) %>% 
  dplyr::select(-c(Material, Methods_Deviation)) %>% 
  mutate(Sample_Name  = str_replace(Sample_Name, "ATP", "INC"))

#atp_summary = read.csv(paste0("ATP/03_ProcessedData/EC_ATP_Summary_ReadyForBoye_",atp.summary,".csv"))

## ECA NPOC/TN ####

npoc_tn_all = read.csv(paste0("Boye_Files/EC/EC_NPOC_TN_Check_for_Duplicates_",npoc.tn.date,"_by_laan208.csv")) %>% 
  mutate(Sample_Name  = str_replace(Sample_ID, "SIR", "INC")) %>% 
  dplyr::select(-c(Date_of_Run, Methods_Deviation, Method_Notes, duplicate, Sample_ID)) %>% 
  relocate(Sample_Name, .before = NPOC_mg_C_per_L)

## ECA C/N

cn_all = read.csv(paste0("CN/02_FormattedData/ECA_CN_",cn.date,".csv")) %>% 
  mutate(Sample_Name  = str_replace(sample_id, "SCN", "INC")) %>% 
  separate(Sample_Name, c("Parent_ID", "Rep"), remove = TRUE, sep = "_INC") %>% 
  unite(Sample_Name, c("parent_id", "Rep"), sep = "_INC") %>% 
  dplyr::select(-c(sample_id, Parent_ID))
##Start Merging Individual data

all_data <- left_join(all_respiration, all_iron, by = "Sample_Name") %>%
  left_join(grav_inc, by = "Sample_Name") %>% 
  left_join(atp_all, by = "Sample_Name") %>% 
  left_join(npoc_tn_all, by = "Sample_Name") %>% 
  left_join(cn_all, by = "Sample_Name") %>% 
  separate(Sample_Name, c("EC", "kit", "INC"), sep = "_", remove = FALSE) %>%
  unite(Sample_ID, c("EC", "kit")) %>% 
  left_join(grain_all, by = "Sample_ID") %>% 
  left_join(mean_ssa, by = "Sample_ID") %>% 
  mutate(Lost_Gravimetric_Water = Initial_Gravimetric_Moisture - Final_Gravimetric_Moisture)%>% 
  mutate(Respiration_Rate_mg_DO_per_L_per_H = if_else(Respiration_Rate_mg_DO_per_L_per_H == -9999, -9999,abs(Respiration_Rate_mg_DO_per_L_per_H))) %>%
  mutate(Respiration_Rate_mg_DO_per_kg_per_H = if_else(Respiration_Rate_mg_DO_per_kg_per_H == -9999, -9999, abs(Respiration_Rate_mg_DO_per_kg_per_H)))%>% 
  filter(!grepl("EC_052|EC_053|EC_057|EC_023", Sample_Name))%>% #remove NEON sites after 1st incubation, no gravimetric moisture (EC_023)
  filter(!grepl("INC_005|INC_Method_002|INC_008|INC_QA_004|INC_Method_001", Methods_Deviation)) %>% #remove samples with too much water (EC_011/012-W), missing replicates (EC_072-W5/D5), spilled sample (EC_041-W4), overexposed samples (EC_027, EC_013, EC_014), less sediment in sample (EC_012-D5), no gravimetric water (EC_023)
  mutate(Fe_mg_per_L = if_else(grepl("SFE_Below", Fe_mg_per_L), "0.001", Fe_mg_per_L)) %>% 
  mutate(Fe_mg_per_L = if_else(grepl("SFE_Above", Fe_mg_per_L), str_extract(Fe_mg_per_L, "(?<=\\|[^|]{1,100}\\|)\\d+\\.\\d+"), Fe_mg_per_L)) %>% 
  mutate(Fe_mg_per_L = as.numeric(Fe_mg_per_L)) %>%
  mutate(Fe_mg_per_kg = if_else(grepl("SFE_Above", Fe_mg_per_kg), str_extract(Fe_mg_per_kg, "(?<=\\|[^|]{1,100}\\|)\\d+\\.\\d+"), Fe_mg_per_kg)) %>% 
  mutate(Fe_mg_per_kg = if_else(grepl("SFE_Below", Fe_mg_per_kg), as.numeric(Fe_mg_per_L * (Incubation_Water_Mass_g/Dry_Sediment_Mass_g)), as.numeric(Fe_mg_per_kg))) %>%
  mutate(Fe_mg_per_kg = as.numeric(Fe_mg_per_kg))

write.csv(all_data,"C:/GitHub/ECA_Multireactor_Incubations/Data/Cleaned Data/All_ECA_Data_05-29-2024.csv", row.names = FALSE)  

# summary_data <- all_data %>% 
#   separate(Sample_Name, c("Sample_Name", "Replicate"), sep = "-") %>% 
#   separate(Replicate, c("Treat", "Replicate"), sep = -1) %>% 
#   unite(Sample_Name, c("Sample_Name", "Treat"), sep = "-") %>% 
#   dplyr::select(-c(Replicate)) %>% 
#   drop_na() %>% 
#   filter(Respiration_Rate_mg_DO_per_L_per_H > -9999) %>% 
#   group_by(Sample_Name) %>% 
#   summarise_if(is.numeric, mean)

# cv = all_data %>% 
#   separate(Sample_Name, c("Sample_ID", "Rep"), sep = "-") %>% 
#   mutate(Rep = if_else(grepl("D", Rep), "Dry", "Wet")) %>%
#   mutate(NPOC_mg_C_per_L = as.numeric(NPOC_mg_C_per_L)) %>% 
#   mutate(TN_mg_N_per_L = as.numeric(TN_mg_N_per_L)) %>% 
#   dplyr::select(-c(Incubation_Water_Mass_g, Dry_Sediment_Mass_g, Final_Water_mass_g, Initial_Water_mass_g)) %>% 
#   group_by(Sample_ID, Rep) %>%
#   summarise(across(where(is.numeric),
#                    list(mean = ~mean(.x, na.rm = TRUE), 
#                         sd = ~sd(.x, na.rm = TRUE), 
#                         cv = ~(sd(.x, na.rm = TRUE)/mean(.x, na.rm = TRUE))*100)))
# 
# ggplot(cv, aes(x = TN_mg_N_per_L_cv)) + 
#   geom_histogram()

medians = all_data %>% 
  separate(Sample_Name, c("Sample_ID", "Rep"), sep = "-") %>% 
  mutate(Rep = if_else(grepl("D", Rep), "Dry", "Wet")) %>%
  mutate(NPOC_mg_C_per_L = as.numeric(NPOC_mg_C_per_L)) %>% 
  mutate(TN_mg_N_per_L = as.numeric(TN_mg_N_per_L)) %>% 
  dplyr::select(-c(Incubation_Water_Mass_g, Dry_Sediment_Mass_g, Final_Water_mass_g, Initial_Water_mass_g)) %>% 
  group_by(Sample_ID, Rep) %>%
  summarise(across(where(is.numeric), ~median(.x, na.rm = TRUE))) %>% 
  rename_with(.cols = c(SpC:Lost_Gravimetric_Water), .fn = ~ paste0("median_", .x)) 
  
write.csv(medians,"C:/GitHub/ECA_Multireactor_Incubations/Data/Cleaned Data/2024-05-29_Medians_ECA_Data.csv") 
  

effect_data <- medians %>% 
  #separate(Sample_Name, c("Sample_Name", "Treat"), sep = "-") %>% 
  filter(Sample_ID != "EC_011_INC") %>% 
  filter(Sample_ID != "EC_012_INC") %>% 
  relocate(median_Initial_Gravimetric_Moisture:median_Final_Gravimetric_Moisture, .after = median_Lost_Gravimetric_Water) %>% 
  #filter(Sample_Name != "EC_021_INC") %>% 
  group_by(Sample_ID) %>% 
  mutate(across(c(median_SpC:median_tn_percent), ~. [Rep == "Wet"] - .[Rep  == "Dry"])) %>% 
  rename_with(.cols = c(median_SpC:median_tn_percent), .fn = ~ paste0("diff_", .x)) %>% 
  distinct(Sample_ID, .keep_all = TRUE) %>% 
  rename(median_Dry_Initial_Gravimetric_Moisture = median_Initial_Gravimetric_Moisture ) %>% 
  rename(median_Dry_Final_Gravimetric_Moisture = median_Final_Gravimetric_Moisture ) %>% 
  rename(median_Dry_Lost_Gravimetric_Moisture = median_Lost_Gravimetric_Water)
  

write.csv(effect_data,"C:/GitHub/ECA_Multireactor_Incubations/Data/Cleaned Data/2024-05-29_Effect_Median_ECA_Data.csv") 
