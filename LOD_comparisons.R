library(tidyverse)

orig = read.csv ("Z:/00_Cross-SFA_ESSDIVE-Data-Package-Upload/01_Study-Data-Package-Folders/ECA_Data_Package/EC_Data_Package/Sample_Data/EC_Sediment_SpC_pH_Temp_Respiration.csv", skip = 2) %>% 
  filter(grepl("EC", Sample_Name)) %>% 
  filter(Respiration_Rate_mg_DO_per_L_per_H != -9999) %>% 
  dplyr::select(c(Sample_Name, Respiration_Rate_mg_DO_per_L_per_H)) %>% 
  rename(Respiration_Rate_mg_L_orig = Respiration_Rate_mg_DO_per_L_per_H)

half_lod = read.csv("C:/Users/laan208/OneDrive - PNNL/Shared Documents - Core Richland and Sequim Lab-Field Team/Data Generation and Files/ECA/INC/03_ProcessedData/ECA_Sediment_Incubations_Respiration_Rates_Half_LOD_2024-07-30.csv")%>% 
  filter(Respiration_Rate_mg_DO_per_L_per_H != -9999) %>% 
  dplyr::select(c(Sample_Name, Respiration_Rate_mg_DO_per_L_per_H))%>% 
  rename(Respiration_Rate_mg_L_half = Respiration_Rate_mg_DO_per_L_per_H)

full_lod  = read.csv("C:/Users/laan208/OneDrive - PNNL/Shared Documents - Core Richland and Sequim Lab-Field Team/Data Generation and Files/ECA/INC/03_ProcessedData/ECA_Sediment_Incubations_Respiration_Rates_Full_LOD_2024-07-30.csv")%>% 
  filter(Respiration_Rate_mg_DO_per_L_per_H != -9999) %>% 
  dplyr::select(c(Sample_Name, Respiration_Rate_mg_DO_per_L_per_H))%>% 
  rename(Respiration_Rate_mg_L_full = Respiration_Rate_mg_DO_per_L_per_H)

all_data = full_join(orig, half_lod, by = "Sample_Name") %>% 
  full_join(full_lod, by = "Sample_Name") %>% 
  mutate(across(c(Respiration_Rate_mg_L_orig:Respiration_Rate_mg_L_full), as.numeric)) %>% 
  mutate(diff_full = Respiration_Rate_mg_L_orig - Respiration_Rate_mg_L_full) %>% 
  mutate(percent_diff_full = abs(diff_full/Respiration_Rate_mg_L_orig)*100)%>% 
  mutate(diff_half = Respiration_Rate_mg_L_orig - Respiration_Rate_mg_L_half) %>% 
  mutate(percent_diff_half = abs(diff_half/Respiration_Rate_mg_L_orig)*100)


half_v_orig = all_data %>% filter(Respiration_Rate_mg_L_orig < -50) %>% 
  ggplot() + 
  geom_point(aes(x = Respiration_Rate_mg_L_orig, y = Respiration_Rate_mg_L_half)) +
  geom_abline()

half_v_full = all_data %>% filter(Respiration_Rate_mg_L_orig < -50) %>% 
  ggplot() + 
  geom_point(aes(x = Respiration_Rate_mg_L_full, y = Respiration_Rate_mg_L_half)) +
  geom_abline()

full_v_orig = all_data %>% filter(Respiration_Rate_mg_L_orig < -50) %>% 
  ggplot() + 
  geom_point(aes(x = Respiration_Rate_mg_L_orig, y = Respiration_Rate_mg_L_full)) +
  geom_abline()

ggarrange(half_v_orig, half_v_full, full_v_orig)



