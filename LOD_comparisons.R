library(tidyverse)

orig = read.csv ("Z:/00_Cross-SFA_ESSDIVE-Data-Package-Upload/01_Study-Data-Package-Folders/ECA_Data_Package/EC_Data_Package/Sample_Data/EC_Sediment_SpC_pH_Temp_Respiration.csv", skip = 2) %>% 
  filter(grepl("EC", Sample_Name)) %>% 
  filter(Respiration_Rate_mg_DO_per_L_per_H != -9999) %>% 
  mutate_at(vars(SpC_microsiemens_per_cm:DO_Concentration_At_Incubation_Time_Zero), as.numeric) #%>% 
  #dplyr::select(c(Sample_Name, Respiration_Rate_mg_DO_per_L_per_H)) %>% 
  #rename(Respiration_Rate_mg_L_orig = Respiration_Rate_mg_DO_per_L_per_H)

final_lod = read.csv("C:/Users/laan208/OneDrive - PNNL/Shared Documents - Core Richland and Sequim Lab-Field Team/Data Generation and Files/ECA/INC/03_ProcessedData/EC_Sediment_Incubations_Respiration_Rates_ReadyForBoye_2024-07-31.csv")%>% 
  filter(Respiration_Rate_mg_DO_per_L_per_H != -9999) %>%
  dplyr::select(-c(Material, Methods_Deviation)) %>% 
  mutate_at(vars(SpC:DO_Concentration_At_Incubation_Time_Zero), as.numeric) %>% 
  rename_with(where(is.numeric), .fn = ~ paste0("final_", .x))

half_lod = read.csv("C:/Users/laan208/OneDrive - PNNL/Shared Documents - Core Richland and Sequim Lab-Field Team/Data Generation and Files/ECA/INC/03_ProcessedData/ECA_Sediment_Incubations_Respiration_Rates_Half_LOD_2024-07-30.csv") %>% 
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


compare = full_join(orig, final_lod, by = "Sample_Name") %>% 
  mutate(spc_comp = if_else(final_SpC == SpC_microsiemens_per_cm, TRUE, FALSE)) %>% 
  mutate(pH_comp = if_else(final_pH == pH, TRUE, FALSE)) %>% 
  mutate(temp_comp = if_else(final_Temp == Temperature_degC, TRUE, FALSE)) %>% 
  mutate(resp_mg_L_comp = if_else(final_Respiration_Rate_mg_DO_per_L_per_H == Respiration_Rate_mg_DO_per_L_per_H, TRUE, FALSE)) %>% 
  mutate(resp_mg_kg_comp = if_else(final_Respiration_Rate_mg_DO_per_kg_per_H == Respiration_Rate_mg_DO_per_kg_per_H, TRUE, FALSE)) %>% 
  mutate(r2_comp = if_else(final_Respiration_R_Squared == Respiration_R_Squared, TRUE, FALSE)) %>% 
  mutate(r2_adj_comp = if_else(final_Respiration_R_Squared_Adj == Respiration_R_Squared_Adj, TRUE, FALSE)) %>%
  mutate(p_comp = if_else(final_Respiration_p_Value == Respiration_p_value, TRUE, FALSE)) %>%
  mutate(total_comp = if_else(final_Total_Incubation_Time_Min == Total_Incubation_Time_Min, TRUE, FALSE)) %>%
  mutate(points_comp = if_else(final_Number_Points_In_Respiration_Regression == Number_Points_In_Respiration_Regression, TRUE, FALSE)) %>%
  mutate(rem_comp = if_else(final_Number_Points_Removed_Respiration_Regression == Number_Points_Removed_Respiration_Regression, TRUE, FALSE))

compare %>% 
  filter(Respiration_Rate_mg_DO_per_L_per_H > -90) %>% 
ggplot(aes(x = final_Respiration_Rate_mg_DO_per_L_per_H, y = Respiration_Rate_mg_DO_per_L_per_H)) +
  geom_point()+
  geom_abline()

compare %>% 
  filter(Respiration_Rate_mg_DO_per_kg_per_H > -9999) %>% 
  filter(Respiration_Rate_mg_DO_per_L_per_H > -90) %>% 
  ggplot(aes(x = final_Respiration_Rate_mg_DO_per_kg_per_H, y = Respiration_Rate_mg_DO_per_kg_per_H)) +
  geom_point()+
  geom_abline()

orig_median = read.csv ("Z:/00_Cross-SFA_ESSDIVE-Data-Package-Upload/01_Study-Data-Package-Folders/ECA_Data_Package/EC_Data_Package/Sample_Data/EC_Sediment_Sample_Data_Summary.csv", skip = 2) %>% 
  filter(grepl("EC", Sample_Name)) %>% 
  mutate_at(vars(Median_SpC_microsiemens_per_cm:Median_01397_N_percent_per_mg), as.numeric) %>%  #%>% 
  select(-c(Field_Name, IGSN, Material))

new_median = read.csv ("Z:/00_Cross-SFA_ESSDIVE-Data-Package-Upload/01_Study-Data-Package-Folders/ECA_Data_Package/ECA_EC_Summary_ReadyForBoye_2024-07-31.csv") %>% 
  select(-c(X)) %>% 
  rename_with(where(is.numeric), .fn = ~ paste0("final_", .x))

compare_median = full_join(orig_median, new_median, by = "Sample_Name") %>% 
  mutate(spc_comp = if_else(final_Median_SpC == Median_SpC_microsiemens_per_cm, TRUE, FALSE)) %>% 
  mutate(pH_comp = if_else(final_Median_pH == Median_pH, TRUE, FALSE)) %>% 
  mutate(temp_comp = if_else(final_Median_Temp == Median_Temperature_degC, TRUE, FALSE)) %>% 
  mutate(resp_mg_L_comp = if_else(final_Median_Respiration_Rate_mg_DO_per_L_per_H == Median_Respiration_Rate_mg_DO_per_L_per_H, TRUE, FALSE)) %>% 
  mutate(resp_mg_kg_comp = if_else(final_Median_Respiration_Rate_mg_DO_per_kg_per_H == Median_Respiration_Rate_mg_DO_per_kg_per_H, TRUE, FALSE))
 
compare_median %>% 
  #filter(Respiration_Rate_mg_DO_per_kg_per_H > -9999) %>% 
  filter(Median_Respiration_Rate_mg_DO_per_L_per_H > -90) %>% 
  ggplot(aes(x = final_Median_Respiration_Rate_mg_DO_per_kg_per_H, y = Median_Respiration_Rate_mg_DO_per_kg_per_H)) +
  geom_point()+
  geom_abline()

