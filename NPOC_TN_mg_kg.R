library(tidyverse)

npoc = read.csv("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Boye_Files/EC/EC_NPOC_TN_Check_for_Duplicates_2024-02-19_by_forb086.csv")

eca.moi.path = ("ECA/INC/03_ProcessedData/ECA_Drying_Masses_Summary_ReadyForBoye_01-29-2024.csv")

moisture = read.csv(paste0("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/", eca.moi.path)) %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "INC", "SIR")) %>% 
  rename(Sample_ID = Sample_Name)

npoc_kg = left_join(npoc, moisture, by = "Sample_ID") %>% 
  filter(!grepl("Blk", Sample_ID)) %>% 
  filter(!is.na(NPOC_mg_C_per_L)) %>% 
  filter(Incubation_Water_Mass_g > 0) %>% 
  mutate(NPOC_mg_C_per_L = as.numeric(NPOC_mg_C_per_L)) %>% 
  mutate(TN_mg_N_per_L = as.numeric(TN_mg_N_per_L)) %>% 
  mutate(NPOC_mg_C_per_kg = round(NPOC_mg_C_per_L * (Incubation_Water_Mass_g/Dry_Sediment_Mass_g), 2)) %>% 
  mutate(TN_mg_N_per_kg = round(TN_mg_N_per_L * (Incubation_Water_Mass_g/Dry_Sediment_Mass_g), 2)) %>% 
  select(c(Date_of_Run, Sample_ID, NPOC_mg_C_per_kg, TN_mg_N_per_kg))


npoc_tn_kg = left_join(npoc, npoc_kg, by = c("Date_of_Run", "Sample_ID")) %>% 
  relocate(NPOC_mg_C_per_kg, .after = TN_mg_N_per_L) %>% 
  relocate(TN_mg_N_per_kg, .after = NPOC_mg_C_per_kg)
  
write.csv(npoc_tn_kg, "C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Boye_Files/EC/EC_NPOC_TN_mg_per_kg_2024-03-01_by_laan208.csv", row.names = FALSE)
                       