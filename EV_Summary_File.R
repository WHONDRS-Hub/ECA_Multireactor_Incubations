## EV Data Package Summary File
#11/12024 M.Laan

library(tidyverse)

rm(list=ls());graphics.off()

# Set working directory to data file
#Example:
pnnl.user = 'laan208'

# Respiration -------------------------------------------------------------

all_respiration <- read.csv("C:/Users/laan208/OneDrive - PNNL/Shared Documents - Core Richland and Sequim Lab-Field Team/Data Generation and Files/ECA/INC/03_ProcessedData/EV_Sediment_Incubations_Respiration_Rates_ReadyForBoye_2024-11-01.csv") %>% 
  dplyr::select(c(Sample_Name, SpC, pH, Temp, Respiration_Rate_mg_DO_per_L_per_H, Respiration_Rate_mg_DO_per_kg_per_H, Methods_Deviation)) %>% 
  mutate(across(c(SpC:Respiration_Rate_mg_DO_per_kg_per_H), as.numeric))


median_respiration = all_respiration %>% 
  mutate(Respiration_Rate_mg_DO_per_L_per_H = ifelse(grepl("INC_Method_001|INC_007", Methods_Deviation), NA, Respiration_Rate_mg_DO_per_L_per_H)) %>% 
  #INC_Method_001 - missing replicates 
  #INC_007 - optode disk put on backwards
  mutate(Respiration_Rate_mg_DO_per_kg_per_H = ifelse(grepl("INC_Method_001|INC_007", Methods_Deviation), NA, Respiration_Rate_mg_DO_per_kg_per_H)) %>% 
  mutate(SpC = ifelse(grepl("INC_Method_001", Methods_Deviation), NA, SpC)) %>% 
  mutate(pH = ifelse(grepl("INC_Method_001|PH_000", Methods_Deviation), NA, pH)) %>% 
    #PH_000 - didn't take pH measurement
  mutate(Temp = ifelse(grepl("INC_Method_001", Methods_Deviation), NA, Temp)) %>%
  mutate(Respiration_Rate_mg_DO_per_kg_per_H = ifelse(Respiration_Rate_mg_DO_per_kg_per_H == "-9999", NA, Respiration_Rate_mg_DO_per_kg_per_H)) %>% 
  separate(Sample_Name, c("Sample_ID", "Rep"), sep = "-") %>% 
  mutate(Rep = if_else(grepl("D", Rep), "D", "W")) %>%
  group_by(Sample_ID, Rep) %>%
  summarise(across(where(is.numeric),
                   list(Median = ~median(.x, na.rm = TRUE),
                        cv = ~sd(.x, na.rm = TRUE)/mean(.x, na.rm =TRUE),
                        n = ~sum(!is.na(.x))), 
                   .names = "{.fn}_{.col}")) %>% 
  ungroup() %>% 
  group_by(Sample_ID, Rep) %>%
  mutate(Remove = ifelse(all(c_across(starts_with("n_")) == 5), "FALSE", "TRUE")) %>% ## Check CV's, then remove 
  dplyr::select(c(Sample_ID, Rep, Median_SpC, Median_pH, Median_Temp, Median_Respiration_Rate_mg_DO_per_L_per_H, Median_Respiration_Rate_mg_DO_per_kg_per_H, Remove)) %>% 
  # select(c(Sample_ID, Rep, Median_SpC_microsiemens_per_cm, Median_pH, Median_Temperature_degC, Median_Respiration_Rate_mg_DO_per_L_per_H, Median_Respiration_Rate_mg_DO_per_kg_per_H, Remove)) %>% 
  unite(Sample_Name, c("Sample_ID", "Rep"))

# Gravimetric Moisture ----------------------------------------------------

grav_inc = read.csv("C:/Users/laan208/OneDrive - PNNL/Shared Documents - Core Richland and Sequim Lab-Field Team/Data Generation and Files/ECA/INC/03_ProcessedData/EV_Drying_Masses_Summary_ReadyForBoye_05-07-2024.csv") %>% 
  mutate(across(c(Initial_Water_Mass_g:Incubation_Water_Mass_g), as.numeric))

median_grav = grav_inc %>% 
  mutate(Initial_Gravimetric_Moisture = ifelse(grepl("INC_Method_001", Methods_Deviation), NA, Initial_Gravimetric_Moisture)) %>% 
  mutate(Final_Gravimetric_Moisture = ifelse(grepl("INC_Method_001", Methods_Deviation), NA, Final_Gravimetric_Moisture)) %>% 
  #missing replicates
  separate(Sample_Name, c("Sample_ID", "Rep"), sep = "-") %>% 
  mutate(Rep = if_else(grepl("D", Rep), "D", "W")) %>%
  group_by(Sample_ID, Rep) %>%
  summarise(across(where(is.numeric),
                   list(Median = ~median(.x, na.rm = TRUE),
                        cv = ~sd(.x, na.rm = TRUE)/mean(.x, na.rm =TRUE),
                        n = ~sum(!is.na(.x))), 
                   .names = "{.fn}_{.col}")) %>% 
  ungroup() %>% 
  group_by(Sample_ID, Rep) %>%
  mutate(Remove = ifelse(all(c_across(starts_with("n_")) == 5), "FALSE", "TRUE")) %>% ## Check CV's, then remove 
  dplyr::select(c(Sample_ID, Rep, Median_Initial_Gravimetric_Moisture, Median_Final_Gravimetric_Moisture, Remove)) %>% 
  unite(Sample_Name, c("Sample_ID", "Rep"))

# Iron --------------------------------------------------------------------

all_iron <- read_csv("C:/Users/laan208/OneDrive - PNNL/Shared Documents - Core Richland and Sequim Lab-Field Team/Data Generation and Files/ECA/Fe/03_ProcessedData/EV_SFE_ReadyForBoye_2024-10-30.csv") %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "SFE", "INC")) %>% 
  rename(Dev_Fe = Methods_Deviation) %>% 
  filter(Dev_Fe != "SFE_ACID_001")

# Fe replicates have high CV's, but no apparent batch effects/deviations

median_iron = all_iron  %>% 
  left_join(grav_inc, by = "Sample_Name") %>% 
  mutate(Fe_mg_per_L = ifelse(grepl("INC_Method_001", Dev_Fe), NA, Fe_mg_per_L)) %>% 
  mutate(Fe_mg_per_kg = ifelse(grepl("INC_Method_001", Dev_Fe), NA, Fe_mg_per_kg)) %>% 
  mutate(Fe_mg_per_L = if_else(grepl("SFE_Below", Fe_mg_per_L), "0.001", Fe_mg_per_L)) %>% 
  mutate(Fe_mg_per_L = if_else(grepl("SFE_Above", Fe_mg_per_L), str_extract(Fe_mg_per_L, "(?<=\\|[^|]{1,100}\\|)\\d+\\.\\d+"), Fe_mg_per_L)) %>% 
  mutate(Fe_mg_per_L = as.numeric(Fe_mg_per_L)) %>%
  mutate(Fe_mg_per_kg = if_else(grepl("SFE_Above", Fe_mg_per_kg), str_extract(Fe_mg_per_kg, "(?<=\\|[^|]{1,100}\\|)\\d+\\.\\d+"), Fe_mg_per_kg)) %>% 
  mutate(Fe_mg_per_kg = if_else(grepl("SFE_Below", Fe_mg_per_kg), as.numeric(Fe_mg_per_L * (Incubation_Water_Mass_g/Dry_Sediment_Mass_g)), as.numeric(Fe_mg_per_kg))) %>%
  mutate(Fe_mg_per_kg = as.numeric(Fe_mg_per_kg)) %>% 
  separate(Sample_Name, c("Sample_ID", "Rep"), sep = "-") %>% 
  mutate(Rep = if_else(grepl("D", Rep), "D", "W")) %>%
  group_by(Sample_ID, Rep) %>%
  summarise(across(where(is.numeric),
                   list(Median = ~median(.x, na.rm = TRUE),
                        cv = ~sd(.x, na.rm = TRUE)/mean(.x, na.rm =TRUE),
                        n = ~sum(!is.na(.x))), 
                   .names = "{.fn}_{.col}")) %>% 
  ungroup() %>% 
  group_by(Sample_ID, Rep) %>%
  mutate(Remove = ifelse(all(c_across(starts_with("n_")) == 5), "FALSE", "TRUE")) %>% ## Check CV's, then remove 
  dplyr::select(c(Sample_ID, Rep, Median_Fe_mg_per_L, Median_Fe_mg_per_kg, Remove)) %>% 
  unite(Sample_Name, c("Sample_ID", "Rep"))


# Merge Median Summary File -----------------------------------------------

medians = left_join(median_respiration, median_grav, by = "Sample_Name") %>% 
  left_join(median_iron, by = "Sample_Name") %>% 
  mutate(remove_any_true = if_any(everything(), ~ .x == TRUE)) %>% 
  dplyr::select(-c(Remove.x, Remove.y)) %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "_INC_", "-")) %>% 
  rename(Median_Missing_Reps = remove_any_true) %>% 
  mutate(Median_Missing_Reps = if_else(is.na(Median_Missing_Reps), FALSE, Median_Missing_Reps))

write.csv(medians,"C:/Users/laan208/OneDrive - PNNL/Shared Documents - Core Richland and Sequim Lab-Field Team/Data Generation and Files/ECA/INC/03_ProcessedData/EV_Sediment_Incubations_Respiration_Rates_Summary_ReadyForBoye_2024-11-06.csv") 
