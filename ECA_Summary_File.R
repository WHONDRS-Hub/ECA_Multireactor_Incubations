## ECA Data Package Summary File
# 5/29/2024 M.Laan

library(tidyverse)

rm(list=ls());graphics.off()

# Set working directory to data file
#Example:
pnnl.user = 'laan208'

#Set wd to Boye Files
setwd("Z:/00_Cross-SFA_ESSDIVE-Data-Package-Upload/01_Study-Data-Package-Folders/ECA_Data_Package/EC_Data_Package/Sample_Data/")


# Respiration -------------------------------------------------------------

all_respiration <- read.csv("EC_Sediment_SpC_pH_Temp_Respiration.csv", skip = 2) %>% 
  slice(-1:-11) %>% 
  filter(Field_Name != "#End_Data") %>% 
  dplyr::select(c(Sample_Name, Specific_Conductance_microsiemens_per_centimeter, pH, Temperature_degrees_Celsius, Respiration_Rate_mg_DO_per_L_per_H, Respiration_Rate_mg_DO_per_kg_per_H, Methods_Deviation)) %>% 
  mutate(across(c(Specific_Conductance_microsiemens_per_centimeter:Respiration_Rate_mg_DO_per_kg_per_H), as.numeric))


median_respiration = all_respiration %>% 
  mutate(Respiration_Rate_mg_DO_per_L_per_H = ifelse(grepl("INC_Method_001|INC_Method_002|INC_QA_004", Methods_Deviation), NA, Respiration_Rate_mg_DO_per_L_per_H)) %>% 
  #missing replicates (EC_072-W5/D5),  overexposed samples (EC_027, EC_013, EC_014), less sediment in sample (EC_012-D5)
  mutate(Respiration_Rate_mg_DO_per_kg_per_H = ifelse(grepl("INC_Method_001|INC_Method_002|INC_QA_004", Methods_Deviation), NA, Respiration_Rate_mg_DO_per_kg_per_H)) %>% 
  mutate(Specific_Conductance_microsiemens_per_centimeter = ifelse(grepl("INC_Method_001|INC_Method_002", Methods_Deviation), NA, Specific_Conductance_microsiemens_per_centimeter)) %>% 
  mutate(pH = ifelse(grepl("INC_Method_001|INC_Method_002", Methods_Deviation), NA, pH)) %>% 
  mutate(Temperature_degrees_Celsius = ifelse(grepl("INC_Method_001|INC_Method_002", Methods_Deviation), NA, Temperature_degrees_Celsius)) %>% 
  mutate(Respiration_Rate_mg_DO_per_kg_per_H = ifelse(Respiration_Rate_mg_DO_per_kg_per_H == "-9999", NA, Respiration_Rate_mg_DO_per_L_per_H)) %>% 
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
  select(c(Sample_ID, Rep, Median_Specific_Conductance_microsiemens_per_centimeter, Median_pH, Median_Temperature_degrees_Celsius, Median_Respiration_Rate_mg_DO_per_L_per_H, Median_Respiration_Rate_mg_DO_per_kg_per_H, Remove)) %>% 
  unite(Sample_Name, c("Sample_ID", "Rep"))

# Gravimetric Moisture ----------------------------------------------------

grav_inc = read.csv("Z:/00_Cross-SFA_ESSDIVE-Data-Package-Upload/01_Study-Data-Package-Folders/ECA_Data_Package/EC_Data_Package/Sample_Data/EC_Sediment_Gravimetric_Moisture.csv", skip = 2) %>% 
  slice(-1:-11) %>% 
  filter(Field_Name != "#End_Data") %>% 
  select(-c(Field_Name)) %>% 
  mutate(across(c(Initial_Water_Mass_g:Incubation_Water_Mass_g), as.numeric))

#Some dry reps have high CV: EC_057 (low moisture), EC_081 (one lower sample), EC_063 (low moisture), EC_088 (one lower sample), EC_076 (low moisture), EC_071 (one lower sample), EC_056 (low moisture), EC_069 (one slightly lower) 

median_grav = grav_inc %>% 
  mutate(X62948_Initial_Gravimetric_Moisture_g_per_g = ifelse(grepl("INC_Method_001|INC_Method_002", Methods_Deviation), NA, X62948_Initial_Gravimetric_Moisture_g_per_g)) %>% 
  mutate(X62948_Final_Gravimetric_Moisture_g_per_g = ifelse(grepl("INC_Method_001|INC_Method_002", Methods_Deviation), NA, X62948_Final_Gravimetric_Moisture_g_per_g)) %>% 
  #missing replicates (EC_072-W5/D5), less sediment in sample (EC_012-D5)
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
  select(c(Sample_ID, Rep, Median_X62948_Initial_Gravimetric_Moisture_g_per_g, Median_X62948_Final_Gravimetric_Moisture_g_per_g, Remove)) %>% 
  unite(Sample_Name, c("Sample_ID", "Rep"))

# Iron --------------------------------------------------------------------

all_iron <- read_csv("Z:/00_Cross-SFA_ESSDIVE-Data-Package-Upload/01_Study-Data-Package-Folders/ECA_Data_Package/EC_Data_Package/Sample_Data/EC_Sediment_Fe.csv", skip = 2) %>% 
  slice(-1:-11) %>% 
  filter(Field_Name != "#End_Data") %>% 
  select(-c(Field_Name)) %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "SFE", "INC")) %>% 
  rename(Dev_Fe = Methods_Deviation)

# Fe replicates have high CV's, but no apparent batch effects/deviations
  
median_iron = all_iron  %>% 
  left_join(grav_inc, by = "Sample_Name") %>% 
  mutate(Fe_mg_per_L = ifelse(grepl("INC_Method_001|INC_Method_002", Dev_Fe), NA, Fe_mg_per_L)) %>% 
  mutate(Fe_mg_per_kg = ifelse(grepl("INC_Method_001|INC_Method_002", Dev_Fe), NA, Fe_mg_per_kg)) %>% 
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
  select(c(Sample_ID, Rep, Median_Fe_mg_per_L, Median_Fe_mg_per_kg, Remove)) %>% 
  unite(Sample_Name, c("Sample_ID", "Rep"))


# ATP ---------------------------------------------------------------------
atp_all = read.csv("Z:/00_Cross-SFA_ESSDIVE-Data-Package-Upload/01_Study-Data-Package-Folders/ECA_Data_Package/EC_Data_Package/Sample_Data/EC_Sediment_ATP.csv", skip = 2) %>% 
  slice(-1:-11) %>% 
  filter(Field_Name != "#End_Data") %>% 
  select(-c(Field_Name)) %>% 
    mutate(Sample_Name  = str_replace(Sample_Name, "ATP", "INC")) %>% 
  mutate(across(c(ATP_nanomoles_per_L:ATP_picomoles_per_g), as.numeric))

#ATP has high CV's but no apparent Method Deviations
median_atp = atp_all %>% 
  mutate(ATP_nanomoles_per_L = ifelse(grepl("INC_Method_001", Methods_Deviation), NA, ATP_nanomoles_per_L)) %>% 
  # ATP_002 (don't have sample), INC_Method_001
  mutate(ATP_picomoles_per_g = ifelse(grepl("INC_Method_001", Methods_Deviation), NA, ATP_picomoles_per_g)) %>%
  mutate(ATP_nanomoles_per_L = ifelse(ATP_nanomoles_per_L == -9999, NA, ATP_nanomoles_per_L)) %>% 
  mutate(ATP_picomoles_per_g = ifelse(ATP_picomoles_per_g == -9999, NA, ATP_picomoles_per_g)) %>% 
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
  select(c(Sample_ID, Rep, Median_ATP_nanomoles_per_L, Median_ATP_picomoles_per_g, Remove)) %>% 
  unite(Sample_Name, c("Sample_ID", "Rep"))

# NPOC/TN -----------------------------------------------------------------

npoc_tn_all = read.csv("Z:/00_Cross-SFA_ESSDIVE-Data-Package-Upload/01_Study-Data-Package-Folders/ECA_Data_Package/EC_Data_Package/Sample_Data/EC_Sediment_NPOC_TN.csv", skip = 2) %>% 
  slice(-1:-11) %>% 
  filter(Field_Name != "#End_Data") %>% 
  select(-c(Field_Name)) %>% 
    mutate(Sample_Name  = str_replace(Sample_Name, "SIR", "INC")) %>% 
  mutate(across(c(Extractable_NPOC_mg_per_L:Extractable_TN_mg_per_kg), as.numeric))

#VI_OCN_010 - missing reps, VB_OCN_001 - broken vials 

#remove blanks from median?

median_npoc_tn = npoc_tn_all %>% 
   mutate(Extractable_NPOC_mg_per_L = ifelse(grepl("VI_OCN_010|VB_OCN_001", Methods_Deviation), NA, Extractable_NPOC_mg_per_L)) %>% 
  mutate(Extractable_NPOC_mg_per_kg = ifelse(grepl("VI_OCN_010|VB_OCN_001", Methods_Deviation), NA, Extractable_NPOC_mg_per_kg)) %>% 
  mutate(Extractable_TN_mg_per_L = ifelse(grepl("VI_OCN_010|VB_OCN_001", Methods_Deviation), NA, Extractable_TN_mg_per_L)) %>% 
  mutate(Extractable_TN_mg_per_kg = ifelse(grepl("VI_OCN_010|VB_OCN_001", Methods_Deviation), NA, Extractable_TN_mg_per_kg)) %>% 
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
  select(c(Sample_ID, Rep, Median_Extractable_NPOC_mg_per_L, Median_Extractable_NPOC_mg_per_kg, Median_Extractable_TN_mg_per_L, Median_Extractable_TN_mg_per_kg, Remove)) %>% 
  unite(Sample_Name, c("Sample_ID", "Rep")) 

# C/N ---------------------------------------------------------------------

cn_all = read.csv("Z:/00_Cross-SFA_ESSDIVE-Data-Package-Upload/01_Study-Data-Package-Folders/ECA_Data_Package/EC_Data_Package/Sample_Data/EC_Sediment_CN.csv", skip = 2) %>% 
  slice(-1:-11) %>% 
  filter(Field_Name != "#End_Data") %>% 
  select(-c(Field_Name)) %>% 
  mutate(Sample_Name  = str_replace(Sample_Name, "SCN", "INC")) %>% 
  mutate(across(c(X01395_C_percent_per_mg :X01397_N_percent_per_mg), as.numeric)) 
  
median_cn = cn_all %>% 
  mutate(X01395_C_percent_per_mg = ifelse(grepl("INC_Method_001", Methods_Deviation), NA, X01395_C_percent_per_mg)) %>% 
  mutate(X01397_N_percent_per_mg = ifelse(grepl("INC_Method_001", Methods_Deviation), NA, X01397_N_percent_per_mg)) %>%
  mutate(X01395_C_percent_per_mg = ifelse(X01395_C_percent_per_mg == -9999, NA, X01395_C_percent_per_mg)) %>% 
  mutate(X01397_N_percent_per_mg = ifelse(X01397_N_percent_per_mg == -9999, NA, X01397_N_percent_per_mg)) %>% 
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
  select(c(Sample_ID, Rep, Median_X01395_C_percent_per_mg, Median_X01397_N_percent_per_mg, Remove)) %>% 
  unite(Sample_Name, c("Sample_ID", "Rep"))

# Ions --------------------------------------------------------------------

# Waiting for ReadyForBoye File from VGC


# Merge Median Summary File -----------------------------------------------

medians = left_join(median_respiration, median_grav, by = "Sample_Name") %>% 
  left_join(median_iron, by = "Sample_Name") %>% 
  left_join(median_atp, by = "Sample_Name") %>% 
  left_join(median_npoc_tn, by = "Sample_Name") %>% 
  left_join(median_cn, by = "Sample_Name") %>% 
  mutate(remove_any_true = if_any(everything(), ~ .x == TRUE)) %>% 
  select(-c(Remove.x, Remove.y, Remove.x.x, Remove.y.y, Remove.x.x.x, Remove.y.y.y)) %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "_INC_", "-")) %>% 
  rename(Median_Missing_Reps = remove_any_true) %>% 
  mutate(Median_Missing_Reps = if_else(is.na(Median_Missing_Reps), FALSE, Median_Missing_Reps))
  
write.csv(medians,"Z:/00_Cross-SFA_ESSDIVE-Data-Package-Upload/01_Study-Data-Package-Folders/ECA_Data_Package/EC_Data_Package/ECA_EC_Summary.csv") 

# Effect Size -------------------------------------------------------------

# Remove Gravimetric Moisture from this?

effect_data <- medians %>% 
  separate(Sample_Name, c("Sample_Name", "Treat"), sep = "-") %>% 
  mutate(Sample_Name = paste0(Sample_Name, "_all")) %>% 
  group_by(Sample_Name) %>% 
  mutate(across(where(is.numeric), ~. [Treat == "W"] - .[Treat  == "D"])) %>% 
  rename_with(~ str_replace_all(., "Median_", "Effect_Size_")) %>% 
  distinct(Sample_Name, .keep_all = TRUE) %>% 
  mutate(Methods_Deviation = if_else(Effect_Size_Missing_Reps == TRUE, "EFFECT_001", "N/A")) %>% 
  select(-c(Treat, Effect_Size_Missing_Reps))
  
write.csv(effect_data,"Z:/00_Cross-SFA_ESSDIVE-Data-Package-Upload/01_Study-Data-Package-Folders/ECA_Data_Package/EC_Data_Package/ECA_EC_Effect_Size.csv") 
