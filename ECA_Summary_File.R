## ECA Data Package

library(tidyverse)

# 5/29/2024 M.Laan

rm(list=ls());graphics.off()

# Set working directory to data file
#Example:
pnnl.user = 'laan208'

# choose file dates to read in 

respiration.date = '2024-05-29'
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


# Respiration -------------------------------------------------------------

all_respiration <- read.csv(paste0("INC/03_ProcessedData/ECA_Sediment_Incubations_Respiration_Rates_ReadyForBoye_",respiration.date,".csv")) %>% 
  dplyr::select(c(Sample_Name, SpC, pH, Temp, Respiration_Rate_mg_DO_per_L_per_H, Respiration_Rate_mg_DO_per_kg_per_H, Methods_Deviation)) 


median_respiration = all_respiration %>% 
  mutate(Respiration_Rate_mg_DO_per_L_per_H = ifelse(grepl("INC_Method_001|INC_Method_002|INC_QA_004", Methods_Deviation), NA, Respiration_Rate_mg_DO_per_L_per_H)) %>% 
  #missing replicates (EC_072-W5/D5),  overexposed samples (EC_027, EC_013, EC_014), less sediment in sample (EC_012-D5)
  mutate(Respiration_Rate_mg_DO_per_kg_per_H = ifelse(grepl("INC_Method_001|INC_Method_002|INC_QA_004", Methods_Deviation), NA, Respiration_Rate_mg_DO_per_kg_per_H)) %>% 
  mutate(SpC = ifelse(grepl("INC_Method_001|INC_Method_002", Methods_Deviation), NA, SpC)) %>% 
  mutate(pH = ifelse(grepl("INC_Method_001|INC_Method_002", Methods_Deviation), NA, pH)) %>% 
  mutate(Temp = ifelse(grepl("INC_Method_001|INC_Method_002", Methods_Deviation), NA, Temp)) %>% 
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
  mutate(Remove = ifelse(all(c_across(ends_with("_n")) == 5), "FALSE", "TRUE")) %>% ## Check CV's, then remove 
  select(c(Sample_ID, Rep, Median_SpC, Median_pH, Median_Temp, Median_Respiration_Rate_mg_DO_per_L_per_H, Median_Respiration_Rate_mg_DO_per_kg_per_H, Remove)) %>% 
  unite(Sample_Name, c("Sample_ID", "Rep"))

# Gravimetric Moisture ----------------------------------------------------

grav_inc = read.csv(paste0("INC/03_ProcessedData/EC_Drying_Masses_Summary_ReadyForBoye_on_",grav.summary,".csv")) 

#Some dry reps have high CV: EC_057 (low moisture), EC_081 (one lower sample), EC_063 (low moisture), EC_088 (one lower sample), EC_076 (low moisture), EC_071 (one lower sample), EC_056 (low moisture), EC_069 (one slightly lower) 

median_grav = grav_inc %>% 
  mutate(Initial_Gravimetric_Moisture = ifelse(grepl("INC_Method_001|INC_Method_002", Methods_Deviation), NA, Initial_Gravimetric_Moisture)) %>% 
  mutate(Final_Gravimetric_Moisture = ifelse(grepl("INC_Method_001|INC_Method_002", Methods_Deviation), NA, Final_Gravimetric_Moisture)) %>% 
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
  mutate(Remove = ifelse(all(c_across(ends_with("_n")) == 5), "FALSE", "TRUE")) %>% ## Check CV's, then remove 
  select(c(Sample_ID, Rep,Median_Initial_Gravimetric_Moisture, Median_Final_Gravimetric_Moisture, Remove)) %>% 
  unite(Sample_Name, c("Sample_ID", "Rep"))

# Iron --------------------------------------------------------------------

all_iron <- read_csv(paste0("Fe/03_ProcessedData/EC_ReadyForBoye_",fe.date,".csv")) %>% 
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
  mutate(Remove = ifelse(all(c_across(ends_with("_n")) == 5), "FALSE", "TRUE")) %>% ## Check CV's, then remove 
  select(c(Sample_ID, Rep, Median_Fe_mg_per_L, Median_Fe_mg_per_kg, Remove)) %>% 
  unite(Sample_Name, c("Sample_ID", "Rep"))


# ATP ---------------------------------------------------------------------
atp_all = read.csv(paste0("ATP/03_ProcessedData/EC_ATP_ReadyForBoye_",atp.date,".csv")) %>% 
    mutate(Sample_Name  = str_replace(Sample_Name, "ATP", "INC")) 

#ATP has high CV's but no apparent Method Deviations

median_atp = atp_all %>% 
  mutate(ATP_nanomol_per_L = ifelse(grepl("INC_Method_001", Methods_Deviation), NA, ATP_nanomol_per_L)) %>% 
  # ATP_002 (don't have sample), INC_Method_001
  mutate(ATP_picomol_per_g = ifelse(grepl("INC_Method_001", Methods_Deviation), NA, ATP_picomol_per_g)) %>%
  mutate(ATP_nanomol_per_L = ifelse(ATP_nanomol_per_L == -9999, NA, ATP_nanomol_per_L)) %>% 
  mutate(ATP_picomol_per_g = ifelse(ATP_picomol_per_g == -9999, NA, ATP_picomol_per_g)) %>% 
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
  mutate(Remove = ifelse(all(c_across(ends_with("_n")) == 5), "FALSE", "TRUE")) %>% ## Check CV's, then remove 
  select(c(Sample_ID, Rep, Median_ATP_nanomol_per_L, Median_ATP_picomol_per_g, Remove)) %>% 
  unite(Sample_Name, c("Sample_ID", "Rep"))

# NPOC/TN -----------------------------------------------------------------

npoc_tn_all = read.csv(paste0("Boye_Files/EC/EC_NPOC_TN_Check_for_Duplicates_",npoc.tn.date,"_by_laan208.csv")) %>% 
  mutate(Sample_Name  = str_replace(Sample_ID, "SIR", "INC")) %>% 
    dplyr::select(-c(Date_of_Run, #Methods_Deviation, 
                   Method_Notes, duplicate, Sample_ID)) %>% 
  relocate(Sample_Name, .before = NPOC_mg_C_per_L)

#VI_OCN_010 - missing reps, VB_OCN_001 - broken vials 

#remove blanks from median?

median_npoc_tn = npoc_tn_all %>% 
  filter(!grepl("Blk", Sample_Name))  %>% 
  mutate(NPOC_mg_C_per_L = ifelse(grepl("VI_OCN_010|VB_OCN_001", Methods_Deviation), NA, NPOC_mg_C_per_L)) %>% 
  mutate(NPOC_mg_C_per_kg = ifelse(grepl("VI_OCN_010|VB_OCN_001", Methods_Deviation), NA, NPOC_mg_C_per_kg)) %>% 
  mutate(TN_mg_N_per_L = ifelse(grepl("VI_OCN_010|VB_OCN_001", Methods_Deviation), NA, TN_mg_N_per_L)) %>% 
  mutate(TN_mg_N_per_kg = ifelse(grepl("VI_OCN_010|VB_OCN_001", Methods_Deviation), NA, TN_mg_N_per_kg)) %>% 
  separate(Sample_Name, c("Sample_ID", "Rep"), sep = "-") %>% 
  mutate(Rep = if_else(grepl("D", Rep), "D", "W")) %>%
  mutate(NPOC_mg_C_per_L = as.numeric(NPOC_mg_C_per_L)) %>% 
  mutate(TN_mg_N_per_L = as.numeric(TN_mg_N_per_L)) %>% 
  group_by(Sample_ID, Rep) %>%
  summarise(across(where(is.numeric),
                   list(Median = ~median(.x, na.rm = TRUE),
                        cv = ~sd(.x, na.rm = TRUE)/mean(.x, na.rm =TRUE),
                        n = ~sum(!is.na(.x))), 
                   .names = "{.fn}_{.col}")) %>% 
  ungroup() %>% 
  group_by(Sample_ID, Rep) %>%
  mutate(Remove = ifelse(all(c_across(ends_with("_n")) == 5), "FALSE", "TRUE")) %>% ## Check CV's, then remove 
  select(c(Sample_ID, Rep, Median_NPOC_mg_C_per_L, Median_NPOC_mg_C_per_kg, Median_TN_mg_N_per_L, Median_TN_mg_N_per_kg, Remove)) %>% 
  unite(Sample_Name, c("Sample_ID", "Rep")) 



# C/N ---------------------------------------------------------------------

# waiting for ReadyForBoye File from Sophia

cn_all = read.csv(paste0("CN/02_FormattedData/ECA_CN_",cn.date,".csv")) %>% 
  mutate(Sample_Name  = str_replace(sample_id, "SCN", "INC")) %>% 
  separate(Sample_Name, c("Parent_ID", "Rep"), remove = TRUE, sep = "_INC") %>% 
  unite(Sample_Name, c("parent_id", "Rep"), sep = "_INC") %>% 
  dplyr::select(-c(sample_id, Parent_ID))


# Ions --------------------------------------------------------------------

# Waiting for ReadyForBoye File from VGC


# Merge Median Summary File -----------------------------------------------

medians = left_join(median_respiration, median_grav, by = "Sample_Name") %>% 
  left_join(median_iron, by = "Sample_Name") %>% 
  left_join(median_atp, by = "Sample_Name") %>% 
  left_join(median_npoc_tn, by = "Sample_Name") %>% 
  mutate(remove_any_true = if_any(everything(), ~ .x == TRUE)) %>% 
  select(-c(Remove.x, Remove.y, Remove.x.x, Remove.y.y, Remove)) %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "_INC_", "-"))
  


# Effect Size -------------------------------------------------------------

# Remove Gravimetric Moisture from this?

effect_data <- medians %>% 
  separate(Sample_Name, c("Sample_Name", "Treat"), sep = "-") %>% 
  group_by(Sample_Name) %>% 
  mutate(across(where(is.numeric), ~. [Treat == "W"] - .[Treat  == "D"])) %>% 
  rename_with(~ str_replace_all(., "Median_", "Effect_Size_")) %>% 
  distinct(Sample_Name, .keep_all = TRUE) %>% 
  select(-c(Treat))
  


write.csv(effect_data,"C:/GitHub/ECA_Multireactor_Incubations/Data/Cleaned Data/2024-05-29_Effect_Median_ECA_Data.csv") 
