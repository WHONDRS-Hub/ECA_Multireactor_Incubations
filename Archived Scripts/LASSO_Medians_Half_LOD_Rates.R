#### Sensitivity Analysis For ECA removals ####

# This script makes figures for ECA physical manuscript and performs PCA Analysis and LASSO regression

library(tidyverse)
library(corrplot)
library(ggpubr)
library(ggpmisc)
library(factoextra)
library(stringr)
library(glmnet)
library(magick)

rm(list=ls());graphics.off()

## Set image export

print.images = F

# Functions ---------------------------------------------------------------

# Transformation for normalization is cube root - have to cube root then add sign back to value to make it positive or negative
cube_root <- function(x) sign(x) * (abs(x))^(1/3)

# Read in/Merge Data ------------------------------------------------------------

## Individual Rate data for histograms ####


## Respiration ####
# change 0s to half of minimum rate
all_data = read.csv("C:/GitHub/ECA_Multireactor_Incubations/Data/EC_Sediment_SpC_pH_Temp_Respiration.csv", skip = 2) %>% 
  filter(grepl("EC", Sample_Name)) %>% 
  filter(!grepl("EC_011|EC_012|EC_023|EC_052|EC_053|EC_057", Sample_Name)) %>%   # remove samples with too much water (EC_011, EC_012), sample with no mg/kg (EC_023), duplicated NEON sites (EC_052, EC_053, EC_057)
  filter(Respiration_Rate_mg_DO_per_kg_per_H != -9999) %>% 
  mutate(Respiration_Rate_mg_DO_per_L_per_H = ifelse(Respiration_Rate_mg_DO_per_L_per_H == 0, -0.0175, Respiration_Rate_mg_DO_per_L_per_H)) %>% 
  mutate(Respiration_Rate_mg_DO_per_kg_per_H = ifelse(Respiration_Rate_mg_DO_per_kg_per_H == 0, -0.056, Respiration_Rate_mg_DO_per_kg_per_H)) %>% 
  mutate(across(c(SpC_microsiemens_per_cm:Respiration_Rate_mg_DO_per_kg_per_H), as.numeric))

median_chem = all_data %>% 
  select(c(Sample_Name, SpC_microsiemens_per_cm, pH, Temperature_degC, Methods_Deviation)) %>%
  mutate(SpC_microsiemens_per_cm = ifelse(grepl("INC_Method_001|INC_Method_002", Methods_Deviation), NA, SpC_microsiemens_per_cm)) %>% 
  mutate(pH = ifelse(grepl("INC_Method_001|INC_Method_002", Methods_Deviation), NA, pH)) %>% 
  mutate(Temperature_degC = ifelse(grepl("INC_Method_001|INC_Method_002", Methods_Deviation), NA, Temperature_degC)) %>%
  separate(Sample_Name, c("Sample_ID", "Rep"), sep = "-") %>% 
  mutate(Rep = if_else(grepl("D", Rep), "D", "W")) %>%
  group_by(Sample_ID) %>%
  summarise(across(where(is.numeric),
                   list(Median = ~median(.x, na.rm = TRUE),
                        cv = ~sd(.x, na.rm = TRUE)/mean(.x, na.rm =TRUE),
                        n = ~sum(!is.na(.x))), 
                   .names = "{.fn}_{.col}")) %>% 
  ungroup() %>% 
  group_by(Sample_ID) %>%
  mutate(Remove = ifelse(all(c_across(starts_with("n_")) == 5), "FALSE", "TRUE")) %>% ## Check CV's, then remove 
  select(c(Sample_ID, Median_SpC_microsiemens_per_cm, Median_pH, Median_Temperature_degC)) %>% 
  ungroup()
  
median_respiration = all_data %>% 
  select(c(Sample_Name, Respiration_Rate_mg_DO_per_L_per_H, Respiration_Rate_mg_DO_per_kg_per_H, Methods_Deviation)) %>% 
  mutate(Respiration_Rate_mg_DO_per_L_per_H = ifelse(grepl("INC_Method_001|INC_Method_002|INC_QA_004", Methods_Deviation), NA, Respiration_Rate_mg_DO_per_L_per_H)) %>% 
  #missing replicates (EC_072-W5/D5),  overexposed samples (EC_027, EC_013, EC_014), less sediment in sample (EC_012-D5)
  mutate(Respiration_Rate_mg_DO_per_kg_per_H = ifelse(grepl("INC_Method_001|INC_Method_002|INC_QA_004", Methods_Deviation), NA, Respiration_Rate_mg_DO_per_kg_per_H)) %>% 
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
   select(c(Sample_ID, Rep, Median_Respiration_Rate_mg_DO_per_L_per_H, Median_Respiration_Rate_mg_DO_per_kg_per_H)) %>% 
  unite(Sample_Name, c("Sample_ID", "Rep"), sep = "-", remove = FALSE) %>% 
  ungroup()

effect_respiration = median_respiration %>% 
  separate(Sample_Name, c("Sample_Name", "Treat"), sep = "-") %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "INC", "all")) %>% 
  select(c(Sample_Name, Treat, Median_Respiration_Rate_mg_DO_per_L_per_H, Median_Respiration_Rate_mg_DO_per_kg_per_H)) %>% 
  group_by(Sample_Name) %>% 
  mutate(across(where(is.numeric), ~. [Treat == "W"] - .[Treat  == "D"])) %>% 
  rename_with(~ str_replace_all(., "Median_", "Effect_Size_")) %>% 
  distinct(Sample_Name, .keep_all = TRUE) %>% 
  select(-c(Treat)) %>% 
  ungroup()

# cube_respiration = all_data %>% 
#   select(c(Sample_Name, Respiration_Rate_mg_DO_per_kg_per_H)) %>% 
#   filter(Respiration_Rate_mg_DO_per_kg_per_H != -9999) %>% 
#   mutate(cube_Respiration_mg_kg = cube_root(abs(as.numeric(Respiration_Rate_mg_DO_per_kg_per_H)))) %>% 
#   mutate(Treat = if_else(grepl("D", Sample_Name), "Dry", "Wet"))

## Calculate wet and dry medians ####

# ATP ####

atp = read.csv("C:/GitHub/ECA_Multireactor_Incubations/Data/EC_Sediment_ATP.csv", skip = 2) %>% 
  filter(grepl("EC", Sample_Name)) %>% 
  filter(!grepl("EC_011|EC_012|EC_023|EC_052|EC_053|EC_057", Sample_Name)) %>%  
  select(-c(Field_Name, IGSN, Material)) %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "ATP", "INC"))

median_atp = atp %>% 
  mutate(ATP_nanomoles_per_L = ifelse(grepl("INC_Method_001", Methods_Deviation), NA, ATP_nanomoles_per_L)) %>% 
  # ATP_002 (don't have sample), INC_Method_001
  mutate(ATP_picomoles_per_g = ifelse(grepl("INC_Method_001", Methods_Deviation), NA, ATP_picomoles_per_g)) %>%
  mutate(ATP_nanomoles_per_L = ifelse(ATP_nanomoles_per_L == -9999, NA, ATP_nanomoles_per_L)) %>% 
  mutate(ATP_picomoles_per_g = ifelse(ATP_picomoles_per_g == -9999, NA, ATP_picomoles_per_g)) %>% 
  mutate(across(c(ATP_nanomoles_per_L:ATP_picomoles_per_g), as.numeric)) %>% 
  separate(Sample_Name, c("Sample_ID", "Rep"), sep = "-") %>% 
  mutate(Rep = if_else(grepl("D", Rep), "D", "W")) %>%
  group_by(Sample_ID) %>%
  summarise(across(where(is.numeric),
                   list(Median = ~median(.x, na.rm = TRUE),
                        cv = ~sd(.x, na.rm = TRUE)/mean(.x, na.rm =TRUE),
                        n = ~sum(!is.na(.x))), 
                   .names = "{.fn}_{.col}")) %>% 
  ungroup() %>% 
  group_by(Sample_ID) %>%
  mutate(Remove = ifelse(all(c_across(starts_with("n_")) == 10), "FALSE", "TRUE")) %>% ## Check CV's, then remove 
  select(c(Sample_ID, Median_ATP_nanomoles_per_L, Median_ATP_picomoles_per_g)) %>% 
  ungroup()

## CN ####

cn = read.csv("C:/GitHub/ECA_Multireactor_Incubations/Data/EC_Sediment_CN.csv", skip = 2) %>% 
  filter(grepl("EC", Sample_Name)) %>% 
  filter(!grepl("EC_011|EC_012|EC_023|EC_052|EC_053|EC_057", Sample_Name)) %>%  
  select(-c(Field_Name, IGSN, Material)) %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "SCN", "INC"))

median_cn = cn %>% 
  mutate(X01395_C_percent_per_mg = ifelse(grepl("INC_Method_001", Methods_Deviation), NA, X01395_C_percent_per_mg)) %>% 
  mutate(X01397_N_percent_per_mg = ifelse(grepl("INC_Method_001", Methods_Deviation), NA, X01397_N_percent_per_mg)) %>%
  mutate(X01395_C_percent_per_mg = ifelse(X01395_C_percent_per_mg == -9999, NA, X01395_C_percent_per_mg)) %>% 
  mutate(X01397_N_percent_per_mg = ifelse(X01397_N_percent_per_mg == -9999, NA, X01397_N_percent_per_mg)) %>% 
  mutate(across(c(X01395_C_percent_per_mg:X01397_N_percent_per_mg), as.numeric)) %>% 
  separate(Sample_Name, c("Sample_ID", "Rep"), sep = "-") %>% 
  mutate(Rep = if_else(grepl("D", Rep), "D", "W")) %>%
  group_by(Sample_ID) %>%
  summarise(across(where(is.numeric),
                   list(Median = ~median(.x, na.rm = TRUE),
                        cv = ~sd(.x, na.rm = TRUE)/mean(.x, na.rm =TRUE),
                        n = ~sum(!is.na(.x))), 
                   .names = "{.fn}_{.col}")) %>% 
  ungroup() %>% 
  group_by(Sample_ID) %>%
  mutate(Remove = ifelse(all(c_across(starts_with("n_")) == 10), "FALSE", "TRUE")) %>% ## Check CV's, then remove 
  select(c(Sample_ID, Median_X01395_C_percent_per_mg, Median_X01397_N_percent_per_mg))%>% 
  ungroup()

## NPOC/TN ####

npoc_tn = read.csv("C:/GitHub/ECA_Multireactor_Incubations/Data/EC_Sediment_NPOC_TN.csv", skip = 2) %>% 
  filter(grepl("EC", Sample_Name)) %>% 
  filter(!grepl("EC_011|EC_012|EC_023|EC_052|EC_053|EC_057", Sample_Name)) %>%  
  select(-c(Field_Name, IGSN, Material)) %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "SIR", "INC"))

median_npoc_tn = npoc_tn %>% 
  mutate(Extractable_NPOC_mg_per_L = ifelse(grepl("VI_OCN_010|VB_OCN_001", Methods_Deviation), NA, Extractable_NPOC_mg_per_L)) %>% 
  mutate(Extractable_NPOC_mg_per_kg = ifelse(grepl("VI_OCN_010|VB_OCN_001", Methods_Deviation), NA, Extractable_NPOC_mg_per_kg)) %>% 
  mutate(Extractable_TN_mg_per_L = ifelse(grepl("VI_OCN_010|VB_OCN_001", Methods_Deviation), NA, Extractable_TN_mg_per_L)) %>% 
  mutate(Extractable_TN_mg_per_kg = ifelse(grepl("VI_OCN_010|VB_OCN_001", Methods_Deviation), NA, Extractable_TN_mg_per_kg)) %>% 
  mutate(across(c(Extractable_NPOC_mg_per_L:Extractable_TN_mg_per_kg), as.numeric)) %>% 
  separate(Sample_Name, c("Sample_ID", "Rep"), sep = "-") %>% 
  mutate(Rep = if_else(grepl("D", Rep), "D", "W")) %>%
  group_by(Sample_ID) %>%
  summarise(across(where(is.numeric),
                   list(Median = ~median(.x, na.rm = TRUE),
                        cv = ~sd(.x, na.rm = TRUE)/mean(.x, na.rm =TRUE),
                        n = ~sum(!is.na(.x))), 
                   .names = "{.fn}_{.col}")) %>% 
  ungroup() %>% 
  group_by(Sample_ID) %>%
  mutate(Remove = ifelse(all(c_across(starts_with("n_")) == 10), "FALSE", "TRUE")) %>% ## Check CV's, then remove 
  select(c(Sample_ID, Median_Extractable_NPOC_mg_per_L, Median_Extractable_NPOC_mg_per_kg, Median_Extractable_TN_mg_per_L, Median_Extractable_TN_mg_per_kg))%>% 
  ungroup()

# Gravimetric Moisture ----------------------------------------------------

grav_inc = read.csv("C:/GitHub/ECA_Multireactor_Incubations/Data/EC_Sediment_Gravimetric_Moisture.csv", skip = 2) %>% 
  slice(-1:-11) %>% 
  filter(Field_Name != "#End_Data") %>% 
  select(-c(Field_Name, IGSN, Material)) %>% 
  mutate(across(c(Initial_Water_Mass_g:Incubation_Water_Mass_g), as.numeric))

#Some dry reps have high CV: EC_057 (low moisture), EC_081 (one lower sample), EC_063 (low moisture), EC_088 (one lower sample), EC_076 (low moisture), EC_071 (one lower sample), EC_056 (low moisture), EC_069 (one slightly lower) 

# median_grav = grav_inc %>% 
#   mutate(X62948_Initial_Gravimetric_Moisture_g_per_g = ifelse(grepl("INC_Method_001|INC_Method_002", Methods_Deviation), NA, X62948_Initial_Gravimetric_Moisture_g_per_g)) %>% 
#   mutate(X62948_Final_Gravimetric_Moisture_g_per_g = ifelse(grepl("INC_Method_001|INC_Method_002", Methods_Deviation), NA, X62948_Final_Gravimetric_Moisture_g_per_g)) %>% 
#   #missing replicates (EC_072-W5/D5), less sediment in sample (EC_012-D5)
#   separate(Sample_Name, c("Sample_ID", "Rep"), sep = "-") %>% 
#   mutate(Rep = if_else(grepl("D", Rep), "D", "W")) %>%
#   group_by(Sample_ID, Rep) %>%
#   summarise(across(where(is.numeric),
#                    list(Median = ~median(.x, na.rm = TRUE),
#                         cv = ~sd(.x, na.rm = TRUE)/mean(.x, na.rm =TRUE),
#                         n = ~sum(!is.na(.x))), 
#                    .names = "{.fn}_{.col}")) %>% 
#   ungroup() %>% 
#   group_by(Sample_ID, Rep) %>%
#   mutate(Remove = ifelse(all(c_across(starts_with("n_")) == 5), "FALSE", "TRUE")) %>% ## Check CV's, then remove 
#   select(c(Sample_ID, Rep, Median_X62948_Initial_Gravimetric_Moisture_g_per_g, Median_X62948_Final_Gravimetric_Moisture_g_per_g)) 

## Fe #### 

fe = read.csv("C:/GitHub/ECA_Multireactor_Incubations/Data/EC_Sediment_Fe.csv", skip = 2) %>% 
  filter(grepl("EC", Sample_Name)) %>% 
  filter(!grepl("EC_011|EC_012|EC_023|EC_052|EC_053|EC_057", Sample_Name)) %>%  
  select(-c(Field_Name, IGSN, Material)) %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "SFE", "INC")) %>% 
  rename(Dev_Fe = Methods_Deviation)

median_iron = fe  %>% 
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
  group_by(Sample_ID) %>%
  summarise(across(where(is.numeric),
                   list(Median = ~median(.x, na.rm = TRUE),
                        cv = ~sd(.x, na.rm = TRUE)/mean(.x, na.rm =TRUE),
                        n = ~sum(!is.na(.x))), 
                   .names = "{.fn}_{.col}")) %>% 
  ungroup() %>% 
  group_by(Sample_ID) %>%
  mutate(Remove = ifelse(all(c_across(starts_with("n_")) == 5), "FALSE", "TRUE")) %>% ## Check CV's, then remove 
  select(c(Sample_ID, Median_Fe_mg_per_L, Median_Fe_mg_per_kg))%>% 
  ungroup()

## Median Data ####

## Wet Dry Medians

all_medians = median_chem %>% 
  left_join(median_atp, by = "Sample_ID") %>% 
  left_join(median_iron, by = "Sample_ID") %>% 
  left_join(median_cn, by = "Sample_ID") %>% 
  left_join(median_npoc_tn, by = "Sample_ID") %>% 
  mutate(Sample_Name = str_replace(Sample_ID, "INC", "all")) %>% 
  select(-c(Sample_ID)) %>% 
  relocate(Sample_Name, .before = Median_SpC_microsiemens_per_cm)

## Read in Median Data to get Dry moisture values

median = read.csv("C:/GitHub/ECA_Multireactor_Incubations/Data/EC_Sediment_Sample_Data_Summary.csv", skip = 2) %>% 
  filter(grepl("EC", Sample_Name)) %>% 
  separate(Sample_Name, c("Sample_Name", "Rep"), sep = "-") %>% 
  mutate(Sample_Name = paste0(Sample_Name, "_all")) %>% 
  filter(!grepl("EC_011|EC_012|EC_023|EC_052|EC_053|EC_057", Sample_Name)) %>%  
  select(-c(Field_Name, IGSN, Material, Median_Missing_Reps, Median_Respiration_Rate_mg_DO_per_kg_per_H)) %>% 
  select(-matches("per_L")) %>% 
  rename_with(~ str_remove_all(., "_[0-9]+")) %>% 
  rename_with(~ str_replace(., "^(([^_]*_){2}[^_]*).*", "\\1")) %>%
  rename_with(~ str_replace_all(., "Median", "median")) %>% 
  rename(median_SpC = median_SpC_microsiemens) %>% 
  rename(median_Temp = median_Temperature_degC)

median_dry = median %>% 
  filter(Rep == "D") %>% 
  select(c(Sample_Name, median_Initial_Gravimetric, median_Final_Gravimetric)) %>% 
  mutate(across(c(median_Initial_Gravimetric:median_Final_Gravimetric), as.numeric)) 

## Effect Size Data ####

effect_pub = read.csv("C:/GitHub/ECA_Multireactor_Incubations/Data/EC_Sediment_Effect_Size.csv", skip = 2) %>%
  filter(grepl("EC", Sample_Name)) %>%
  filter(!grepl("EC_011|EC_012|EC_023|EC_052|EC_053|EC_057", Sample_Name)) %>%
 select(-c(IGSN, Field_Name, Material, Methods_Deviation, Effect_Size_Respiration_Rate_mg_DO_per_L_per_H, Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H))

## Read in grain size/ssa variables ####

grain = read.csv("C:/GitHub/ECA_Multireactor_Incubations/Data/v3_CM_SSS_Sediment_Sample_Data_Summary.csv", skip = 2) %>% 
  filter(grepl("CM", Sample_Name)) %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "CM", "EC")) %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "Sediment", "all")) %>% 
  select(c(Sample_Name, Percent_Tot_Sand, Percent_Coarse_Sand, Percent_Med_Sand, Percent_Fine_Sand, Percent_Silt, Percent_Clay, Mean_Specific_Surface_Area_m2_per_g))

## Join all data

effect_data = left_join(effect_respiration, grain, by = "Sample_Name") %>% 
  left_join(effect_pub) %>% 
  left_join(all_medians) %>% 
  left_join(median_dry) %>% 
  mutate(across(c(Effect_Size_SpC_microsiemens_per_cm:Mean_Specific_Surface_Area_m2_per_g), as.numeric)) %>%  # make data numeric 
    select(-c(Effect_Size_Initial_Gravimetric_Moisture_g_per_g, Effect_Size_Final_Gravimetric_Moisture_g_per_g)) %>% 
  rename(median_Dry_Initial_Gravimetric = median_Initial_Gravimetric) %>% 
  rename(median_Dry_Final_Gravimetric = median_Final_Gravimetric) %>% 
  mutate(Effect_Size_Respiration_Rate_mg_DO_per_L_per_H = abs(Effect_Size_Respiration_Rate_mg_DO_per_L_per_H)) %>% 
  mutate(Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H = abs(Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H)) %>% 
  mutate(across(c(Effect_Size_Respiration_Rate_mg_DO_per_L_per_H:median_Dry_Final_Gravimetric), as.numeric))
  
# Transform Data ----------------------------------------------------------

## Cube PCA for LASSO####

# Fe outlier not in analysis
cube_effect = effect_data %>% 
  mutate(across(where(is.numeric), cube_root)) %>% # cube root transform data
  rename_with(where(is.numeric), .fn = ~ paste0("cube_", .x)) %>% 
  column_to_rownames("Sample_Name") %>%
  filter(cube_Effect_Size_Fe_mg_per_kg > -1) %>%  # remove Fe outlier for analysis 
  select(-contains("per_L")) %>% 
  relocate(cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H, .before = cube_Effect_Size_SpC_microsiemens_per_cm)

# Data frame that includes Fe outlier
fe_cube_effect = effect_data %>% 
  mutate(across(where(is.numeric), cube_root)) %>% # cube root transform data
  rename_with(where(is.numeric), .fn = ~ paste0("cube_", .x)) %>% 
  column_to_rownames("Sample_Name") %>%
  select(-contains("per_L")) %>% 
  filter(cube_Effect_Size_Fe_mg_per_kg < -1) %>%
  relocate(cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H, .before = cube_Effect_Size_SpC_microsiemens_per_cm)

# dry_cube_effect = cube_effect_data %>% 
#   filter(Rep == "D") %>% 
#   column_to_rownames("Sample_Name") %>% 
#   select(-c(Rep))
# 
# wet_cube_effect = cube_effect_data %>% 
#   filter(Rep == "W")%>% 
#   column_to_rownames("Sample_Name") %>% 
#   select(-c(Rep))

## Pearson Correlation Matrix ####

## All Medians
# scale data before it goes into correlation matrix

scale_cube_effect = as.data.frame(scale(cube_effect, center = T, scale = T))

scale_cube_effect_pearson <- cor(scale_cube_effect, method = "pearson")

if (print.images == T){
  
  png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Cube_Median_Effect_Pearson_Correlation_Matrix_Half_Rates.png"), width = 12, height = 12, units = "in", res = 300)
  
  corrplot(scale_cube_effect_pearson,type = "upper", method = "number", tl.col = "black", tl.cex = 0.8, cl.cex = 0.5, number.cex = 0.5,  title = "Effect Samples Pearson Correlation")
  
}

dev.off()

# make one line correlation matrix with just effect size

corr_effect = matrix(scale_cube_effect_pearson[1, ], nrow = 1)


colnames(corr_effect) = colnames(scale_cube_effect_pearson)

rownames(corr_effect) = rownames(scale_cube_effect_pearson)[1]

# Try to plot pearson and LASSO together as two lines

corr_effect_df = as.data.frame(corr_effect) %>% 
  reshape2::melt() %>% 
  rename(Coefficients = value) %>% 
  filter(Coefficients != 1) %>% 
  mutate(y = "Pearson")


## Dry Medians 
# scale data before it goes into correlation matrix

scale_cube_dry_effect = as.data.frame(scale(dry_cube_effect))

scale_cube_dry_effect_pearson <- cor(scale_cube_dry_effect, method = "pearson")

if (print.images == T){
  
  png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Cube_Dry_Median_Effect_Pearson_Correlation_Matrix.png"), width = 12, height = 12, units = "in", res = 300)
  
  corrplot(scale_cube_dry_effect_pearson,type = "upper", method = "number", tl.col = "black", tl.cex = 1, cl.cex = 1,  title = "Effect Samples Pearson Correlation")
  
}

dev.off()

## Wet Medians
# scale data before it goes into correlation matrix

scale_cube_wet_effect = as.data.frame(scale(wet_cube_effect))

scale_cube_wet_effect_pearson <- cor(scale_cube_wet_effect, method = "pearson")

if (print.images == T){
  
  png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Cube_Wet_Median_Effect_Pearson_Correlation_Matrix.png"), width = 12, height = 12, units = "in", res = 300)
  
  corrplot(scale_cube_wet_effect_pearson,type = "upper", method = "number", tl.col = "black", tl.cex = 1, cl.cex = 1,  title = "Effect Samples Pearson Correlation")
  
}

dev.off()


## Downselected LASSO - Loop through coefficients to choose for LASSO ####

## All Medians

all_pearson_df = as.data.frame(scale_cube_effect_pearson)

row_names_pearson <- rownames(all_pearson_df)

all_pearson_df$Variable <- row_names_pearson

# Melt the dataframe for plotting
all_pearson_melted <- reshape2::melt(all_pearson_df, id.vars = "Variable") %>% 
  filter(value != 1) %>% 
  mutate(value = abs(value)) %>% 
  filter(!grepl("Respiration", Variable))

all_effect_melted <- all_pearson_melted %>% 
  filter(grepl("Effect_Size_Respiration", variable)) %>%
  filter(!grepl("Silt", Variable)) # remove silt from variables, 0 values so not using

all_choose_melted <- all_pearson_melted %>% 
  filter(!grepl("Respiration", variable)) %>%
  filter(!grepl("Silt", variable)) %>%
  filter(!grepl("Silt", Variable)) %>% #try removing silt (0 values)
  #distinct(value, .keep_all = TRUE) %>% 
  left_join(all_effect_melted, by = "Variable") %>% 
  rename(Variable_1 = Variable) %>% 
  rename(Variable_2 = variable.x) %>% 
  rename(Correlation = value.x) %>% 
  rename(Variable_1_Effect_Correlation = value.y) %>% 
  select(-c(variable.y)) %>% 
  left_join(all_effect_melted, by = c("Variable_2" = "Variable")) %>% 
  rename(Variable_2_Effect_Correlation = value) %>% 
  select(-c(variable))

all_loop_melt = all_choose_melted %>% 
  arrange(desc(Correlation))

# Pearson correlation coefficient to remove above
correlation = 0.7

## Start loop to remove highly correlated (> 0.5)
all_effect_filter = function(all_loop_melt) {
  
  rows_to_keep = rep(TRUE, nrow(all_loop_melt))
  
  for (i in seq_len(nrow(all_loop_melt))) {
    
    if (!rows_to_keep[i]) next
    
    row = all_loop_melt[i, ]
    
    if (row$Correlation < correlation) next
    
    if(row$Variable_1_Effect_Correlation >= row$Variable_2_Effect_Correlation) {
      
      var_to_keep = row$Variable_1
      var_to_remove = row$Variable_2
      
    } else {
      
      var_to_keep = row$Variable_2
      var_to_remove = row$Variable_1
      
    }
    
    all_loop_melt$Variable_to_Keep[i] = var_to_keep
    all_loop_melt$Variable_to_Remove[i] = var_to_remove
    
    for (j in seq(i + 1, nrow(all_loop_melt))) {
      
      if(all_loop_melt$Variable_1[j] == var_to_remove || all_loop_melt$Variable_2[j] == var_to_remove) {
        
        rows_to_keep[j] = FALSE
        
      }
      
    }
    
    
  }
  
  return(all_loop_melt[rows_to_keep, ])
  
}

all_filtered_data = all_effect_filter(all_loop_melt) 

# pull out variables to remove
all_removed_variables = all_filtered_data %>% 
  distinct(Variable_to_Remove)

# pull out all variables 
all_variables = all_effect_melted %>% 
  select(c(Variable))

# remove variables from all variables to get variables to keep for LASSO 
all_kept_variables = all_effect_melted[!(all_effect_melted$Variable %in% all_removed_variables$Variable_to_Remove), ] #keeps SpC, Temp, pH, ATP, NPOC, TN (ext), TOC, TN (solid), med sand, silt

# if silt is removed, keeps SpC, Temp, pH, ATP, NPOC, TOC, TN (solid), med sand, fine sand, lost grav. moisture

## LASSO VARIABLES ####

# Keep variables selected from down-selected correlation matrix and add Cube_Effect_Size
col_to_keep = unique(all_kept_variables$Variable)
col_to_keep = c(col_to_keep, "cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H")

all_cube_variables = cube_effect[, col_to_keep, drop = FALSE]

## LASSO with Correlation Matrix Selected Variables ####
set.seed(42)
## Set response variable (Cube_Effect_Size) and scale
yvar <- data.matrix(scale(all_cube_variables$cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H, center = TRUE, scale = TRUE))
mean(yvar)
sd(yvar)

## Set predictor variables and scale
exclude_col = "cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H"

all_x_cube_variables = as.data.frame(scale(all_cube_variables[, !(names(all_cube_variables) %in% exclude_col)], center = T, scale = T))
#mean(x_cube_variables$Cube_SpC_Diff)
#sd(x_cube_variables$Cube_SpC_Diff)

xvars <- data.matrix(all_x_cube_variables)

lasso = cv.glmnet(xvars, yvar, alpha = 1, nfolds = 5,
                  ,standardize = FALSE, standardize.response = FALSE, intercept = FALSE
                  #,standardize = TRUE, standardize.response = TRUE, intercept = FALSE
                  # , standardize = TRUE, standardize.response = FALSE, intercept = FALSE
)

best_lambda <- lasso$lambda.min
best_lambda

plot(lasso)

best_lasso_model <- glmnet(xvars, yvar, alpha = 1, lambda = best_lambda, family = "gaussian",
                           , standardize = FALSE, standardize.response = FALSE, intercept = FALSE
                           #  , standardize = TRUE, standardize.response = TRUE, intercept = FALSE
                           #, standardize = TRUE, standardize.response = FALSE, intercept = FALSE
)

#coef(best_lasso_model)
all_lasso_coefs = coef(best_lasso_model)
#Fine Sand (0.26), ATP (0.19), %N (0.07), Effect Fe (0.07), Effect SpC (0.05)
yvar_predict <- predict(best_lasso_model, s = best_lambda, newx = xvars)

sst <- sum((yvar - mean(yvar))^2)
sse <- sum((yvar_predict - yvar)^2)

rsq = 1 - sse/sst

rsq #0.45

ds_lasso_df = as.data.frame(as.matrix(all_lasso_coefs))

colnames(ds_lasso_df) = c("Coefficients")

ds_lasso_df = ds_lasso_df %>% 
  rownames_to_column(var = "variable") %>% 
  slice(-1) 

ds_lasso_df$y = "LASSO"


## Dry Medians 

# 1) Pivot data frame and sort highest to lowest

dry_pearson_df <- as.data.frame(scale_cube_dry_effect_pearson)

row_names_dry_pearson <- rownames(dry_pearson_df)

dry_pearson_df$Variable <- row_names_dry_pearson

# Melt the dataframe for plotting
dry_pearson_melted <- reshape2::melt(dry_pearson_df, id.vars = "Variable") %>% 
  filter(value != 1) %>% 
  mutate(value = abs(value)) %>% 
  filter(!grepl("Respiration", Variable))

dry_effect_melted <- dry_pearson_melted %>% 
  filter(grepl("Effect_Size", variable)) %>%
  filter(!grepl("Silt", Variable)) # remove silt from variables, lots of 0 values so not using

dry_choose_melted <- dry_pearson_melted %>% 
  filter(!grepl("Respiration", variable)) %>%
  filter(!grepl("Silt", variable)) %>%
  filter(!grepl("Silt", Variable)) %>% #try removing silt (0 values)
  #distinct(value, .keep_all = TRUE) %>% 
  left_join(dry_effect_melted, by = "Variable") %>% 
  rename(Variable_1 = Variable) %>% 
  rename(Variable_2 = variable.x) %>% 
  rename(Correlation = value.x) %>% 
  rename(Variable_1_Effect_Correlation = value.y) %>% 
  select(-c(variable.y)) %>% 
  left_join(dry_effect_melted, by = c("Variable_2" = "Variable")) %>% 
  rename(Variable_2_Effect_Correlation = value) %>% 
  select(-c(variable))

dry_loop_melt = dry_choose_melted %>% 
  arrange(desc(Correlation))

# Pearson correlation coefficient to remove above
correlation = 0.7

## Start loop to remove highly correlated (> 0.5)
dry_effect_filter = function(dry_loop_melt) {
  
  rows_to_keep = rep(TRUE, nrow(dry_loop_melt))
  
  for (i in seq_len(nrow(dry_loop_melt))) {
    
    if (!rows_to_keep[i]) next
    
    row = dry_loop_melt[i, ]
    
    if (row$Correlation < correlation) next
    
    if(row$Variable_1_Effect_Correlation >= row$Variable_2_Effect_Correlation) {
      
      var_to_keep = row$Variable_1
      var_to_remove = row$Variable_2
      
    } else {
      
      var_to_keep = row$Variable_2
      var_to_remove = row$Variable_1
      
    }
    
    dry_loop_melt$Variable_to_Keep[i] = var_to_keep
    dry_loop_melt$Variable_to_Remove[i] = var_to_remove
    
    for (j in seq(i + 1, nrow(dry_loop_melt))) {
      
      if(dry_loop_melt$Variable_1[j] == var_to_remove || dry_loop_melt$Variable_2[j] == var_to_remove) {
        
        rows_to_keep[j] = FALSE
        
      }
      
    }
    
    
  }
  
  return(dry_loop_melt[rows_to_keep, ])
  
}

dry_filtered_data = dry_effect_filter(dry_loop_melt) 

# pull out variables to remove
dry_removed_variables = dry_filtered_data %>% 
  distinct(Variable_to_Remove)

# pull out all variables 
dry_all_variables = dry_effect_melted %>% 
  select(c(Variable))

# remove variables from all variables to get variables to keep for LASSO 
dry_kept_variables = dry_effect_melted[!(dry_effect_melted$Variable %in% dry_removed_variables$Variable_to_Remove), ] #keeps SpC, Temp, pH, ATP, NPOC, TN (ext), TOC, TN (solid), med sand, silt

# if silt is removed, keeps SpC, Temp, pH, ATP, NPOC, TOC, TN (solid), med sand, fine sand, lost grav. moisture

## LASSO VARIABLES ####

# Keep variables selected from down-selected correlation matrix and add Cube_Effect_Size
dry_col_to_keep = unique(dry_kept_variables$Variable)
dry_col_to_keep = c(dry_col_to_keep, "Effect_Size")

dry_cube_variables = dry_cube_effect[, dry_col_to_keep, drop = FALSE]

## LASSO with Correlation Matrix Selected Variables ####
set.seed(42)
## Set response variable (Cube_Effect_Size) and scale
yvar <- data.matrix(scale(dry_cube_variables$Effect_Size, center = TRUE, scale = TRUE))
mean(yvar)
sd(yvar)

## Set predictor variables and scale
exclude_col = "Effect_Size"

dry_x_cube_variables = as.data.frame(scale(dry_cube_variables[, !(names(dry_cube_variables) %in% exclude_col)], center = T, scale = T))
#mean(x_cube_variables$Cube_SpC_Diff)
#sd(x_cube_variables$Cube_SpC_Diff)

xvars <- data.matrix(dry_x_cube_variables)

lasso = cv.glmnet(xvars, yvar, alpha = 1, nfolds = 5,
                  ,standardize = FALSE, standardize.response = FALSE, intercept = FALSE
                  #,standardize = TRUE, standardize.response = TRUE, intercept = FALSE
                  # , standardize = TRUE, standardize.response = FALSE, intercept = FALSE
)

best_lambda <- lasso$lambda.min
best_lambda

plot(lasso)

best_lasso_model <- glmnet(xvars, yvar, alpha = 1, lambda = best_lambda, family = "gaussian",
                           , standardize = FALSE, standardize.response = FALSE, intercept = FALSE
                           #  , standardize = TRUE, standardize.response = TRUE, intercept = FALSE
                           #, standardize = TRUE, standardize.response = FALSE, intercept = FALSE
)

#lasso_coefs = coef(best_lasso_model)
ds_lasso_coefs = coef(best_lasso_model)
#Fine Sand, ATP, SpC, N% Tot Sand
yvar_predict <- predict(best_lasso_model, s = best_lambda, newx = xvars)

sst <- sum((yvar - mean(yvar))^2)
sse <- sum((yvar_predict - yvar)^2)

rsq = 1 - sse/sst

rsq #0.455

## Wet Medians 

# 1) Pivot data frame and sort highest to lowest

wet_pearson_df <- as.data.frame(scale_cube_wet_effect_pearson)

row_names_wet_pearson <- rownames(wet_pearson_df)

wet_pearson_df$Variable <- row_names_wet_pearson

# Melt the dataframe for plotting
wet_pearson_melted <- reshape2::melt(wet_pearson_df, id.vars = "Variable") %>% 
  filter(value != 1) %>% 
  mutate(value = abs(value)) %>% 
  filter(!grepl("Respiration", Variable))

wet_effect_melted <- wet_pearson_melted %>% 
  filter(grepl("Effect_Size", variable)) %>%
  filter(!grepl("Silt", Variable)) # remove silt from variables, lots of 0 values so not using

wet_choose_melted <- wet_pearson_melted %>% 
  filter(!grepl("Respiration", variable)) %>%
  filter(!grepl("Silt", variable)) %>%
  filter(!grepl("Silt", Variable)) %>% #try removing silt (0 values)
  #distinct(value, .keep_all = TRUE) %>% 
  left_join(wet_effect_melted, by = "Variable") %>% 
  rename(Variable_1 = Variable) %>% 
  rename(Variable_2 = variable.x) %>% 
  rename(Correlation = value.x) %>% 
  rename(Variable_1_Effect_Correlation = value.y) %>% 
  select(-c(variable.y)) %>% 
  left_join(wet_effect_melted, by = c("Variable_2" = "Variable")) %>% 
  rename(Variable_2_Effect_Correlation = value) %>% 
  select(-c(variable))

wet_loop_melt = wet_choose_melted %>% 
  arrange(desc(Correlation))

# Pearson correlation coefficient to remove above
correlation = 0.7

## Start loop to remove highly correlated (> 0.5)
wet_effect_filter = function(wet_loop_melt) {
  
  rows_to_keep = rep(TRUE, nrow(wet_loop_melt))
  
  for (i in seq_len(nrow(wet_loop_melt))) {
    
    if (!rows_to_keep[i]) next
    
    row = wet_loop_melt[i, ]
    
    if (row$Correlation < correlation) next
    
    if(row$Variable_1_Effect_Correlation >= row$Variable_2_Effect_Correlation) {
      
      var_to_keep = row$Variable_1
      var_to_remove = row$Variable_2
      
    } else {
      
      var_to_keep = row$Variable_2
      var_to_remove = row$Variable_1
      
    }
    
    wet_loop_melt$Variable_to_Keep[i] = var_to_keep
    wet_loop_melt$Variable_to_Remove[i] = var_to_remove
    
    for (j in seq(i + 1, nrow(wet_loop_melt))) {
      
      if(wet_loop_melt$Variable_1[j] == var_to_remove || wet_loop_melt$Variable_2[j] == var_to_remove) {
        
        rows_to_keep[j] = FALSE
        
      }
      
    }
    
    
  }
  
  return(wet_loop_melt[rows_to_keep, ])
  
}

wet_filtered_data = wet_effect_filter(wet_loop_melt) 

# pull out variables to remove
wet_removed_variables = wet_filtered_data %>% 
  distinct(Variable_to_Remove)

# pull out all variables 
wet_all_variables = wet_effect_melted %>% 
  select(c(Variable))

# remove variables from all variables to get variables to keep for LASSO 
wet_kept_variables = wet_effect_melted[!(wet_effect_melted$Variable %in% wet_removed_variables$Variable_to_Remove), ] #keeps SpC, Temp, pH, ATP, NPOC, TN (ext), TOC, TN (solid), med sand, silt

# if silt is removed, keeps SpC, Temp, pH, ATP, NPOC, TOC, TN (solid), med sand, fine sand, lost grav. moisture

## LASSO VARIABLES ####

# Keep variables selected from down-selected correlation matrix and add Cube_Effect_Size
wet_col_to_keep = unique(wet_kept_variables$Variable)
wet_col_to_keep = c(wet_col_to_keep, "Effect_Size")

wet_cube_variables = wet_cube_effect[, wet_col_to_keep, drop = FALSE]

## LASSO with Correlation Matrix Selected Variables ####
set.seed(42)
## Set response variable (Cube_Effect_Size) and scale
yvar <- data.matrix(scale(wet_cube_variables$Effect_Size, center = TRUE, scale = TRUE))
mean(yvar)
sd(yvar)

## Set predictor variables and scale
exclude_col = "Effect_Size"

wet_x_cube_variables = as.data.frame(scale(wet_cube_variables[, !(names(wet_cube_variables) %in% exclude_col)], center = T, scale = T))
#mean(x_cube_variables$Cube_SpC_Diff)
#sd(x_cube_variables$Cube_SpC_Diff)

xvars <- data.matrix(wet_x_cube_variables)

lasso = cv.glmnet(xvars, yvar, alpha = 1, nfolds = 5,
                  ,standardize = FALSE, standardize.response = FALSE, intercept = FALSE
                  #,standardize = TRUE, standardize.response = TRUE, intercept = FALSE
                  # , standardize = TRUE, standardize.response = FALSE, intercept = FALSE
)

best_lambda <- lasso$lambda.min
best_lambda

plot(lasso)

best_lasso_model <- glmnet(xvars, yvar, alpha = 1, lambda = best_lambda, family = "gaussian",
                           , standardize = FALSE, standardize.response = FALSE, intercept = FALSE
                           #  , standardize = TRUE, standardize.response = TRUE, intercept = FALSE
                           #, standardize = TRUE, standardize.response = FALSE, intercept = FALSE
)

#lasso_coefs = coef(best_lasso_model)
ds_lasso_coefs = coef(best_lasso_model)
#Fine Sand, ATP, SpC, N% Tot Sand
yvar_predict <- predict(best_lasso_model, s = best_lambda, newx = xvars)

sst <- sum((yvar - mean(yvar))^2)
sse <- sum((yvar_predict - yvar)^2)

rsq = 1 - sse/sst

rsq #0.42

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Cube_Median_Effect_vs_ATP_Scatter.png"), width = 6, height = 6, units = "in", res = 300)

ggplot(cube_effect, aes(y = cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H, x = cube_Median_ATP_picomoles_per_g)) +
  geom_point() +
  theme_bw() +
  #stat_cor(data = fe_cube_out, size = 5, digits = 2, aes(label = paste(..rr.label.., ..p.label.., sep = "~`;`~")))+
  stat_cor(data = cube_effect, label.x = 1.1, label.y = 11, size = 4, digits = 2, aes(label = paste(..rr.label..)))+
  stat_cor(data = cube_effect, label.x = 1.1, label.y = 10.25, size = 4, digits = 2, aes(label = paste(..p.label..)))+
  stat_poly_line(data = cube_effect, se = FALSE) + 
  ylab("Effect Size Respiration Rate (mg/kg)") +
  xlab("Median ATP (pmol/g)")

dev.off()

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Cube_Median_Effect_vs_TOC_Scatter.png"), width = 6, height = 6, units = "in", res = 300)

ggplot(cube_effect, aes(y = cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H, x = cube_Median_X01395_C_percent_per_mg)) +
  geom_point() +
  theme_bw() +
  #stat_cor(data = fe_cube_out, size = 5, digits = 2, aes(label = paste(..rr.label.., ..p.label.., sep = "~`;`~")))+
  stat_cor(data = cube_effect, label.x = 0.55, label.y = 11, size = 4, digits = 2, aes(label = paste(..rr.label..)))+
  stat_cor(data = cube_effect, label.x = 0.55, label.y = 10.25, size = 4, digits = 2, aes(label = paste(..p.label..)))+
  stat_poly_line(data = cube_effect, se = FALSE) + 
  ylab("Effect Size Respiration Rate (mg/kg)") +
  xlab("Median TOC (%)")

dev.off()

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Cube_Median_Effect_vs_TN_Scatter.png"), width = 6, height = 6, units = "in", res = 300)

ggplot(cube_effect, aes(y = cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H, x = cube_Median_X01397_N_percent_per_mg)) +
  geom_point() +
  theme_bw() +
  #stat_cor(data = fe_cube_out, size = 5, digits = 2, aes(label = paste(..rr.label.., ..p.label.., sep = "~`;`~")))+
  stat_cor(data = cube_effect, label.x = 0.225, label.y = 11, size = 4, digits = 2, aes(label = paste(..rr.label..)))+
  stat_cor(data = cube_effect, label.x = 0.225, label.y = 10.25, size = 4, digits = 2, aes(label = paste(..p.label..)))+
  stat_poly_line(data = cube_effect, se = FALSE)+ 
  ylab("Effect Size Respiration Rate (mg/kg)") +
  xlab("Median TN (%)")

dev.off()

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Cube_Median_Effect_vs_Median_SpC_Scatter.png"), width = 6, height = 6, units = "in", res = 300)

ggplot(cube_effect, aes(y = cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H, x = cube_Median_SpC_microsiemens_per_cm)) +
  geom_point() +
  theme_bw() +
  #stat_cor(data = fe_cube_out, size = 5, digits = 2, aes(label = paste(..rr.label.., ..p.label.., sep = "~`;`~")))+
  stat_cor(data = cube_effect, label.x = 2.5, label.y = 11, size = 4, digits = 2, aes(label = paste(..rr.label..)))+
  stat_cor(data = cube_effect, label.x = 2.5, label.y = 10.25, size = 4, digits = 2, aes(label = paste(..p.label..)))+
  stat_poly_line(data = cube_effect, se = FALSE)+ 
  ylab("Effect Size Respiration Rate (mg/kg)") +
  xlab("Median SpC (\u03BCS/cm)")

dev.off()

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Cube_Median_Effect_vs_Effect_Fe_Scatter.png"), width = 6, height = 6, units = "in", res = 300)

ggplot(cube_effect, aes(y = cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H, x = cube_Effect_Size_Fe_mg_per_kg)) +
  geom_point() +
  geom_point(fe_cube_effect, mapping = aes(y = cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H, x = cube_Effect_Size_Fe_mg_per_kg), color = "red") +
  theme_bw() +
  #stat_cor(data = fe_cube_out, size = 5, digits = 2, aes(label = paste(..rr.label.., ..p.label.., sep = "~`;`~")))+
  stat_cor(data = cube_effect, label.x = -2.5, label.y = 11, size = 4, digits = 2, aes(label = paste(..rr.label..)))+
  stat_cor(data = cube_effect, label.x = -2.5, label.y = 10.25, size = 4, digits = 2, aes(label = paste(..p.label..)))+
  stat_poly_line(data = cube_effect, se = FALSE)+ 
  ylab("Effect Size Respiration Rate (mg/kg)") +
  xlab("Effect Size Fe (II) (mg/kg)")

dev.off()

## Combine Heat Map

lasso_pear_df = bind_rows(ds_lasso_df, corr_effect_df) %>% 
  rename(type = y) %>% 
  filter(!grepl("Silt", variable)) %>% 
  mutate(variable = ifelse(variable == "cube_Effect_Size_SpC_microsiemens_per_cm", " Effect Size Specific Conductivity (\u03BCS/cm)", ifelse(variable == "cube_Effect_Size_pH", "Effect Size pH", ifelse(variable == "cube_Effect_Size_Temperature_degC", " Effect Size Temperature (\u00B0C)",  ifelse(variable == "cube_Effect_Size_Fe_mg_per_kg", "Effect Size Fe (II) (mg/kg)", ifelse(variable == "cube_Effect_Size_ATP_picomoles_per_g", "Effect Size ATP (pmol/g)", ifelse(variable == "cube_Effect_Size_Extractable_NPOC_mg_per_kg", "Effect Size Extractable NPOC (mg/kg)", ifelse(variable == "cube_Effect_Size_Extractable_TN_mg_per_kg", "Effect Size Extractable TN (mg/kg)", ifelse(variable == "cube_Effect_Size_C_percent_per_mg", "Effect Size TOC (%)", ifelse(variable == "cube_Effect_Size_N_percent_per_mg", "Effect Size TN (%)", ifelse(variable == "cube_Percent_Tot_Sand", "Total Sand (%)", ifelse(variable == "cube_Percent_Med_Sand", "Medium Sand (%)", ifelse(variable == "cube_Percent_Fine_Sand", "Fine Sand (%)", ifelse(variable == "cube_Median_SpC_microsiemens_per_cm", "Median Specific Conductivity (\u03BCS/cm)", ifelse(variable == "cube_Median_pH", "Median pH", ifelse(variable == "cube_Median_Temperature_degC", "Median Temperature (\u00B0C)", ifelse(variable == "cube_Median_ATP_picomoles_per_g", "Median ATP (pmol/g)",  ifelse(variable == "cube_Median_Fe_mg_per_kg", "Median Fe (II) (mg/kg)", ifelse(variable == "cube_Median_X01395_C_percent_per_mg", "Median TOC (%)",   ifelse(variable ==  "cube_Median_X01397_N_percent_per_mg", "Median TN (%)", ifelse(variable == "cube_Median_Extractable_NPOC_mg_per_kg", "Median Extractable NPOC (mg/kg)", ifelse(variable == "cube_Median_Extractable_TN_mg_per_kg", "Median Extractable TN (mg/kg)", ifelse(variable == "cube_median_Dry_Initial_Gravimetric", "Median Initial Dry Gravimetric Moisture (g/g)", ifelse(variable == "cube_Percent_Coarse_Sand", "Coarse Sand (%)", ifelse(variable == "cube_Percent_Clay", "Clay (%)", ifelse(variable == "cube_Mean_Specific_Surface_Area_m2_per_g", "Specific Surface Area (m\u00B2/g)", ifelse(variable == "cube_median_Dry_Final_Gravimetric", "Median Final Dry Gravimetric Moisture (g/g)",
 variable)))))))))))))))))))))))))))

color_palette <- colorRampPalette(c("#B2182B", "#F7F7F7", "#2166AC"))(200)

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Combined_Heat_Matrix.png"), width = 12, height = 4, units = "in", res = 300)

ggplot(lasso_pear_df, aes(variable, type)) +
  geom_tile(fill = "white", color = "black") +
  geom_text(aes(label = round(Coefficients, 2), color = Coefficients), size = 3, fontface = "bold") + 
  scale_color_gradientn(colors = color_palette, 
                        limit = c(-1, 1),
                        guide = "none") +
  theme_bw() + 
  theme(aspect.ratio = 0.1, 
        axis.text.x = element_text(angle = 90, hjust = 0, size = 10), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())+#, 
  #axis.text.y = element_blank()) +
  scale_x_discrete(position = "top") #+
# labs(y = "LASSO Coefficients")

dev.off()

## Merge Matrices + Scatter Plots

combine_hm_image = image_read("C:/GitHub/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/2024-09-17_Combined_Heat_Matrix.png")

combine_label_image = image_annotate(combine_hm_image, "A", size = 100, location = "+25+50", color = "black")

fine_sand_image = image_read("C:/GitHub/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/2024-08-30_Cube_Median_Effect_vs_Fine_Sand_Scatter.png")

fine_sand_label_image = image_annotate(fine_sand_image, "B", size = 65, location = "+30+20", color = "black")

fine_sand_scale_image = image_scale(fine_sand_label_image, "65%")

atp_image = image_read("C:/GitHub/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/2024-09-17_Cube_Median_Effect_vs_ATP_Scatter.png")

atp_label_image = image_annotate(atp_image, "C", size = 65, location = "+30+20", color = "black")

atp_scale_image = image_scale(atp_label_image, "65%")

fe_image = image_read("C:/GitHub/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/2024-09-17_Cube_Median_Effect_vs_Effect_Fe_Scatter.png")

fe_label_image = image_annotate(fe_image, "E", size = 65, location = "+30+20", color = "black")

fe_scale_image = image_scale(fe_label_image, "65%")

tn_image = image_read("C:/GitHub/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/2024-09-17_Cube_Median_Effect_vs_TN_Scatter.png")

tn_label_image = image_annotate(tn_image, "D", size = 65, location = "+30+20", color = "black")

tn_scale_image = image_scale(tn_label_image, "65%")

spc_image = image_read("C:/GitHub/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/2024-09-17_Cube_Median_Effect_vs_Median_SpC_Scatter.png")

spc_label_image = image_annotate(spc_image, "F", size = 65, location = "+30+20", color = "black")

spc_scale_image = image_scale(spc_label_image, "65%")

toc_image = image_read("C:/GitHub/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/2024-09-17_Cube_Median_Effect_vs_TOC_Scatter.png")

toc_label_image = image_annotate(toc_image, "G", size = 65, location = "+30+20", color = "black")

toc_scale_image = image_scale(toc_label_image, "65%")

first_row_image = image_append(c(fine_sand_scale_image, atp_scale_image, tn_scale_image))

second_row_image = image_append(c(fe_scale_image, spc_scale_image, toc_scale_image))

whole_image = image_append(c(combine_label_image, first_row_image, second_row_image), stack = TRUE)

image_write(whole_image, path = "C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/2024-09-17_Combined_Heat_Map.png")



## Control Point Influence

cpi = all_data %>% 
  mutate(Respiration = abs(as.numeric(Respiration_Rate_mg_DO_per_kg_per_H))) %>% 
  filter(Respiration < 9000) %>% 
  select(c(Sample_Name, Respiration))

sum_cpi = sum(cpi$Respiration)

median_cpi = median(cpi$Respiration)

high_cpi = cpi %>% 
  filter(Respiration > median_cpi)

sum_median_cpi = sum(high_cpi$Respiration)

tot_cpi = sum_median_cpi/sum_cpi
