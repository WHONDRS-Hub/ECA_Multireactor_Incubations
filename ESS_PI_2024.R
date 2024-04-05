library(tidyverse)
library(glmnet)
library(factoextra)
library(janitor)
library(corrplot)

rm(list=ls());graphics.off()

atp_means = read.csv("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/ATP/03_ProcessedData/EC_ATP_Summary_ReadyForBoye_03-05-2024.csv") %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "ATP", "INC"))

fe_means = read.csv("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Fe/03_ProcessedData/EC_SFE_Summary_ReadyForBoye_03-05-2024.csv") %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "SFE", "INC"))

resp_means = read.csv("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/INC/03_ProcessedData/ECA_Sediment_Incubations_Respiration_Rates_Summary_ReadyForBoye_2024-03-05.csv")

moi_means = read.csv("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/INC/03_ProcessedData/ECA_Drying_Masses_Summary_ReadyForBoye_01-29-2024.csv") %>% 
  separate(Sample_Name, c("Sample_Name", "Rep"), sep = -1) %>% 
  mutate(Initial_Gravimetric_Water = Initial_Water_Mass_g/Dry_Sediment_Mass_g) %>% 
  mutate(Final_Gravimetric_Water = Final_Water_Mass_g/Dry_Sediment_Mass_g) %>% 
  mutate(Lost_Gravimetric_Water = Initial_Gravimetric_Water - Final_Gravimetric_Water) %>% 
  filter(!grepl("023", Sample_Name)) %>% 
  group_by(Sample_Name) %>% 
  summarise(across(c(Initial_Gravimetric_Water:Lost_Gravimetric_Water), mean))

grain <- read.csv("C:/Github/ECA_Multireactor_Incubations/Data/v3_CM_SSS_Sediment_Grain_Size.csv", skip = 2, header = TRUE) %>% 
  filter(!row_number() %in% c(1:11)) %>% 
  dplyr::select(-c(Field_Name, Material)) %>% 
  filter(!grepl("SSS", Sample_Name)) %>% 
  filter(Sample_Name != "") %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "CM", "EC")) %>% 
  mutate(Sample_ID = str_remove(Sample_Name, "_GRN")) %>% 
  mutate_at(c("Percent_Fine_Sand", "Percent_Med_Sand", "Percent_Coarse_Sand", "Percent_Tot_Sand", "Percent_Silt", "Percent_Clay"), as.numeric) %>% 
  select(-c(IGSN, Methods_Deviation, Sample_Name))

ssa <- read.csv("C:/GitHub/ECA_Multireactor_Incubations/Data/CM_SSS_Sediment_Specific_Surface_Area.csv", skip = 2, header = TRUE) %>% 
  filter(!row_number() %in% c(1:11)) %>% 
  dplyr::select(-c(Field_Name, Material, IGSN)) %>% 
  filter(!grepl("SSS", Sample_Name)) %>% 
  filter(Sample_Name != "") %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "CM", "EC")) %>% 
  separate(Sample_Name, c("Sample_ID", "Rep"), sep = -6) %>% 
  filter(Specific_Surface_Area_m2_per_g != -9999) %>% 
  filter(!grepl("Negative", Specific_Surface_Area_m2_per_g)) %>%
  mutate(Specific_Surface_Area_m2_per_g = as.numeric(Specific_Surface_Area_m2_per_g)) %>%
  group_by(Sample_ID) %>%
  summarise(mean_ssa = mean(Specific_Surface_Area_m2_per_g, na.rm = TRUE))

icr_thermo =  read.csv("C:/GitHub/ECA_Multireactor_Incubations/Data/ratios_thermodynamics.csv") %>%
  rename(Sample_ID = Location)

npoc_dir = "C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/NPOC_TN/03_ProcessedData"

npoc_rec <- list.files(path = npoc_dir, recursive = TRUE, pattern = "Compiled", full.names = TRUE)

remove_rows <- function(df) {
  df %>% filter(!grepl("Blk", Sample_ID))  # Remove rows where ID is 2
}

change_column_type <- function(df) {
  df %>% mutate(NPOC_mg_C_per_L = as.double(NPOC_mg_C_per_L))  # Change Value column to double
}

npoc_list = map(npoc_rec, read.csv) %>% 
  map(remove_rows) %>% 
  map(change_column_type)

npoc = npoc_list %>% 
  bind_rows() %>% 
  separate(Sample_ID, c("Sample_ID", "Treat"), sep = "-") %>% 
  separate(Treat, c("Treat", "Rep"), sep = -1) %>% 
  mutate(TN_mg_N_per_L = sub(".*\\|([0-9.]+)_ppm_Final_Corrected", "\\1", TN_mg_N_per_L)) %>% 
  mutate(TN_mg_N_per_L = as.numeric(TN_mg_N_per_L))

npoc_means = npoc %>% 
  group_by(Sample_ID, Treat) %>% 
  summarise(mean_npoc = mean(NPOC_mg_C_per_L), mean_tn = mean(TN_mg_N_per_L)) %>% 
  unite(Sample_ID, c("Sample_ID", "Treat"), sep = "-") %>% 
  mutate(Sample_Name = str_replace(Sample_ID, "SIR", "INC")) %>% 
  select(-c(Sample_ID))


means = left_join(resp_means, moi_means, by = "Sample_Name") %>% 
  left_join(atp_means) %>% 
  left_join(fe_means) %>% 
  left_join(npoc_means) %>% 
  select(-contains("SD")) %>% 
  separate(Sample_Name, c("Sample_ID", "Rep"), sep = -6, remove = FALSE) %>% 
  separate(Rep, c("INC", "Treat"), sep = -1) %>% 
  left_join(grain, by = "Sample_ID") %>% 
  left_join(ssa, by = "Sample_ID") %>% 
  filter(!grepl("052", Sample_ID)) %>% 
  filter(!grepl("053", Sample_ID)) %>% 
  filter(!grepl("057", Sample_ID)) %>% 
  filter(!grepl("023", Sample_ID)) %>% 
  select(-c(Material, INC, Sample_Name)) %>% 
  mutate(Mean_OutliersRemoved_Respiration_Rate_mg_DO_per_L_per_H = abs(Mean_OutliersRemoved_Respiration_Rate_mg_DO_per_L_per_H)) %>% 
  mutate(Mean_OutliersRemoved_Respiration_Rate_mg_DO_per_kg_per_H = abs(Mean_OutliersRemoved_Respiration_Rate_mg_DO_per_kg_per_H)) %>% 
  mutate(Mean_WithOutliers_Respiration_Rate_mg_DO_per_L_per_H = abs(Mean_WithOutliers_Respiration_Rate_mg_DO_per_L_per_H))%>% 
  mutate(Mean_WithOutliers_Respiration_Rate_mg_DO_per_kg_per_H = abs(Mean_WithOutliers_Respiration_Rate_mg_DO_per_kg_per_H))

ggplot(means, aes(x = Mean_Fe_mg_per_kg, y = Mean_ATP_picomol_per_g)) + 
  geom_point(aes(color = Treat))

## EFFECT SIZE CLEANING ####

## Calculate Effect Size 
best_effect = means %>% 
  filter(!grepl("011", Sample_ID)) %>% 
  filter(!grepl("012", Sample_ID)) %>% 
  group_by(Sample_ID) %>% 
  mutate(across(c(Mean_SpC:mean_tn), ~. [Treat == "W"] - .[Treat  == "D"])) %>% 
  rename_with(.cols = c(Mean_SpC:mean_tn), .fn = ~ paste0("diff_", .x)) %>% 
  distinct(Sample_ID, .keep_all = TRUE) %>% left_join(icr_thermo)  

eff_all = best_effect %>% 
  drop_na() %>% 
  select(-c(Treat, diff_Mean_WithOutliers_Respiration_Rate_mg_DO_per_L_per_H, diff_Mean_WithOutliers_Respiration_Rate_mg_DO_per_kg_per_H, diff_Mean_OutliersRemoved_Respiration_Rate_mg_DO_per_L_per_H, diff_Initial_Gravimetric_Water, diff_Lost_Gravimetric_Water, diff_Mean_ATP_nanomol_per_L, diff_Mean_Fe_mg_per_L, diff_Mean_SpC, diff_Mean_pH, diff_Mean_Temp, Percent_Clay, Percent_Med_Sand, Percent_Coarse_Sand, Percent_Tot_Sand, Percent_Silt)) %>% 
  column_to_rownames("Sample_ID") %>% 
  #rename(`SpC Diff.` = diff_Mean_SpC) %>% 
  #rename(`pH Diff.` = diff_Mean_pH) %>% 
  #rename(`Temp. Diff.` = diff_Mean_Temp) %>% 
  rename(`Effect Size` = diff_Mean_OutliersRemoved_Respiration_Rate_mg_DO_per_kg_per_H) %>% 
  rename(`Final Grav. Water Diff.` = diff_Final_Gravimetric_Water) %>% 
  rename(`ATP Diff.` = diff_Mean_ATP_picomol_per_g) %>% 
  rename(`Fe Diff.` = diff_Mean_Fe_mg_per_kg) %>% 
  rename(`NPOC Diff.` = diff_mean_npoc) %>% 
  rename(`TN Diff.` = diff_mean_tn) %>% 
  rename(`% Fine Sand` = Percent_Fine_Sand) %>% 
  rename(`Mean SSA` = mean_ssa) %>% 
  relocate(`Effect Size`, .before = `Final Grav. Water Diff.`)

cor_matrix <- cor(eff_all, method = "spearman")

corrplot(cor_matrix, method = "circle", type = "upper", tl.cex = 0.6, number.cex = 0.4, diag = FALSE, tl.col = "black")


ggplot(eff_all, aes(x = `ATP Diff.`, y = `Effect Size`)) + 
  geom_point()

ggplot(means, aes(x = Mean_ATP_picomol_per_g, fill = Treat)) +
  geom_histogram()
