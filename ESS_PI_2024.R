library(tidyverse)
library(glmnet)
library(factoextra)
library(janitor)
library(corrplot)
library(patchwork)

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

#write.csv(means,"C:/GitHub/ECA_Multireactor_Incubations/Data/Cleaned Data/ECA_Means_ESS_PI.csv", row.names = FALSE )


cm_atp = read.csv("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ICON_ModEx_SSS/11_ATP/03_ProcessedData/CM_ATP_ReadyForBoye_04-05-2024.csv") %>% 
  separate(Sample_Name, c("Sample_ID", "Rep"), sep = "_ATP-") %>% 
  group_by(Sample_ID) %>% 
  summarise(mean_atp_pico = mean(ATP_picomol_per_g), mean_atp_nm = mean(ATP_nanomol_per_L)) %>% 
  mutate(Sample_ID = str_replace(Sample_ID, "CM", "EC"))


atp_comp = atp_means %>% 
  separate(Sample_Name, c("Sample_ID", "Treat"), sep = "_INC-") %>% 
    left_join(cm_atp) %>% 
  filter(mean_atp_pico != -9999)

ggplot(atp_comp, aes(x = mean_atp_pico, y = Mean_ATP_picomol_per_g)) +
  geom_point(aes(color = Treat)) + 
  xlab("ICON ATP") + 
  ylab("ECA ATP") +
  geom_abline(intercept = 0, slope = 1)

ggplot(means, aes(x = Mean_Fe_mg_per_kg, y = Mean_ATP_picomol_per_g)) + 
  geom_point(aes(color = Treat)) 

ggplot(means, aes(x = ))

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


eff_all$rank_atp = rank(eff_all$`ATP Diff.`)

eff_all$rank_effect = rank(eff_all$`Effect Size`)

cor_matrix <- cor.test(x = eff_all$`ATP Diff.`, y = eff_all$`Effect Size`, method = "spearman")

print(cor_matrix)

ggplot(eff_all, aes(x = rank_atp, y = rank_effect)) +
  geom_point() + 
  annotate("text", 
           x = max(eff_all$rank_atp), 
           y = max(eff_all$rank_effect), 
      label = expression(paste(rho, " = -0.46, p = 0.008")), hjust = 1, vjust = 1)


cor_matrix <- cor(eff_all, method = "spearman")

corrplot(cor_matrix, method = "circle", type = "upper", tl.cex = 0.6, number.cex = 0.4, diag = FALSE, tl.col = "black")

wet = means[means$Treat == "W", ]

dry = means[means$Treat == "D", ]

means$rank_atp = rank(means$Mean_ATP_picomol_per_g)

means$rank_resp = rank(means$Mean_OutliersRemoved_Respiration_Rate_mg_DO_per_kg_per_H)


cor_atp <- cor(x = means$rank_atp, y = means$rank_resp, method = "spearman")

dry = means[means$Treat == "D", ]

cor_dry = cor(x = dry$rank_atp, y = dry$rank_resp, method = "spearman") 

ggplot(eff_all, aes(x = `% Fine Sand`, y = `Effect Size`)) + 
  geom_point()

ggplot(means, aes(x = Mean_ATP_picomol_per_g, fill = Treat)) +
  geom_histogram()

## Cube Effect Size ####
cube_root <- function(x) sign(x) * (abs(x))^(1/3)

cube_effect = eff_all %>% 
  mutate(across(where(is.numeric), cube_root))

cube_matrix <- cor(cube_effect, method = "pearson")

corrplot(cube_matrix, method = "circle", type = "upper", tl.cex = 0.6, number.cex = 0.4, diag = FALSE, tl.col = "black")

fs_p = ggplot(cube_effect, aes(x = `% Fine Sand`, y = `Effect Size`)) + 
  geom_point() +
  theme_minimal() +
  xlab("Cube Root %Fine Sand") +
  ylab("Cube Root Effect Size")

atp_p = ggplot(cube_effect, aes(x = `ATP Diff.`, y = `Effect Size`)) +
  geom_point() +
  theme_minimal() +
  xlab("Cube Root ATP (pmol/g) Difference") +
  ylab("Cube Root Effect Size")

fe_p = ggplot(cube_effect, aes(x = `Fe Diff.`, y = `Effect Size`)) +
  geom_point() +
  theme_minimal() +
  xlab("Cube Root Fe (mg/kg) Difference") +
  ylab("Cube Root Effect Size")

merged_p = fe_p + atp_p + fs_p

layout <- plot_layout(
  ncol = 3,  # Number of columns in the layout
  widths = c(2,2)  # Width ratio for each column
)

# Display the merged plot with custom layout
merged_p + layout

time_zero_do = read.csv("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/INC/03_ProcessedData/ECA_Sediment_Incubations_Respiration_Rates_ReadyForBoye_2024-04-05.csv") %>% 
  select(c(Sample_Name, DO_Concentration_At_Incubation_Time_Zero)) %>% 
  filter(DO_Concentration_At_Incubation_Time_Zero != -9999)

atp_all = read.csv("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/ATP/03_ProcessedData/EC_ATP_ReadyForBoye_01-26-2024.csv") %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "ATP", "INC")) %>% 
  filter(ATP_picomol_per_g != -9999)

do_atp = left_join(time_zero_do, atp_all, by = "Sample_Name") %>% 
  mutate(Treat = case_when(grepl("W", Sample_Name) ~ "W", 
                           grepl("D", Sample_Name) ~ "D"))


ggplot(do_atp, aes(x = ATP_picomol_per_g, y = DO_Concentration_At_Incubation_Time_Zero)) +
  geom_point(aes(color = Treat))
