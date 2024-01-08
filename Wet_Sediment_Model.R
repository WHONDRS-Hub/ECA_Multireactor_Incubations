library(ggplot2)
library(tidyverse)
library(dplyr)
library(readxl)

rm(list=ls());graphics.off()

pnnl.user = 'laan208'

setwd("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files")

## GRAIN SIZE ####

grain <- read.csv("C:/Github/ECA_Multireactor_Incubations/Data/v2_CM_SSS_Sediment_Grain_Size.csv", skip = 2, header = TRUE)

grain <- grain %>% 
  filter(!row_number() %in% c(1:11)) %>% 
  dplyr::select(-c(Field_Name, Material)) 

grain_new <- read.csv(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ICON_ModEx_SSS/09_Grain_Size/03_ProcessedData/20230721_Grain_Size_SBR_RC4_CM_R21/20230721_Data_Processed_Grain_Size_SBR_RC4_CM_R21.csv"))

grain <- grain %>% 
  mutate(Sample_ID = Sample_Name) %>% 
  dplyr::select(-c(Sample_Name))

grain_all <- rbind(grain, grain_new)

grain_all$Sample_ID <- str_replace(grain_all$Sample_ID, "CM", "EC") 

grain_all <- grain_all %>% 
  filter(!grepl("SSS", Sample_ID)) %>% 
  separate(Sample_ID, into = c("EC", "Kit", "GRN"), sep = "_") %>% 
  unite("Sample_ID", EC:Kit, sep = "_") %>% 
  dplyr::select(-c(GRN)) %>% 
  filter(!grepl("NA", Sample_ID))

## KNOWN INCUBATION WEIGHTS ####

inc <- read_xlsx(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/INC/01_RawData/2022_Data_Raw_INC_Tube_Weight_ECA_EC.xlsx"))

##fix mistyped samples??

inc <- inc %>% 
  mutate(`INC_tube_50mL_wet_Sediment_with_Water_g (After Subsampling for SFE + SpC)`= 
           ifelse(`INC_tube_50mL_wet_Sediment_with_Water_g (After Subsampling for SFE + SpC)` >= 80, `INC_tube_50mL_wet_Sediment_with_Water_g (After Subsampling for SFE + SpC)` - 10,
                  `INC_tube_50mL_wet_Sediment_with_Water_g (After Subsampling for SFE + SpC)`)) %>% 
  mutate(mass_sed_water = (`INC_tube_50mL_wet_Sediment_with_Water_g (After Subsampling for SFE + SpC)` - INC_tube_50ml_empty_g) + 6) %>% 
  mutate(Sample_Name = SampleID) %>% 
  dplyr::select(-c(SampleID)) %>% 
  filter(!grepl("EV", Sample_Name)) %>% 
  #na.omit(mass_sed_water) %>% 
  relocate(Sample_Name, .before = INC_tube_50ml_empty_g)

## Compare to SSA
ssa <- read_csv(paste0("C:/GitHub/ECA_Multireactor_Incubations/Data/eca_ssa_predatapackage.csv"))

ssa <- ssa %>% 
  separate(Parent_ID, c("EC", "Site"), sep = "_")

for (i in 1:nrow(ssa)){
  
  if (str_count(ssa$Site[i], "[0-9]") <= 2){
    
    ssa$Site[i] = paste0("0", ssa$Site[i])
    
  }
  
  else {
    
    ssa$Site[i] = ssa$Site[i]
  }
  
}

ssa <- ssa %>% 
  unite(Sample_ID, c("EC", "Site"), sep = "_")

mean_ssa <- ssa %>% 
  group_by(Sample_ID) %>% 
  mutate(average_ssa = mean(ssa_m2_g)) %>% 
  distinct(Sample_ID, .keep_all = TRUE) %>% 
  dplyr::select(Sample_ID, average_ssa)

## ORIGINAL INCUBATION WEIGHTS ####

og_inc <- read.csv(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/INC/03_ProcessedData/ECA_Drying_Masses_Summary_merged_by_laan208_on_2023-12-01.csv"))

## MOISTURE TIN DATA ####

moisture <- read_csv(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/MOI/03_ProcessedData/ECA_Moisture_Outliers_Removed.csv"))

moisture_clean <- moisture %>%
  filter(is.na(flag)) %>% 
  separate(Sample_Name, c("EC", "Site", "Rep"), sep = "_", remove = FALSE) %>% 
  unite(Sample_ID, c("EC", "Site")) %>% 
  group_by(Sample_ID) %>% 
  mutate(average_grav = mean(percent_water_content_dry)) %>% 
  mutate(cv = (sd(percent_water_content_dry)/mean(percent_water_content_dry))*100)

moisture_average <- moisture_clean %>% 
  distinct(Sample_ID, .keep_all = TRUE) %>% 
  select(c(Sample_ID, average_grav))

## Regression of Moisture by Grain Size

moi_grn <- left_join(moisture_average, grain_all, by = c("Sample_ID"))

moi_grn$Percent_Clay <- as.numeric(moi_grn$Percent_Clay)

moi_grn$Percent_Silt <- as.numeric(moi_grn$Percent_Silt)

## MERGE DATA FOR MASS WATER####

inc_data <- left_join(og_inc, inc, by = c("Sample_Name")) %>% 
  mutate(mass_water = mass_sed_water - Dry_Sediment_Mass_g) %>% separate(Sample_Name, into = c("EC", "Kit", "REP"), sep = "_", remove = FALSE) %>% 
  unite("Sample_ID", EC:Kit, sep = "_") %>% 
  drop_na(mass_water) 

ggplot(inc_data, aes (x = mass_water))+
  geom_histogram()

mean_water <- mean(inc_data$mass_water)
median_water <- median(inc_data$mass_water)
sd_water <- sd(inc_data$mass_water)

## MERGE DATA FOR GRAIN SIZE MODEL ####

all_data <- left_join(inc_data, grain_all, by = c("Sample_ID")) %>% 
  left_join(moisture_clean, by = c("Sample_ID")) %>% 
  left_join(mean_ssa, by = c("Sample_ID")) %>% 
  dplyr::select(c(Sample_Name.x, Sample_ID, Initial_Water_mass_g, Final_Water_mass_g, Dry_Sediment_Mass_g, mass_sed_water, mass_water, average_ssa, Percent_Clay, Percent_Silt, Percent_Fine_Sand, Percent_Med_Sand, Percent_Coarse_Sand, Percent_Tot_Sand,  percent_water_content_dry, average_grav)) %>% 
  distinct(Sample_Name.x, .keep_all = TRUE) %>% 
  rename(Sample_Name = Sample_Name.x)

all_data$Percent_Clay <- as.numeric(all_data$Percent_Clay)

all_data$Percent_Silt <- as.numeric(all_data$Percent_Silt)

all_data$Percent_Fine_Sand <- as.numeric(all_data$Percent_Fine_Sand)

all_data$Percent_Med_Sand <- as.numeric(all_data$Percent_Med_Sand)

all_data$Percent_Coarse_Sand <- as.numeric(all_data$Percent_Coarse_Sand)

all_data$Percent_Tot_Sand <- as.numeric(all_data$Percent_Tot_Sand)


all_data <- all_data %>% 
  mutate(Percent_Mud = Percent_Clay + Percent_Silt) %>% 
  mutate(log_ssa = log10(average_ssa)) %>% 
  mutate(log_mud = log10(Percent_Mud))

## MLM using GRAIN SIZE and MASS WATER ####

#R2 0.5045 - don't have mass_sed_water for all
model_sed <- lm(mass_water ~ Dry_Sediment_Mass_g + Percent_Clay + Percent_Silt, data = all_data)

summary(model_sed)
summary(model_sed)$coefficient

model_residuals_sed = model_sed$residuals

hist(model_residuals_sed)
qqnorm(model_residuals_sed)
qqline(model_residuals_sed)

#R2 0.3994
model <- lm(mass_water ~ Percent_Clay + Percent_Silt + Percent_Coarse_Sand + Percent_Med_Sand + Percent_Fine_Sand + average_ssa, data = all_data)

summary(model)
summary(model)$coefficient

model_residuals = model$residuals

hist(model_residuals)
qqnorm(model_residuals)
qqline(model_residuals)

## R2 0.3866
model.5 <- lm(mass_water ~ Percent_Clay + Percent_Silt + Percent_Coarse_Sand + average_ssa, data = all_data)

summary(model.5)

## R2 0.2612
model_1.5 <- lm(mass_water ~  Percent_Coarse_Sand + average_ssa, data = all_data)

summary(model_1.5)

## ssa, % mud, R2 0.2043

model2 <- lm(mass_water ~  Percent_Mud + average_ssa, data = all_data)

summary(model2)
confint(model2)

sigma(model2)/mean(all_data$mass_water)


model2.5 <- lm(mass_water ~  log_mud + log_ssa, data = all_data)

summary(model2.5)
confint(model2)

#R2 0.1532
model3 <- lm(mass_water ~  average_ssa, data = all_data)

summary(model3)

## 50.90 - 0.17*Percent_Clay + 0.03*Percent_Silt - 0.007*Coarse_Sand + 0.006*Percent_Med_Sand + 0.021*average_ssa
# -4, +2% error in Mass Water Estimate

#57.34 -(0.42*Dry_Sediment_Mass_g) - (0.073*Percent_Clay) + (0.011*Percent_Silt))
# -3, +3% error in Mass Water Estimate

model_data <- all_data %>% 
  #mutate(modelled_water = 50.90 - (0.017*Percent_Clay) + (0.03*Percent_Silt) - (0.007*Percent_Coarse_Sand) + (0.006*Percent_Med_Sand) + (0.021*average_ssa)) %>% 
  mutate(modelled_water = 57.34 -(0.42*Dry_Sediment_Mass_g) - (0.073*Percent_Clay) + (0.011*Percent_Silt)) %>% 
  select(c(Sample_Name, Sample_ID, mass_water, modelled_water, Dry_Sediment_Mass_g)) %>% 
  mutate(percent_error = ((mass_water - modelled_water)/mass_water)*100) %>% 
  mutate(mean_water_error = ((mass_water - mean_water)/mass_water)*100) %>% 
  mutate(median_water_error = ((mass_water - median_water)/mass_water)*100)

## Resp MODEL

all_respiration <- read.csv("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/rates/ReadyForBoye/ECA_Sediment_Incubations_Respiration_Rates_ReadyForBoye_2023-11-08.csv")

all_respiration <- all_respiration %>% 
  select(c(Sample_Name, Respiration_Rate_mg_DO_per_L_per_H)) %>% separate(Sample_Name, c("EC", "Site", "INC"), sep = "_", remove = FALSE) %>% 
  unite(Sample_ID, c(EC:Site), sep = "_", remove = TRUE)

#leads to same %error range as modelled water, range in respiration rates within ~1% of total rates 

model_resp_data <- left_join(all_respiration, model_data, by = c( "Sample_Name", "Sample_ID")) %>% 
  left_join(grain_all, by = "Sample_ID") %>% 
  left_join(og_inc, by = c("Sample_Name")) %>% 
  left_join(mean_ssa, by = "Sample_ID") %>% 
  dplyr::select(-c(Dry_Sediment_Mass_g.x)) %>% 
  rename(Dry_Sediment_Mass_g = Dry_Sediment_Mass_g.y) %>% 
  mutate(Percent_Clay = as.numeric(Percent_Clay)) %>% 
  mutate(Percent_Silt = as.numeric(Percent_Silt)) %>% 
  mutate(Percent_Coarse_Sand = as.numeric(Percent_Coarse_Sand)) %>% 
  mutate(Percent_Med_Sand = as.numeric(Percent_Med_Sand)) %>% 
  mutate(mass_water = if_else(is.na(modelled_water), (57.34 -(0.42*Dry_Sediment_Mass_g) - (0.073*Percent_Clay) + (0.011*Percent_Silt) ), mass_water)) %>% 
  mutate(flag = if_else(is.na(modelled_water), "Modelled Water Mass", "N/A")) %>%
  mutate(Respiration_Rate_mg_DO_per_kg_per_H = if_else(Respiration_Rate_mg_DO_per_L_per_H > -9999, Respiration_Rate_mg_DO_per_L_per_H * (mass_water/Dry_Sediment_Mass_g), -9999)) %>% 
  mutate(rate_mg_per_L_per_kg_model = if_else(Respiration_Rate_mg_DO_per_L_per_H > -9999,  Respiration_Rate_mg_DO_per_L_per_H * (modelled_water/Dry_Sediment_Mass_g), -9999)) %>% 
  mutate(percent_error_resp = ((Respiration_Rate_mg_DO_per_kg_per_H - rate_mg_per_L_per_kg_model)/Respiration_Rate_mg_DO_per_kg_per_H)*100) %>% 
  mutate(range_model = Respiration_Rate_mg_DO_per_kg_per_H - rate_mg_per_L_per_kg_model)


respiration_rates <- model_resp_data %>% 
  dplyr::select(c(Sample_Name, Respiration_Rate_mg_DO_per_L_per_H, Respiration_Rate_mg_DO_per_kg_per_H,  mass_water, modelled_water, Dry_Sediment_Mass_g, flag)) %>%
  mutate(across(c("Respiration_Rate_mg_DO_per_kg_per_H", "mass_water", "modelled_water"), round, 3))

write.csv(respiration_rates,paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/rates/ReadyForBoye/ECA_Sediment_Incubations_mg_kg_rates_",pnnl.user,"_on_",Sys.Date(),".csv"), row.names = F)  

solid <- respiration_rates %>% 
  dplyr::select(c(Sample_Name, mass_water, Dry_Sediment_Mass_g, flag))

write.csv(solid,paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/INC/03_ProcessedData/ECA_Sediment_Incubations_solid_solution_",pnnl.user,"_on_",Sys.Date(),".csv"), row.names = F)  
