library(ggplot2)
library(tidyverse)
library(dplyr)
library(readxl)

rm(list=ls());graphics.off()

pnnl.user = 'laan208'

setwd("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files")

## Pull in Data for Model

## GRAIN SIZE ####

grain <- read.csv("C:/Github/ECA_Multireactor_Incubations/Data/v3_CM_SSS_Sediment_Grain_Size.csv", skip = 2, header = TRUE)

grain_all <- grain %>% 
  filter(!row_number() %in% c(1:11)) %>% 
  dplyr::select(-c(Field_Name, Material)) %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "CM", "EC")) %>% 
  filter(!grepl("SSS", Sample_Name)) %>% 
  separate(Sample_Name, into = c("EC", "Kit", "GRN"), sep = "_") %>% 
  unite("Parent_ID", EC:Kit, sep = "_") %>% 
  dplyr::select(-c(GRN, IGSN)) %>% 
  filter(!grepl("NA", Parent_ID))%>% 
  mutate(Percent_Clay = as.numeric(Percent_Clay)) %>% 
  mutate(Percent_Silt = as.numeric(Percent_Silt)) %>% 
  mutate(Percent_Fine_Sand = as.numeric(Percent_Fine_Sand)) %>% 
  mutate(Percent_Med_Sand = as.numeric(Percent_Med_Sand)) %>% 
  mutate(Percent_Coarse_Sand = as.numeric(Percent_Coarse_Sand)) %>% 
  mutate(Percent_Tot_Sand = as.numeric(Percent_Tot_Sand))

## KNOWN INCUBATION WEIGHTS ####

inc <- read_xlsx(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/INC/01_RawData/2022_Data_Raw_INC_Tube_Weight_ECA_EC.xlsx"))

#Fix typo in samples. Weights > 80 cannot exist, so subtract 10 g

#calculate mass_sed_water (Mass of sediment + water) by subtracting the empty 50 mL tube weight from the total weight. Add 6 g (Mass of slurry removed for pH, SpC, Fe measurements) 

inc_clean <- inc %>% 
  mutate(`INC_tube_50mL_wet_Sediment_with_Water_g (After Subsampling for SFE + SpC)`= 
           ifelse(`INC_tube_50mL_wet_Sediment_with_Water_g (After Subsampling for SFE + SpC)` >= 80, `INC_tube_50mL_wet_Sediment_with_Water_g (After Subsampling for SFE + SpC)` - 10,
                  `INC_tube_50mL_wet_Sediment_with_Water_g (After Subsampling for SFE + SpC)`)) %>% 
  mutate(mass_sed_water = (`INC_tube_50mL_wet_Sediment_with_Water_g (After Subsampling for SFE + SpC)` - INC_tube_50ml_empty_g) + 6) %>% 
  mutate(Sample_Name = SampleID) %>% 
  dplyr::select(-c(SampleID)) %>% 
  filter(grepl("EC", Sample_Name)) %>% 
  #na.omit(mass_sed_water) %>% 
  relocate(Sample_Name, .before = INC_tube_50ml_empty_g)

## Specific Surface Area Data
ssa <- read.csv("C:/GitHub/ECA_Multireactor_Incubations/Data/CM_SSS_Sediment_Specific_Surface_Area.csv", skip = 2, header = TRUE)

#Chose only ICON SSA 
ssa_clean <- ssa %>% 
  filter(!row_number() %in% c(1:11)) %>% 
  dplyr::select(-c(Field_Name, Material)) %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "CM", "EC")) %>% 
  filter(!grepl("SSS", Sample_Name)) %>% 
  separate(Sample_Name, into = c("EC", "Kit", "GRN"), sep = "_") %>% 
  unite("Parent_ID", EC:Kit, sep = "_") %>% 
  dplyr::select(-c(GRN, IGSN)) %>% 
  filter(!grepl("NA", Parent_ID)) %>% 
  filter(!grepl("-9999", Specific_Surface_Area_m2_per_g)) %>% 
  mutate(Specific_Surface_Area_m2_per_g = as.numeric(Specific_Surface_Area_m2_per_g))

#Take means of SSA measurements
mean_ssa <- ssa_clean %>% 
  group_by(Parent_ID) %>% 
  mutate(average_ssa = mean(Specific_Surface_Area_m2_per_g)) %>% 
  distinct(Parent_ID, .keep_all = TRUE) %>% 
  dplyr::select(Parent_ID, average_ssa)

## WEIGHTS from 21 day incubation  ####

# Gives dry sediment mass calculated from moisture tins, and initial/final water masses

og_inc <- read.csv(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/INC/03_ProcessedData/ECA_Drying_Masses_Summary_ReadyForBoye_2024-01-04.csv"))

## MERGE DATA ####

## Incubation weights

# calculate mass water in final samples by subtracting Mass of sediment + water - dry sediment mass 

#dropping NA's in mass water removes samples that we don't have final incubation weights for (110 total masses)

inc_data <- left_join(og_inc, inc_clean, by = c("Sample_Name")) %>% 
  mutate(mass_water = mass_sed_water - Dry_Sediment_Mass_g) %>% 
  separate(Sample_Name, into = c("EC", "Kit", "REP"), sep = "_", remove = FALSE) %>% 
  unite("Parent_ID", EC:Kit, sep = "_") %>% 
  drop_na(mass_water) 

#ggplot(inc_data, aes (x = mass_water))+
  #geom_histogram()

# Final water masses from ~49 mL to 53 mL

mean_water <- mean(inc_data$mass_water) #50.8
median_water <- median(inc_data$mass_water) #50.7
sd_water <- sd(inc_data$mass_water) # 0.85

## MERGE DATA FOR GRAIN SIZE MODEL ####

all_data <- left_join(inc_data, grain_all, by = c("Parent_ID")) %>% 
  left_join(mean_ssa, by = c("Parent_ID")) %>% 
  dplyr::select(c(Sample_Name, Parent_ID, Initial_Water_Mass_g, Final_Water_Mass_g, Dry_Sediment_Mass_g, mass_sed_water, mass_water, average_ssa, Percent_Clay, Percent_Silt, Percent_Fine_Sand, Percent_Med_Sand, Percent_Coarse_Sand, Percent_Tot_Sand)) %>% 
  mutate(Percent_Mud = Percent_Clay + Percent_Silt) %>% 
  mutate(log_ssa = log10(average_ssa)) %>% 
  mutate(log_mud = log10(Percent_Mud))

## MLM using GRAIN SIZE and MASS WATER ####

#Adj.R2 0.4893 - don't have mass_sed_water for all, only dry_sed_mass sig., %clay/silt both 0.1 (%Silt, %Clay sometimes sig. in other models), removing grn size variables R2 = 0.4759
model_sed <- lm(mass_water ~ Dry_Sediment_Mass_g + Percent_Clay + Percent_Silt, data = all_data)

summary(model_sed)
summary(model_sed)$coefficient

model_residuals_sed = model_sed$residuals

hist(model_residuals_sed)
qqnorm(model_residuals_sed)
qqline(model_residuals_sed)

#57.34 -(0.42*Dry_Sediment_Mass_g) - (0.073*Percent_Clay) + (0.011*Percent_Silt))
# -3, +3% error in Mass Water Estimate using model 
# using mean data -3.2 - +3.8% 
# using median data -3 - +4.0%

model_data <- all_data %>% 
  mutate(modelled_water = 57.34 -(0.42*Dry_Sediment_Mass_g) - (0.073*Percent_Clay) + (0.011*Percent_Silt)) %>% 
  select(c(Sample_Name, Parent_ID, mass_water, modelled_water)) %>% 
  mutate(percent_error = ((mass_water - modelled_water)/mass_water)*100) %>% 
  mutate(mean_water_error = ((mass_water - mean_water)/mass_water)*100) %>% 
  mutate(median_water_error = ((mass_water - median_water)/mass_water)*100)

## Final Solid:Solution Ratio
solid <- left_join(og_inc, model_data, by = "Sample_Name") %>% 
  separate(Sample_Name, into = c("EC", "Kit", "REP"), sep = "_", remove = FALSE) %>% 
  unite("Parent_ID", EC:Kit, sep = "_") %>% 
  left_join(grain_all, by = "Parent_ID") %>% 
  dplyr::select(c(Sample_Name, Initial_Water_Mass_g, Final_Water_Mass_g, Dry_Sediment_Mass_g, mass_water, modelled_water, Percent_Clay, Percent_Silt, Methods_Deviation.x)) %>%
  rename(Methods_Deviation = Methods_Deviation.x) %>% 
  mutate(Incubation_Water_Mass_g = if_else(is.na(mass_water), 57.34 -(0.42*Dry_Sediment_Mass_g) - (0.073*Percent_Clay) + (0.011*Percent_Silt), mass_water)) %>% 
  mutate(flag = if_else(is.na(modelled_water), "Modelled Water Mass", "N/A")) %>% 
  mutate(Incubation_Water_Mass_g = if_else(Incubation_Water_Mass_g > 1000, -9999, Incubation_Water_Mass_g)) %>% 
  mutate(Incubation_Water_Mass_g = round(Incubation_Water_Mass_g, 2)) %>% 
  select(c(Sample_Name, Initial_Water_Mass_g, Final_Water_Mass_g, Dry_Sediment_Mass_g, Incubation_Water_Mass_g, Methods_Deviation, flag))

write.csv(solid,paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/INC/03_ProcessedData/ECA_Drying_Masses_Summary_ReadyForBoye_01-17-2024.csv"), row.names = F)  

## Resp MODEL

all_respiration <- read.csv("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/rates/ReadyForBoye/ECA_Sediment_Incubations_Respiration_Rates_ReadyForBoye_2023-11-08.csv")

all_respiration <- all_respiration %>% 
  select(c(Sample_Name, Respiration_Rate_mg_DO_per_L_per_H)) %>%
  rename(Sample_ID = Sample_Name) %>% 
  separate(Sample_ID, c("EC", "Site", "INC"), sep = "_", remove = FALSE) %>% 
  unite(Sample_Name, c(EC:Site), sep = "_", remove = TRUE)

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





## Other Models Tried ####
model <- lm(mass_water ~ Dry_Sediment_Mass_g + Percent_Mud + Percent_Tot_Sand, data = all_data)

summary(model)
summary(model)$coefficient

#R2 0.3653
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


model2.5 <- lm(mass_water ~  log_mud + log_ssa, data = all_data)

summary(model2.5)
confint(model2)


