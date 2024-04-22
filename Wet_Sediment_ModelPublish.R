## ECA Model for Wet/Dry Sediment Masses
# Maggi Laan, maggi.laan@pnnl.gov
# 2024-01-19

#### To update upon publishing: 
#Boye Ready for Drying Masses on Github
#Change working directory, remove PNNL user
#Update Typical Code to have Incubation Water Mass


# Read in Libraries
library(ggplot2)
library(tidyverse)

rm(list=ls());graphics.off()

# Set PNNL User
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

## KNOWN Incubation Water Masses and Dry Sediment Mass ####

# Gives dry sediment mass calculated from moisture tins, and initial/final water masses

inc <- read.csv(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/INC/03_ProcessedData/ECA_Drying_Masses_Summary_ReadyForBoye_01-17-2024.csv"))

# Remove -9999 values (EC_023), remove samples that were modelled to get samples with measured Incubation Water Masses (110 total samples)

inc_model_data <- inc %>% 
  filter(grepl("N/A", flag)) %>% 
  separate(Sample_Name, into = c("EC", "Kit", "REP"), sep = "_", remove = FALSE) %>% 
  unite("Parent_ID", EC:Kit, sep = "_") %>% 
  rename(Incubation_Water_Mass_g_Actual = Incubation_Water_Mass_g)

#ggplot(inc_model_data, aes (x = Incubation_Water_Mass_g))+
#geom_histogram()

# Incubation_Water_Mass_g from ~49 mL to 53 mL
#mean_water <- mean(inc_model_data$Incubation_Water_Mass_g) #50.8
#median_water <- median(inc_model_data$Incubation_Water_Mass_g) #50.7
#sd_water <- sd(inc_model_data$Incubation_Water_Mass_g) # 0.85

## Specific Surface Area Data
ssa <- read.csv("C:/GitHub/ECA_Multireactor_Incubations/Data/CM_SSS_Sediment_Specific_Surface_Area.csv", skip = 2, header = TRUE)

#Chose only ICON SSA 
ssa_clean <- ssa %>% 
  filter(!row_number() %in% c(1:11)) %>% 
  dplyr::select(-c(Field_Name, Material, IGSN)) %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "CM", "EC")) %>% 
  filter(!grepl("SSS", Sample_Name)) %>% 
  separate(Sample_Name, into = c("EC", "Kit", "GRN"), sep = "_") %>% 
  unite("Parent_ID", EC:Kit, sep = "_") %>% 
  dplyr::select(-c(GRN)) %>% 
  filter(!grepl("NA", Parent_ID)) %>% 
  filter(!grepl("-9999", Specific_Surface_Area_m2_per_g)) %>% 
  mutate(Specific_Surface_Area_m2_per_g = as.numeric(Specific_Surface_Area_m2_per_g))

#Take means of SSA measurements
mean_ssa <- ssa_clean %>% 
  group_by(Parent_ID) %>% 
  mutate(average_ssa = mean(Specific_Surface_Area_m2_per_g)) %>% 
  distinct(Parent_ID, .keep_all = TRUE) %>% 
  dplyr::select(Parent_ID, average_ssa)

## MERGE DATA  ####

all_data <- left_join(inc_model_data, grain_all, by = c("Parent_ID")) %>% 
  left_join(mean_ssa, by = c("Parent_ID")) %>% 
  dplyr::select(c(Sample_Name, Parent_ID, Initial_Water_Mass_g, Final_Water_Mass_g, Dry_Sediment_Mass_g, Incubation_Water_Mass_g_Actual, average_ssa, Percent_Clay, Percent_Silt, Percent_Fine_Sand, Percent_Med_Sand, Percent_Coarse_Sand, Percent_Tot_Sand)) %>% 
  mutate(Percent_Mud = Percent_Clay + Percent_Silt) %>% 
  mutate(log_ssa = log10(average_ssa)) %>% 
  mutate(log_mud = log10(Percent_Mud))

## MLM using GRAIN SIZE and MASS WATER ####

#Adj.R2 0.4893 -  %clay/silt both 0.1 (%Silt, %Clay sometimes sig. in other models), removing grn size variables R2 = 0.4759
model_sed <- lm(Incubation_Water_Mass_g_Actual ~ Dry_Sediment_Mass_g + Percent_Clay + Percent_Silt, data = all_data)

summary(model_sed)
summary(model_sed)$coefficient

#Check model residuals
model_residuals_sed = model_sed$residuals

hist(model_residuals_sed)
qqnorm(model_residuals_sed)
qqline(model_residuals_sed)

#57.34 -(0.42*Dry_Sediment_Mass_g) - (0.073*Percent_Clay) + (0.011*Percent_Silt))
# -3, +3% error in Mass Water Estimate using model 
# using mean data -3.2 - +3.8% 
# using median data -3 - +4.0%

model_data <- all_data %>% 
  mutate(Incubation_Water_Mass_g_Modelled = 57.34 -(0.42*Dry_Sediment_Mass_g) - (0.073*Percent_Clay) + (0.011*Percent_Silt)) %>% 
  select(c(Sample_Name, Parent_ID, Incubation_Water_Mass_g_Actual, Incubation_Water_Mass_g_Modelled)) %>% 
  mutate(percent_error = ((Incubation_Water_Mass_g_Actual - Incubation_Water_Mass_g_Modelled)/Incubation_Water_Mass_g_Actual)*100) #%>% 
 # mutate(mean_water_error = ((mass_water - mean_water)/mass_water)*100) %>% 
 # mutate(median_water_error = ((mass_water - median_water)/mass_water)*100)

## Final Solid:Solution Ratio
solid <- left_join(inc, model_data, by = "Sample_Name") %>% 
  separate(Sample_Name, into = c("EC", "Kit", "REP"), sep = "_", remove = FALSE) %>% 
  unite("Parent_ID", EC:Kit, sep = "_") %>% 
  left_join(grain_all, by = "Parent_ID") %>% 
  dplyr::select(c(Sample_Name, Initial_Water_Mass_g, Final_Water_Mass_g, Dry_Sediment_Mass_g, Incubation_Water_Mass_g, Incubation_Water_Mass_g_Actual, Incubation_Water_Mass_g_Modelled, Percent_Clay, Percent_Silt, Methods_Deviation.x)) %>%
  rename(Methods_Deviation = Methods_Deviation.x) %>% 
  mutate(Incubation_Water_Mass_g_Modelled = if_else(is.na(Incubation_Water_Mass_g_Actual), 57.34 -(0.42*Dry_Sediment_Mass_g) - (0.073*Percent_Clay) + (0.011*Percent_Silt), -9999)) %>% 
  mutate(flag = if_else(is.na(Incubation_Water_Mass_g_Actual), "Modelled Water Mass", "N/A")) %>%
   mutate(Incubation_Water_Mass_g_Modelled = round(Incubation_Water_Mass_g_Modelled, 2)) %>% 
  mutate(Incubation_Water_Mass_g = if_else(flag == "N/A", Incubation_Water_Mass_g_Actual, Incubation_Water_Mass_g_Modelled)) %>% 
  mutate(Incubation_Water_Mass_g = if_else(Incubation_Water_Mass_g > 1000, -9999, Incubation_Water_Mass_g)) %>% 
mutate(Incubation_Water_Mass_g = round(Incubation_Water_Mass_g, 2)) %>% 
  mutate(Methods_Deviation = if_else(Methods_Deviation != "N/A" & flag == "Modelled Water Mass", paste(Methods_Deviation, "; INC_QA_005"), Methods_Deviation)) %>% 
  mutate(Methods_Deviation = if_else(Methods_Deviation == "N/A" & flag == "Modelled Water Mass", "INC_QA_005", Methods_Deviation)) %>% 
  select(c(Sample_Name, Initial_Water_Mass_g, Final_Water_Mass_g, Dry_Sediment_Mass_g, Incubation_Water_Mass_g,Methods_Deviation))

#write.csv(solid,paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/INC/03_ProcessedData/ECA_Drying_Masses_Summary_ReadyForBoye_01-29-2024.csv"), row.names = F)  

## EV Incubation Water Masses ####

inc <- read_xlsx(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/INC/01_RawData/2022_Data_Raw_INC_Tube_Weight_ECA_EC.xlsx")) %>% 
  filter(grepl("EV", SampleID)) 

mean_tube_g = inc %>%
  filter(INC_tube_50ml_empty_g != -9999) %>% 
  mutate(mean_tube_g = mean(INC_tube_50ml_empty_g)) %>% 
  right_join(inc, by = "SampleID") %>% 
  fill(mean_tube_g, .direction = c("updown")) %>% 
  mutate(INC_tube_50ml_empty_g = ifelse(is.na(INC_tube_50ml_empty_g.x), mean_tube_g, INC_tube_50ml_empty_g.x)) %>% 
  select(c(SampleID, INC_tube_50ml_empty_g, `INC_tube_50mL_wet_Sediment_with_Water_g (After Subsampling for SFE + SpC).y`))

#calculate mass_sed_water (Mass of sediment + water) by subtracting the empty 50 mL tube weight from the total weight. Add 6 g (Mass of slurry removed for pH, SpC, Fe measurements) 

inc_clean <- mean_tube_g %>% 
  mutate(mass_sed_water = (`INC_tube_50mL_wet_Sediment_with_Water_g (After Subsampling for SFE + SpC).y` - INC_tube_50ml_empty_g) + 6) %>% 
  mutate(Sample_Name = SampleID) %>% 
  dplyr::select(-c(SampleID)) %>% 
  #na.omit(mass_sed_water) %>% 
  relocate(Sample_Name, .before = INC_tube_50ml_empty_g) %>% 
  rename(`INC_tube_50mL_wet_Sediment_with_Water_g (After Subsampling for SFE + SpC)` = `INC_tube_50mL_wet_Sediment_with_Water_g (After Subsampling for SFE + SpC).y`)

og_inc <- read.csv(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/INC/03_ProcessedData/EV_Drying_Masses_Summary_merged_by_laan208_on_2024-04-19.csv"))

## MERGE DATA ####

## Incubation weights

# calculate mass water in final samples by subtracting Mass of sediment + water - dry sediment mass 

#dropping NA's in mass water removes samples that we don't have final incubation weights for (110 total masses)

inc_data <- left_join(og_inc, inc_clean, by = c("Sample_Name")) %>% 
  mutate(Incubation_Water_Mass_g = mass_sed_water - Dry_Sediment_Mass_g) %>% 
  mutate(Incubation_Water_Mass_g = if_else(Incubation_Water_Mass_g < 0, -9999, Incubation_Water_Mass_g)) %>% 
  mutate(Incubation_Water_Mass_g = if_else(is.na(Incubation_Water_Mass_g), -9999, Incubation_Water_Mass_g)) %>% 
  mutate(Incubation_Water_Mass_g = round(Incubation_Water_Mass_g, 2)) %>% 
  select(c(Sample_Name, Initial_Water_Mass_g, Final_Water_Mass_g, Dry_Sediment_Mass_g, Incubation_Water_Mass_g))




## Other Models Tried ####
# model <- lm(mass_water ~ Dry_Sediment_Mass_g + Percent_Mud + Percent_Tot_Sand, data = all_data)
# 
# summary(model)
# summary(model)$coefficient
# 
# #R2 0.3653
# model <- lm(mass_water ~ Percent_Clay + Percent_Silt + Percent_Coarse_Sand + Percent_Med_Sand + Percent_Fine_Sand + average_ssa, data = all_data)
# 
# summary(model)
# summary(model)$coefficient
# 
# model_residuals = model$residuals
# 
# hist(model_residuals)
# qqnorm(model_residuals)
# qqline(model_residuals)
# 
# ## R2 0.3866
# model.5 <- lm(mass_water ~ Percent_Clay + Percent_Silt + Percent_Coarse_Sand + average_ssa, data = all_data)
# 
# summary(model.5)
# 
# ## R2 0.2612
# model_1.5 <- lm(mass_water ~  Percent_Coarse_Sand + average_ssa, data = all_data)
# 
# summary(model_1.5)
# 
# ## ssa, % mud, R2 0.2043
# 
# model2 <- lm(mass_water ~  Percent_Mud + average_ssa, data = all_data)
# 
# summary(model2)
# 
# 
# model2.5 <- lm(mass_water ~  log_mud + log_ssa, data = all_data)
# 
# summary(model2.5)
# confint(model2)
