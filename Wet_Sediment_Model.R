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
                  `INC_tube_50mL_wet_Sediment_with_Water_g (After Subsampling for SFE + SpC)`))

inc <- inc %>% 
  mutate(mass_sed_water = (`INC_tube_50mL_wet_Sediment_with_Water_g (After Subsampling for SFE + SpC)` - INC_tube_50ml_empty_g) + 6) %>% 
  mutate(Sample_Name = SampleID) %>% 
  dplyr::select(-c(SampleID)) %>% 
  filter(!grepl("EV", Sample_Name)) %>% 
  #na.omit(mass_sed_water) %>% 
  relocate(Sample_Name, .before = INC_tube_50ml_empty_g)

## Compare to SSA
ssa <- read.csv(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/SSA/03_ProcessedData/SSA_Means_2023.csv"))

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


## ORIGINAL INCUBATION WEIGHTS ####

og_inc <- read.csv(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/INC/03_ProcessedData/ECA_Drying_Masses_Summary_merged_by_laan208_on_2023-10-24.csv"))

#og_inc <- og_inc %>% 
  #distinct(Sample_Name, .keep_all = TRUE)

## MOISTURE TIN DATA ####

moisture <- read_csv(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/MOI/03_ProcessedData/EC_Moisture_Content_2022.csv"))

moisture <- moisture %>% 
  filter(Sample_Name != "EC_057_MOI-3") %>% 
  filter(Sample_Name != "EC_084_MOI-2") %>% 
  separate(Sample_Name, c("EC", "Site", "Rep"), sep = "_") %>% 
  group_by(Site) %>% 
  mutate(average = mean(percent_water_content_dry)) %>% 
  mutate(cv = sd(percent_water_content_dry)/mean(percent_water_content_dry))

## Remove EC_082_MOI-3: confirmed with ICON average
moisture_clean <- moisture %>% 
  filter(!(Site == "082" & Rep == "MOI-3")) %>% 
  unite("Sample_ID", EC:Site, sep = "_") 

moisture_average <- moisture_clean %>% 
  distinct(Sample_ID, .keep_all = TRUE)

## Regression of Moisture by Grain Size

moi_grn <- left_join(moisture_clean, grain_all, by = c("Sample_ID"))

moi_grn <- moi_grn %>% 
  drop_na(true_dry_weight_g) %>% 
  rename(average_perc_water = average)

moi_grn$Percent_Clay <- as.numeric(moi_grn$Percent_Clay)

moi_grn$Percent_Silt <- as.numeric(moi_grn$Percent_Silt)

ggplot(moi_grn, aes(x = average, y = (Percent_Clay + Percent_Silt))) +
  geom_point(aes(color = Sample_ID))

## MERGE DATA FOR MASS WATER####

inc_data <- left_join(og_inc, inc, by = c("Sample_Name"))


inc_data <- inc_data %>% 
  mutate(mass_water = mass_sed_water - Dry_Sediment_Mass_g) %>% separate(Sample_Name, into = c("EC", "Kit", "REP"), sep = "_", remove = FALSE) %>% 
  unite("Sample_ID", EC:Kit, sep = "_") %>% 
  drop_na(mass_water) 

ggplot(inc_data, aes (x = mass_water))+
  geom_histogram()

mean_water <- mean(inc_data$mass_water)
median_water <- median(inc_data$mass_water)
sd(inc_data$mass_water)

## MERGE DATA FOR GRAIN SIZE MODEL ####

all_data <- left_join(inc_data, grain_all, by = c("Sample_ID"))

all_data <- left_join(all_data, moisture_clean, by = c("Sample_ID"))

all_data <- left_join(all_data, ssa, by = c("Sample_ID"))

all_data <- all_data %>% 
  dplyr::select(c(Sample_Name, Sample_ID, Initial_Water_mass_g, Final_Water_mass_g, Dry_Sediment_Mass_g, mass_water, mean, Percent_Clay, Percent_Silt, Percent_Fine_Sand, Percent_Med_Sand, Percent_Coarse_Sand, Percent_Tot_Sand,  percent_water_content_dry, average)) %>% 
  rename(average_perc_water = average) %>% 
  distinct(Sample_Name, .keep_all = TRUE)

all_data$Percent_Clay <- as.numeric(all_data$Percent_Clay)

all_data$Percent_Silt <- as.numeric(all_data$Percent_Silt)

all_data$Percent_Fine_Sand <- as.numeric(all_data$Percent_Fine_Sand)

all_data$Percent_Med_Sand <- as.numeric(all_data$Percent_Med_Sand)

all_data$Percent_Coarse_Sand <- as.numeric(all_data$Percent_Coarse_Sand)

all_data$Percent_Tot_Sand <- as.numeric(all_data$Percent_Tot_Sand)


all_data <- all_data %>% 
  mutate(Percent_Mud = Percent_Clay + Percent_Silt)

## MLM using GRAIN SIZE and MASS WATER ####

model <- lm(mass_water ~ Percent_Mud + Percent_Coarse_Sand + Percent_Med_Sand + Percent_Fine_Sand + Dry_Sediment_Mass_g + average_perc_water, data = all_data)

summary(model)
summary(model)$coefficient

model_residuals = model$residuals

hist(model_residuals)
qqnorm(model_residuals)
qqline(model_residuals)

## dry wet sed. % mud, average grav water

model2 <- lm(mass_water ~  mean + Percent_Mud, data = all_data)

summary(model2)
confint(model2)

sigma(model2)/mean(all_data$mass_water)

## 50.2 + 0.00184*ssa + 0.0422*Percent_Mud
# +/- 4% error in Mass Water Estimate

model_data <- all_data %>% 
  mutate(modelled_water = 50.2 + (0.00184*mean) + (0.0422*Percent_Mud)) %>% 
  select(c(Sample_Name, mass_water, modelled_water)) %>% 
  mutate(percent_error = ((mass_water - modelled_water)/mass_water)*100)


## FE MODEL ####

fe_files <- list.files(path = "ECA/Fe/03_ProcessedData/", pattern = "SFE_ECA_EC", recursive = TRUE, full.names = T) %>% 
  lapply(read_csv) %>% 
  bind_rows %>% 
  rename(Sample_Name = Sample_ID) %>% 
  select(-c(...1, sample_label)) 

fe_files$Sample_Name <- str_replace(fe_files$Sample_Name, "SFE", "INC") 


model_verf <- left_join(all_data, fe_files, by = c("Sample_Name"))

## 50.2 + 0.00184*ssa + 0.0422*Percent_Mud
# +/- 4% error in Mass Water Estimate

model_sfe <- model_verf %>% 
  relocate(Technical, .after = Sample_Name) %>% 
  filter(model_verf$mg_Fe_per_L >= 0) %>% 
  mutate(Percent_Mud = Percent_Clay + Percent_Silt) %>% 
  select(c(Sample_Name, Technical, Dry_Sediment_Mass_g, mass_water, mean, Percent_Mud, mg_Fe_per_L, average_perc_water)) %>%
  mutate(mass_model_water = 50.2 + (0.00184*mean) + (0.0422*Percent_Mud)) %>%  
  relocate(mass_model_water, .after = mass_water) %>% 
  mutate(model_water_error = ((mass_model_water - mass_water)/mass_water)*100) %>% 
  relocate(model_water_error, .after = mass_model_water) %>% 
  mutate(actual_water_mg_kg = (mg_Fe_per_L * (mass_water/Dry_Sediment_Mass_g)*1000)) %>% 
  mutate(median_water_mg_kg = mg_Fe_per_L * (median_water/Dry_Sediment_Mass_g)*1000) %>% 
  mutate(median_water_error = ((median_water_mg_kg - actual_water_mg_kg)/actual_water_mg_kg)*100) %>% 
  mutate(mean_water_mg_kg = mg_Fe_per_L * (mean_water/Dry_Sediment_Mass_g)*1000) %>% 
  mutate(mean_water_error = ((mean_water_mg_kg - actual_water_mg_kg)/actual_water_mg_kg)*100) %>% 
  mutate(model_water_mg_kg = mg_Fe_per_L * (mass_model_water/Dry_Sediment_Mass_g)*1000) %>% 
    mutate(model_water_error_mg_kg = ((model_water_mg_kg - actual_water_mg_kg)/actual_water_mg_kg)*100)

hist(model_verf$mass_model_water)
hist(model_verf$mass_water)

ggplot(model_verf, aes(x = mass_water, y = mass_model_water))+
  geom_point() +
  geom_smooth()

### mean error: -7 to 3%
### median error: -7 to 3%
### model error: -3 to 3% (+/- 1.5 mL)

### issues from one day (EC_093 and EC_094 W1,2,3 samples where 80 written instead of 70 instead of the rest, typo?)
