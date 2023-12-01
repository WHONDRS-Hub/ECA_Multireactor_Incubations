##ECA Data Analysis 

#### Read in Data ####

#read in libraries

library(lubridate);library(writexl);library(raster);library(tidyverse);library(devtools)
library(readxl)
library(corrplot)
library(corrr)
library(vegan)
library(FactoMineR)
library(factoextra)
library(ggpmisc)

rm(list=ls());graphics.off()

# Set working directory to data file
#Example:
pnnl.user = 'laan208'

# choose file dates to read in 

effect.date = '2023-11-08'
respiration.date = '2023-11-08'
removed.respiration.date = '2023-11-13'
mg.kg.respiration.date = '2023-12-01'
grav.date = '2023-12-01'

#Read in all data
setwd(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/"))

#effect size - change date to most recent
effect_size <- read_csv(paste0("Optode multi reactor/Optode_multi_reactor_incubation/rates/ReadyForBoye/ECA_Effect_Size_ReadyForBoye_",effect.date,".csv"))


#Respiration rates with removals from dist matrix to calculate effect size 

#respiration <- read.csv(paste0("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/rates/ReadyForBoye/ECA_Sediment_Incubations_Removed_Respiration_Rates_",removed.respiration.date,".csv"))

#all_respiration <- read.csv(paste0("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/rates/ReadyForBoye/ECA_Sediment_Incubations_Respiration_Rates_ReadyForBoye_",respiration.date,".csv"))

all_respiration <- read.csv(paste0("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/rates/ReadyForBoye/ECA_Sediment_Incubations_mg_kg_rates_laan208_on_",mg.kg.respiration.date,".csv"))

all_respiration <- all_respiration %>% 
  dplyr::select(c(Sample_Name, Respiration_Rate_mg_DO_per_L_per_H, Respiration_Rate_mg_DO_per_kg_per_H))

#ECA Iron
iron <- read_csv(paste0("Fe/03_ProcessedData/EC_SFE_ReadyForBoye_06-29-2023.csv"))

iron <- iron %>% 
  dplyr::select(c(sample_label, Fe_mg_per_L))

iron_ex <- read_csv(paste0("Fe/03_ProcessedData/20230711_Data_Processed_SFE_ECA_EC_1-153/20230711_Data_Processed_SFE_ECA_EC_1-153.csv"))

iron_ex <- iron_ex %>% 
  rename(Fe_mg_per_L = mg_Fe_per_L) %>%
  filter(is.na(flag)) %>% 
  dplyr::select(c(sample_label, Fe_mg_per_L))

iron <- rbind(iron, iron_ex)

rbind()


#ICON Grain Size
grain <- read_csv("C:/GitHub/ECA_Multireactor_Incubations/Data/v2_CM_SSS_Sediment_Grain_Size.csv")

ssa <- read_csv(paste0("C:/GitHub/ECA_Multireactor_Incubations/Data/eca_ssa_predatapackage.csv"))


#All incubation pH, SpC, temp
chemistry <- read_csv("INC/03_ProcessedData/SpC_pH_Temp.csv")


#Gravimetric Moisture

grav_inc <- read.csv(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/INC/03_ProcessedData/ECA_Drying_Masses_Summary_merged_by_laan208_on_",grav.date,".csv"))
