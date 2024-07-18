## This script only written for ECA moisture tins and samples incubated in 40 mL vials (although jars included in csv file, they were not analyzed for outliers at time of this analysis)

##### Load Library #####
library(tidyverse)
library(dplyr)
library(readxl)
library(EnvStats)

##### Load data ######
rm(list=ls());graphics.off()

pnnl.user = 'laan208'
date = '01-17-2024'
study.code = 'CM'

## Load Data

# gravimetric moisture calculated from tins

eca <- read_csv(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/MOI/03_ProcessedData/",study.code,"_MOI_ReadyForBoye_", date, ".csv"))

# compare to ICON gravimetric water masses from sediment incubations
icon <- read_csv(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ICON_ModEx_SSS/Subsampling/Sediment_Mass_Water_Calculations/ICON_Sediment_Incubations_Sediment_g_and_Water_mL_merged_by_laan208_on_2023-10-16.csv"))

# calculate gravimetric moisture from icon
icon_grav <- icon %>% 
  mutate(grav_icon = (((Wet_Sediment_Mass_g - Dry_Sediment_Mass_g)/Dry_Sediment_Mass_g)*100)) %>% 
  relocate(grav_icon, .after = Sample_ID) %>% 
  separate(Sample_ID, c("EC", "Site", "Rep"), sep = "_") %>% 
  filter(!grepl("-4|-5|-6",Rep)) %>% 
  separate(Rep, c("INC", "Rep"), sep = "-") %>% 
  group_by(Site) %>% 
  mutate(average_icon = mean(grav_icon)) %>% 
  mutate(cv_icon = (sd(grav_icon)/mean(grav_icon))*100) %>% 
  ungroup() %>% 
  relocate(average_icon, .after = grav_icon) %>% 
  relocate(cv_icon, .after = average_icon) %>% 
  dplyr::select(-c(EC, INC, Methods_Deviation, Notes, Wet_Sediment_Mass_g, Dry_Sediment_Mass_g, Water_Mass_g, Wet_Sediment_mL))

# Calculate Mean gravimetric moisture from ECA and remove previously defined outliers:
  #EC_057_MOI-3 spilled before weighing
  #EC_084_MOI-2 typo? negative value

eca_grav <- eca %>% 
  filter(Sample_Name != "CM_057_MOI-3") %>% 
  filter(Sample_Name != "CM_084_MOI-2") %>% 
  separate(Sample_Name, c("Sample_Site", "Rep"), sep = "-", remove = FALSE) %>% 
  separate(Sample_Site, c("EC", "Site", "MOI"), sep = "_") %>% 
  group_by(Site) %>% 
  mutate(average_eca = mean(Gravimetric_Moisture)) %>% 
  mutate(cv_eca = (sd(Gravimetric_Moisture)/mean(Gravimetric_Moisture))*100) %>%   select(c("Sample_Name", "Site", "Rep", "grav_eca", "average_eca", "cv_eca"))

# Merge icon and eca gravimetric moistures

icon_eca <- merge(eca_grav, icon_grav, by = c("Site", "Rep"), all = TRUE) %>% 
  drop_na(average_eca) %>% # remove CM_023
  mutate(grav_moi_diff = (average_eca - average_icon)) %>% # look at difference in icon vs. eca grav. moisture
  relocate(grav_moi_diff, .after = Site) %>% 
  unite(Sample_Name, c("Site", "Rep"), sep = "_", remove = FALSE)

ggplot(icon_eca)+
  geom_point(aes(x = Site, y = Gravimetric_Moisture, color = "ECA"))+
  geom_point(aes(x = Site, y = grav_icon, color = "ICON"))+
  scale_fill_discrete(labels = c("ECA", "ICON"))


### ABOVE 15% CV (based on histogram)

#Site 082 - EC_082_MOI-3: Large outlier, confirmed with ICON (+/- 10%)

#Site 028 - not confirmed, one large outlier

#Site 106 - EC_106_MOI-2: Large outlier, confirmed with ICON (+/- 10%)

#Site 060 - EC_060_MOI-3: Large outlier, confirmed with ICON (+/- 10%)

#Site_035 - not confirmed - keep all?

#Site_012 - not confirmed - keep all?

#Site 062 - average close to ICON, keep all

#Site 068 - EC_068_MOI-1: Large outlier, confirmed with ICON (+/- 15%)

#Site 097 - average close to ICON

#Site 064 - not confirmed - keep all?

#Site 002 - close, don't remove anything

# Remove samples that could be confirmed as outliers using icon gravimetric moisture, calculate new average gravimetric moisture and CV
outliers_eca <- eca_grav %>%
  filter(!grepl("082_MOI-3", Sample_Name)) %>% 
  filter(!grepl("106_MOI-2", Sample_Name)) %>% 
  filter(!grepl("060_MOI-3", Sample_Name)) %>% 
  filter(!grepl("028_MOI-3", Sample_Name)) %>% 
  filter(!grepl("068_MOI-1", Sample_Name)) %>% 
  group_by(Site) %>% 
  mutate(new_average_eca = mean(Gravimetric_Moisture)) %>% 
  mutate(new_cv_eca = (sd(Gravimetric_Moisture)/mean(Gravimetric_Moisture))*100) %>% 
  ungroup() %>% 
  relocate(new_average_eca, .after = average_eca) %>% 
  relocate(new_cv_eca, .after = cv_eca) 

#make data frame of flagged outliers 

flag <- left_join(eca, outliers_eca, by = "Sample_Name") %>% 
  mutate(flag = if_else(grepl("CM_023", Sample_Name), "Missing Wet Weight",if_else(is.na(new_average_eca), "removed outlier", "NA"))) %>% 
  rename(Gravimetric_Moisture = Gravimetric_Moisture.x) %>% 
  select(c(Sample_Name, flag))

# Remake dataframe and flag outlier samples
clean_eca <- left_join(eca, flag)

write.csv(clean_eca, paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/MOI/03_ProcessedData/CM_MOI_ReadyForBoye_01-19-2024.csv"), row.names = F)
 
