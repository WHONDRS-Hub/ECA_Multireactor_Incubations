##### Load Library #####
library(tidyverse)
library(dplyr)
library(readxl)
library(EnvStats)

##### Load data ######
rm(list=ls());graphics.off()

pnnl.user = 'laan208'

## Comparisons to ICON DWT

eca <- read_csv(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/MOI/03_ProcessedData/EC_Moisture_Content_2022.csv"))

#icon <- read_csv(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ICON_ModEx_SSS/Subsampling/Sediment_Mass_Water_Calculations/ICON_Sediment_Incubations_Sediment_g_and_Water_mL_merged_by_laan208_on_2023-10-16.csv"))

## 21 day dry down

weights <- read_xlsx(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/INC/01_RawData/2022_Data_Raw_INC_ECA_EC.xlsx"))

wt_clean <- weights %>% 
  filter(`Jars or 40 mL vials` != "Jars")

icon_grav <- icon %>% 
  mutate(grav_icon = (((Wet_Sediment_Mass_g - Dry_Sediment_Mass_g)/Dry_Sediment_Mass_g)*100)) %>% 
  relocate(grav_icon, .after = Sample_ID) %>% 
  separate(Sample_ID, c("EC", "Site", "Rep"), sep = "_") %>% 
  filter(!grepl("-4",Rep)) %>% 
  filter(!grepl("-5",Rep)) %>% 
  filter(!grepl("-6",Rep)) %>% 
  separate(Rep, c("INC", "Rep"), sep = "-") %>% 
  group_by(Site) %>% 
  mutate(average_icon = mean(grav_icon)) %>% 
  mutate(cv_icon = (sd(grav_icon)/mean(grav_icon))*100) %>% 
  relocate(average_icon, .after = grav_icon) %>% 
  relocate(cv_icon, .after = average_icon) %>% 
  dplyr::select(-c(EC, INC, Methods_Deviation, Notes, Wet_Sediment_Mass_g, Dry_Sediment_Mass_g, Water_Mass_g, Wet_Sediment_mL))

##EC_057_MOI-3 spilled before weighing
##EC_084_MOI-2 typo? negative value

eca_grav <- eca %>% 
  filter(Sample_Name != "EC_057_MOI-3") %>% 
  filter(Sample_Name != "EC_084_MOI-2") %>% 
  separate(Sample_Name, c("Sample_Site", "Rep"), sep = "-", remove = FALSE) %>% 
  separate(Sample_Site, c("EC", "Site", "MOI"), sep = "_") %>% 
  mutate(grav_eca = ((wet_weight_g - true_dry_weight_g)/true_dry_weight_g)*100) %>% 
  group_by(Site) %>% 
  mutate(average_eca = mean(grav_eca)) %>% 
  mutate(cv_eca = (sd(grav_eca)/mean(grav_eca))*100) %>%   select(c("Sample_Name", "Site", "Rep", "grav_eca", "average_eca", "cv_eca"))

eca_cv <- eca_grav %>% 
  distinct(Site, .keep_all = TRUE)

ggplot(eca_cv, aes(x = cv_eca)) +
  geom_histogram()

icon_eca <- merge(eca_grav, icon_grav, by = c("Site", "Rep"), all = TRUE)

#028 not used for ECA

icon_eca <- icon_eca %>% 
  drop_na(average_eca) %>% 
  mutate(grav_moi_diff = (average_eca - average_icon)) %>% 
  relocate(grav_moi_diff, .after = Site) %>% 
  unite(Sample_Name, c("Site", "Rep"), sep = "_", remove = FALSE)

ggplot(icon_eca)+
  geom_point(aes(x = Site, y = grav_eca, color = "ECA"))+
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


outliers_eca <- eca_grav %>%
  filter(!grepl("082_MOI-3", Sample_Name)) %>% 
  filter(!grepl("106_MOI-2", Sample_Name)) %>% 
  filter(!grepl("060_MOI-3", Sample_Name)) %>% 
  filter(!grepl("028_MOI-3", Sample_Name)) %>% 
  filter(!grepl("068_MOI-1", Sample_Name)) %>% 
  group_by(Site) %>% 
  mutate(new_average_eca = mean(grav_eca)) %>% 
  mutate(new_cv_eca = (sd(grav_eca)/mean(grav_eca))*100) %>% 
  relocate(new_average_eca, .after = average_eca) %>% 
  relocate(new_cv_eca, .after = cv_eca) 

clean_eca <- left_join(eca, outliers_eca, by = "Sample_Name") %>% 
  mutate(flag = if_else(is.na(new_average_eca), "removed outlier", "NA")) %>% 
  mutate(flag = if_else(grepl("EC_023", Sample_Name), "Missing Wet Weight", flag)) %>% 
  mutate(Gravimetric_Moisture = ((wet_weight_g - true_dry_weight_g)/true_dry_weight_g)*100) %>% 
  dplyr::select(c(Sample_Name, wet_weight_g, true_dry_weight_g, percent_water_content_dry, percent_water_content_wet, Gravimetric_Moisture, flag)) %>% 
  rename(Wet_Sediment_Mass_g = wet_weight_g) %>% 
  rename(Dry_Sediment_Mass_g = true_dry_weight_g) %>% 
  mutate(Water_Mass_g = Wet_Sediment_Mass_g - Dry_Sediment_Mass_g) %>% 
  mutate(Gravimetric_Moisture = round(Gravimetric_Moisture, 2)) %>% 
  mutate(Methods_Deviation = "N/A") %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "EC", "CM")) %>% 
  dplyr::select(c(Sample_Name, Wet_Sediment_Mass_g, Dry_Sediment_Mass_g, Water_Mass_g, Gravimetric_Moisture, flag, Methods_Deviation))

write.csv(clean_eca, paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/MOI/03_ProcessedData/CM_MOI_ReadyForBoye_01-19-2024.csv"), row.names = F)
 
