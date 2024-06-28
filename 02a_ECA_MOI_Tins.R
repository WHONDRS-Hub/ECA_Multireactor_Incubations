##### Load Library #####
library(tidyverse)
library(readxl)

##### Load data ######
rm(list=ls());graphics.off()

pnnl.user = 'laan208'
study.code = 'EC'
year = '2022'

# Set working directory to data file
  #For Macbook Users
#data.in.path <- file.path("/Users/delg580/Library/CloudStorage/OneDrive-SharedLibraries-PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/MOI/01_RawData")

#data.out.path <- file.path("/Users/delg580/Library/CloudStorage/OneDrive-SharedLibraries-PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/MOI/03_ProcessedData")

#For Windows Users

data.out.path <- file.path("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/MOI/03_ProcessedData")

# Load Data
moisture <- read_excel(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/MOI/01_RawData/",year,"_Data_Raw_MOI_ECA_",study.code,".xlsx")) %>% janitor::clean_names()

# Subsetting data by wet weight
wet_weight <- moisture %>% 
  dplyr::select(sample_name, date, tare_weight_g, sample_weight_g) %>% 
  filter(!grepl("-9999" , tare_weight_g))

colnames(wet_weight) [4] <- "wet_weight_g"

# Subsetting data by dry weight
dry_weight <- moisture %>% 
  dplyr::select(sample_name, date, tare_weight_g, sample_weight_g) %>% 
  filter(grepl("-9999", tare_weight_g))

colnames(dry_weight) [4] <- "dry_weight_g"


dim(wet_weight); length(unique(wet_weight$sample_name)) #number of unique entries
dim(dry_weight); length(unique(dry_weight$sample_name)) #number of unique entries
length(which(wet_weight$sample_name %in% dry_weight$sample_name)) #which ID's in wet weight file are found in the dry weight file

# Merged wet and dry weights by sample name
merged_weights <- merge(x = wet_weight, y = dry_weight, by = "sample_name", all = TRUE) %>% janitor::clean_names()
dim(merged_weights)
head(merged_weights)
names(merged_weights)

# Calculate True dry weight 
merged_weights$true_dry_weight_g <- merged_weights$dry_weight_g - merged_weights$tare_weight_g_x 
range(merged_weights$true_dry_weight_g)

# Calculatee Gravimetric Water Content
merged_weights$percent_water_content_dry <- ((merged_weights$wet_weight_g - merged_weights$true_dry_weight_g)/merged_weights$true_dry_weight_g)*100


# Clean up data
cols_to_drop <- c ("date_x","date_y", "tare_weight_g_y","tare_weight_g_x", "dry_weight_g" )
merged_weights <- merged_weights[ , !(names(merged_weights) %in% cols_to_drop)]

range(merged_weights$percent_water_content_dry)

merged_weights$flag[is.na(merged_weights$wet_weight_g)] = "Missing Wet Weight"

merged_weights <- merged_weights %>% 
  separate(sample_name, c("EC", "Site", "Rep"), sep = "_")

# Add extra 0 to missing ones in ECA Names
for (i in 1:nrow(merged_weights)){
  
  if (str_count(merged_weights$Site[i], "[0-9]") <= 2){
    
    merged_weights$Site[i] = paste0("0", merged_weights$Site[i])
    
  }
  
  else {
    
    merged_weights$Site[i] = merged_weights$Site[i]
  }
  
}

## Clean Data 
  ## EC: 
  
  ## EV_003, EV_004 - not publishing this data
  ## EV_020_MOI-3: Spilled, needs to be removed

eca_grav <- merged_weights %>% 
  unite("Sample_Name", c("EC", "Site", "Rep"), sep = "_") %>% 
  #filter(!grepl("EV_003|EV_004|EV_020_MOI-3", Sample_Name)) %>% 
  separate(Sample_Name, c("Sample_Site", "Rep"), sep = "-", remove = FALSE) %>% 
  separate(Sample_Site, c("EC", "Site", "MOI"), sep = "_") %>% 
  mutate(grav_eca = ((wet_weight_g - true_dry_weight_g)/true_dry_weight_g)*100) %>% 
  group_by(Site) %>% 
  mutate(average_eca = mean(grav_eca)) %>% 
  mutate(cv_eca = (sd(grav_eca)/mean(grav_eca))*100) %>%   dplyr::select(c("Sample_Name", "Site", "Rep", "grav_eca", "average_eca", "cv_eca"))

## Remove Outliers if Needed ####

## EC: 

## EC_057_MOI-3 spilled before weighing
## EC_084_MOI-2 typo? negative value

## EV: highest CV is 45%, but no distinct outliers

merged_weights_average = merged_weights %>%
  #filter(Sample_Name != "EC_057_MOI-3") %>% 
  #filter(Sample_Name != "EC_084_MOI-2") %>% 
  separate(Sample_Name, c("EC", "Site", "Rep"), sep = "_") %>% 
  group_by(Site) %>% 
  mutate(average = mean(percent_water_content_dry)) %>% 
  mutate(cv = sd(percent_water_content_dry)/mean(percent_water_content_dry))


eca_cv <- eca_grav %>% 
  distinct(Site, .keep_all = TRUE)

ggplot(eca_cv, aes(x = cv_eca)) +
  geom_histogram()


clean_eca <- merged_weights %>% 
  filter(!grepl("EV_003|EV_004|EV_020_MOI-3", Sample_Name)) %>% 
  mutate(Gravimetric_Moisture = ((wet_weight_g - true_dry_weight_g)/true_dry_weight_g)*100) %>% 
  dplyr::select(c(Sample_Name, wet_weight_g, true_dry_weight_g, percent_water_content_dry, Gravimetric_Moisture, flag)) %>% 
  rename(Wet_Sediment_Mass_g = wet_weight_g) %>% 
  rename(Dry_Sediment_Mass_g = true_dry_weight_g) %>% 
  mutate(Water_Mass_g = Wet_Sediment_Mass_g - Dry_Sediment_Mass_g) %>% 
  mutate(Gravimetric_Moisture = round(Gravimetric_Moisture, 2)) %>% 
  mutate(Methods_Deviation = "N/A") %>% 
 # mutate(Sample_Name = str_replace(Sample_Name, "EC", "CM")) %>% 
  dplyr::select(c(Sample_Name, Wet_Sediment_Mass_g, Dry_Sediment_Mass_g, Water_Mass_g, Gravimetric_Moisture, flag, Methods_Deviation))

write.csv(clean_eca, paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/MOI/03_ProcessedData/",study.code,"_MOI_ReadyForBoye_", Sys.Date(),".csv"), row.names = F)

### Run ECA data through ECA Moisture Outliers Code to compare to ICON dry weights and remove outliers
