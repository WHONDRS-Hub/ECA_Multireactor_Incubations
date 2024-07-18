##### Load Library #####
library(tidyverse)
library(readxl)

##### Load data ######
rm(list=ls());graphics.off()

pnnl.user = 'laan208'
study.code = 'EC'
year = '2022'

# Set working directory to data file

#For Windows Users
data.out.path <- file.path("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/MOI/03_ProcessedData")

#For Mac Users 
#data.out.path <- file.path("/Users/",pnnl.user,"/Library/CloudStorage/OneDrive-SharedLibraries-PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/MOI/03_ProcessedData")

# Load Data
moisture <- read_excel(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/MOI/01_RawData/",year,"_Data_Raw_MOI_ECA_",study.code,".xlsx")) %>% janitor::clean_names()

#For Macbook Users
#moisture <- read_excel(paste0("/Users/",pnnl.user,"/Library/CloudStorage/OneDrive-SharedLibraries-PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/MOI/01_RawData/",year,"_Data_Raw_MOI_ECA_",study.code,".xlsx")) %>% janitor::clean_names()

# Subsetting data by wet weight
wet_weight <- moisture %>% 
  dplyr::select(sample_name, date, tare_weight_g, sample_weight_g_tared) %>% 
  filter(!grepl("-9999" , tare_weight_g))

colnames(wet_weight) [4] <- "wet_weight_g"

# Subsetting data by dry weight
dry_weight <- moisture %>% 
  dplyr::select(sample_name, date, tare_weight_g, sample_weight_g_tared) %>% 
  filter(grepl("-9999", tare_weight_g))

colnames(dry_weight) [4] <- "dry_weight_g"


dim(wet_weight); length(unique(wet_weight$sample_name)) #number of unique entries
dim(dry_weight); length(unique(dry_weight$sample_name)) #number of unique entries
length(which(wet_weight$sample_name %in% dry_weight$sample_name)) #which ID's in wet weight file are found in the dry weight file (missing EC_23 missing wet weight)

# Merged wet and dry weights by sample name
merged_weights <- merge(x = wet_weight, y = dry_weight, by = "sample_name", all = TRUE) %>% janitor::clean_names()
dim(merged_weights)
head(merged_weights)
names(merged_weights)

# Calculate True dry weight ((dry weight + moisture tin) - moisture tin)
merged_weights$true_dry_weight_g <- merged_weights$dry_weight_g - merged_weights$tare_weight_g_x 


# Calculate Gravimetric Water Content 
#(water mass/dry sediment mass)*100
merged_weights$percent_water_content_dry <- ((merged_weights$wet_weight_g - merged_weights$true_dry_weight_g)/merged_weights$true_dry_weight_g)*100


# Clean up data
cols_to_drop <- c ("date_x","date_y", "tare_weight_g_y","tare_weight_g_x", "dry_weight_g" )
merged_weights <- merged_weights[ , !(names(merged_weights) %in% cols_to_drop)]

# Add extra 0 to missing ones in ECA Names
merged_weights <- merged_weights %>% 
  separate(sample_name, c("EC", "Site", "Rep"), sep = "_")

for (i in 1:nrow(merged_weights)){
  
  if (str_count(merged_weights$Site[i], "[0-9]") <= 2){
    
    merged_weights$Site[i] = paste0("0", merged_weights$Site[i])
    
  }
  
  else {
    
    merged_weights$Site[i] = merged_weights$Site[i]
  }
  
}

merged_weights <- merged_weights %>% 
  unite("Sample_Name", c("EC", "Site", "Rep"), sep = "_")

## Get average moisture content and stdev of moisture content from replicates 
  
eca_grav <- merged_weights %>% 
  separate(Sample_Name, c("Sample_Site", "Rep"), sep = "-", remove = FALSE) %>% 
  separate(Sample_Site, c("EC", "Site", "MOI"), sep = "_") %>% 
  group_by(Site) %>% 
  mutate(average_eca = mean(percent_water_content_dry)) %>% 
  mutate(cv_eca = (sd(percent_water_content_dry)/mean(percent_water_content_dry))*100) %>%  
  dplyr::select(c("Sample_Name", "Site", "Rep", "wet_weight_g", "true_dry_weight_g", "percent_water_content_dry", "average_eca", "cv_eca")) %>% 
  ungroup()

## Remove Obivous Outliers (those with notes) if Needed ####

## EC: 

# EC_057_MOI-3 spilled before weighing
# EC_084_MOI-2 typo? negative value

## EV: highest CV is 45%, but no distinct outliers
# EV_003, EV_004 - not publishing this data, sampled twice 
# EV_020_MOI-3: Spilled, needs to be removed

eca_grav_rem = eca_grav %>%
  filter(Sample_Name != "EC_057_MOI-3") %>% 
  filter(Sample_Name != "EC_084_MOI-2") %>% 
  filter(!grepl("EV_003|EV_004|EV_020_MOI-3", Sample_Name)) %>% # use this when processing EV data
  group_by(Site) %>% 
  mutate(average_rem = mean(percent_water_content_dry)) %>% 
  mutate(cv_rem = sd(percent_water_content_dry)/mean(percent_water_content_dry)*100) %>% 
  ungroup() %>% 
  right_join(eca_grav, by = "Sample_Name") %>%
  mutate(flag = if_else(is.na(wet_weight_g.y),  "Missing Wet Weight", if_else(is.na(average_rem), "removed outlier", "NA"))) # Add flag to sites missing wet weight then flag removed samples

# Keep average site moisture and look at histogram
eca_cv <- eca_grav %>% 
  distinct(Site, .keep_all = TRUE)

ggplot(eca_cv, aes(x = cv_eca)) +
  geom_histogram()

# Make final DF for publishing/running through outlier script
clean_eca <- eca_grav_rem %>% 
 rename(Gravimetric_Moisture = percent_water_content_dry.y) %>% 
  rename(Wet_Sediment_Mass_g = wet_weight_g.y) %>%
  rename(Dry_Sediment_Mass_g = true_dry_weight_g.y) %>% 
  mutate(Water_Mass_g = Wet_Sediment_Mass_g - Dry_Sediment_Mass_g) %>% 
  mutate(Gravimetric_Moisture = round(Gravimetric_Moisture, 2)) %>% 
  mutate(Methods_Deviation = "N/A") %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "EC", "CM")) %>% #change EC to CM for data publishing
  dplyr::select(c(Sample_Name, Gravimetric_Moisture, Wet_Sediment_Mass_g, Dry_Sediment_Mass_g, Water_Mass_g, flag, Methods_Deviation))

write.csv(clean_eca, paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/MOI/03_ProcessedData/",study.code,"_MOI_ReadyForBoye_", Sys.Date(),".csv"), row.names = F)

### Run ECA data through ECA Moisture Outliers Code to compare to ICON dry weights and remove outliers
