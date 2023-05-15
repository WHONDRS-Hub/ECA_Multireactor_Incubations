##### Load Library #####
library(tidyverse)
library(dplyr)
library(readxl)

pnnl.user = 'laan208'

# Set working directory to data file
  #For Macbook Users
data.in.path <- file.path("/Users/delg580/Library/CloudStorage/OneDrive-SharedLibraries-PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/MOI/01_RawData")

data.out.path <- file.path("/Users/delg580/Library/CloudStorage/OneDrive-SharedLibraries-PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/MOI/03_ProcessedData")

  #For Windows Users

data.out.path <- file.path("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/MOI/03_ProcessedData")

# Load Data
moisture <- read_excel(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/MOI/01_RawData/2022_Data_Raw_MOI_ECA_EC.xlsx")) %>% janitor::clean_names()

#Remove data from Jars experiment
moisture_clean <- moisture %>% filter(!grepl("Jars", moisture$from_jars_or_40_m_l_vials))

# Subseting data by wet weight
wet_weight <- moisture_clean %>% select(sample_name, date, tare_weight_g, sample_weight_g) %>% 
  filter(!grepl("-9999" , moisture_clean$tare_weight_g))

colnames(wet_weight) [4] <- "wet_weight_g"

# Subseting data by dry weight
dry_weight <- moisture_clean %>% select(sample_name, date, tare_weight_g, sample_weight_g) %>% 
  filter(grepl("-9999", moisture_clean$tare_weight_g))

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

# Calculating percent water content on a dry basis
merged_weights$percent_water_content_dry <- ((merged_weights$wet_weight_g - merged_weights$true_dry_weight_g)/merged_weights$true_dry_weight_g)*100

# Calculating percent water content on a wet basis
merged_weights$percent_water_content_wet <- ((merged_weights$wet_weight_g - merged_weights$true_dry_weight_g)/merged_weights$wet_weight_g)*100


# Clean up data
cols_to_drop <- c ("date_x","date_y", "tare_weight_g_y","tare_weight_g_x", "dry_weight_g" )
merged_weights <- merged_weights[ , !(names(merged_weights) %in% cols_to_drop)]

range(merged_weights$percent_water_content_wet) 
range(merged_weights$percent_water_content_dry)

merged_weights$flag[is.na(merged_weights$wet_weight_g)] = "Missing Wet Weight"
  

# Export data
write.csv(merged_weights, file.path(data.out.path,"EC_Moisture_Content_2022.csv"), quote = F, row.names = F) #removes "" from data and no need for row names here, puts csv into the data file

#Ploting wet vs dry percent water content
pdf(file.path(data.out.path,"wet_vs_dry_percent_water_content.pdf"))  
plot(merged_weights$percent_water_content_wet ~ merged_weights$percent_water_content_dry)
dev.off()

pdf(file.path(data.out.path,"percent_water_content_dry.pdf"))
hist(merged_weights$percent_water_content_dry)
dev.off()


pdf(file.path(data.out.path,"percent_water_content_wet.pdf")) 
hist(merged_weights$percent_water_content_wet)
dev.off()






