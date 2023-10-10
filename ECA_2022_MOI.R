##### Load Library #####
library(tidyverse)
library(dplyr)
library(readxl)

##### Load data ######
rm(list=ls());graphics.off()

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
wet_weight <- moisture_clean %>% dplyr::select(sample_name, date, tare_weight_g, sample_weight_g) %>% 
  filter(!grepl("-9999" , moisture_clean$tare_weight_g))

colnames(wet_weight) [4] <- "wet_weight_g"

# Subseting data by dry weight
dry_weight <- moisture_clean %>% dplyr::select(sample_name, date, tare_weight_g, sample_weight_g) %>% 
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

merged_weights_average = merged_weights %>%
  filter(Sample_Name != "EC_57_MOI-3") %>% 
  filter(Sample_Name != "EC_84_MOI-2") %>% 
  separate(Sample_Name, c("EC", "Site", "Rep"), sep = "_") %>% 
  group_by(Site) %>% 
  mutate(average = mean(percent_water_content_dry)) %>% 
  mutate(cv = sd(percent_water_content_dry)/mean(percent_water_content_dry))

# Export data
write.csv(merged_weights, file.path(data.out.path,"EC_Moisture_Content_2022.csv"), quote = F, row.names = F) #removes "" from data and no need for row names here, puts csv into the data file

## Comparisons to ICON DWT

icon <- read_csv(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ICON_ModEx_SSS/Subsampling/Sediment_Mass_Water_Calculations/ICON_Sediment_Incubations_Sediment_g_and_Water_mL_merged_by_laan208_on_2023-10-02.csv"))

icon_grav <- icon %>% 
  mutate(grav_moi = (((Wet_Sediment_Mass_g - Dry_Sediment_Mass_g)/Dry_Sediment_Mass_g)*100)) %>% 
  relocate(grav_moi, .after = Sample_ID) %>% 
  separate(Sample_ID, c("EC", "Site", "Rep"), sep = "_") %>% 
  filter(!grepl("-4",Rep)) %>% 
  filter(!grepl("-5",Rep)) %>% 
  filter(!grepl("-6",Rep)) %>% 
  separate(Rep, c("INC", "Rep"), sep = "-") %>% 
  group_by(Site) %>% 
  mutate(average_icon = mean(grav_moi)) %>% 
  mutate(cv_icon = sd(grav_moi)/mean(grav_moi)) %>% 
  relocate(average_icon, .after = grav_moi) %>% 
  relocate(cv_icon, .after = average_icon) %>% 
  dplyr::select(-c(EC, INC, Methods_Deviation))


eca_grav <- merged_weights_average %>% 
  separate(Rep, c("INC", "Rep"), sep = "-") %>% 
  dplyr::select(-c(EC, INC))

icon_eca <- merge(eca_grav, icon_grav, by = c("Site", "Rep"), all = TRUE)

#028 not used for ECA

icon_eca <- icon_eca %>% 
  drop_na(average) %>% 
  relocate(percent_water_content_dry, .after = Site) %>% 
  relocate(average, .after = percent_water_content_dry) %>% 
  relocate(cv, .after = average) %>% 
  relocate(grav_moi, .after = cv) %>% 
  relocate(average_icon, .after = grav_moi) %>% 
  relocate(cv_icon, .after = average_icon) %>% 
  mutate(grav_moi_diff = (average - average_icon)) %>% 
  relocate(grav_moi_diff, .after = Site) %>% 
  filter(Site != "028")

high_cv <- icon_eca %>% 
  filter(cv >= 0.27) %>% 
  select(c(Site, percent_water_content_dry, grav_moi)) %>% 
  pivot_longer(!Site, names_to = "Project", values_to = "Grav_Moi")

ggplot(high_cv)+
  geom_boxplot(aes(x = Site, y = Grav_Moi, fill = Project))+
  scale_fill_discrete(labels = c("ICON", "ECA"))
 
high_cv_outliers <- icon_eca %>% 
  filter(cv >= 0.27) %>% 
  filter(!(Site == "082" & percent_water_content_dry > 100)) %>% 
  filter(!(Site == "106" & percent_water_content_dry > 100)) %>% 
filter(!(Site == "060" & percent_water_content_dry > 50)) %>% 
  select(c(Site, percent_water_content_dry, grav_moi)) %>% 
  pivot_longer(!Site, names_to = "Project", values_to = "Grav_Moi")

ggplot(high_cv_outliers)+
  geom_boxplot(aes(x = Site, y = Grav_Moi, fill = Project))+
  scale_fill_discrete(labels = c("ICON", "ECA"))


avg_icon_eca <- icon_eca %>% 
  distinct(Site, .keep_all = TRUE)

ggplot(avg_icon_eca, aes(x = Site))+
  geom_point(aes(y = average, color = "ECA"))+
  geom_point(aes(y = average_icon, color = "ICON"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(avg_icon_eca, aes(x = Site))+
  geom_point(aes(y = grav_moi_diff))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  


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