## This script is for getting dry sediment mass, wet sediment mass, water mass, and gravimetric water over the 21 day ECA incubation

## This file currently using published moisture tin data - will need to be updated if calculating these values before publishing

## Read in Libraries ####
library(readxl);library(tidyverse)

rm(list=ls());graphics.off()

## Inputs for Data ####
pnnl.user = 'laan208'
study.code = 'EL' #EC, EV, MEL
date = '2024-10-24' #Date of Moisture Tin File, use for analyzing data before publishing
year = '2024' # Year of ECA incubation
  

## Pull in Data ####

## Moisture tin data already summarized
if (study.code == "EC"){
dry_wt <- read.csv("C:/GitHub/ECA_Multireactor_Incubations/Data/v3_CM_SSS_Sediment_Sample_Data_Summary.csv", skip = 2) %>% 
  slice(-1:-11, -495) %>% 
  select(c(Sample_Name, Mean_62948_Gravimetric_Moisture_g_per_g)) 
} else if (study.code == "EV"){
  
  dry_wt <- read.csv(paste0("C:/Users/",pnnl.user,"/OneDrive - PNNL/Shared Documents - Core Richland and Sequim Lab-Field Team/Data Generation and Files/ECA/MOI/03_ProcessedData/",study.code,"_MOI_ReadyForBoye_", date,".csv"))
    
} else if (study.code == "EL"){
  
  dry_wt <- read.csv(paste0("C:/Users/",pnnl.user,"/OneDrive - PNNL/Shared Documents - Core Richland and Sequim Lab-Field Team/Data Generation and Files/ECA/MOI/03_ProcessedData/M",study.code,"_MOI_ReadyForBoye_", date,".csv")) %>% 
    filter(!grepl("removed outlier", flag)) %>% 
    separate(Sample_Name, c("MEL", "Site", "MOI"), sep = "_") %>% 
    separate(MOI, c("MOI", "Rep"), sep = "-") %>% 
    group_by(Site, MOI) %>% 
    mutate(Mean_Gravimetric_Moisture = mean(Gravimetric_Moisture))
}

## 21 day incubation data 
moisture <- read_xlsx(paste0("C:/Users/",pnnl.user,"/OneDrive - PNNL/Shared Documents - Core Richland and Sequim Lab-Field Team/Data Generation and Files/ECA/INC/01_RawData/",year,"_Data_Raw_INC_ECA_",study.code,".xlsx"))

## Clean Data ####

## Moisture tins: 

# Rearrange gravimetric moisture to calculate average gram of wet sediment (water mass + 1 g dry sediment mass) divided by average gram of dry sediment per Site. This makes it possible to calculate dry sediment mass and water mass from wet sediment mass

if (study.code == "EC") {

dry_wt_averages <- dry_wt %>% 
  filter(grepl("CM", Sample_Name)) %>% 
  filter(Mean_62948_Gravimetric_Moisture_g_per_g != -9999) %>% 
  mutate(average_wet_g_by_dry_g = 1 + as.numeric((Mean_62948_Gravimetric_Moisture_g_per_g))) %>% 
  mutate(Sample_ID = str_replace(Sample_Name, "_Sediment", "")) %>% 
  mutate(Sample_ID = str_replace(Sample_ID, "CM", "EC")) %>% 
  select(-c(Sample_Name))

} else if (study.code == "EL") {
  
  dry_wt_averages <- dry_wt %>% 
    filter(Gravimetric_Moisture != -9999) %>% 
    mutate(average_wet_g_by_dry_g = 1 + as.numeric((Mean_Gravimetric_Moisture))) %>% 
    mutate(MOI = str_replace(MOI, "MS", "IN")) %>% 
    mutate(MEL = str_replace(MEL, "MEL", "EL")) %>% 
    unite(Sample_ID, c(MEL, Site, MOI), sep = "_") %>% 
    distinct(Sample_ID, .keep_all = T) %>% 
    select(c(Sample_ID,average_wet_g_by_dry_g, Mean_Gravimetric_Moisture))
  
}
  

## 21 Day Incubation

# For ECA: 37,38 was incubated on 11/23 so needs to be removed from sheet, 72-5 were not incubated (not enough sediment), 12-D5 didn't have 20 g of sediment, 21 and 33 incubated on 9/28. Remove Jar samples

all_moisture = moisture %>% 
  mutate(Date = as.Date(Date)) %>% 
  mutate(Sample_weight_Fill_g = as.numeric(Sample_weight_Fill_g)) %>% 
  select(c(Sample_Name, Date, Tare_weight_g, Sample_weight_g, Sample_weight_Fill_g)) %>% 
  filter(!grepl("EC_01|EC_02|EC_03|EC_06|EC_07|EC_04|EC_08|EC_10|EC_15", Sample_Name)) %>% # these are jar samples
 filter(!(Date == "2022-11-23" & grepl("EC_37", Sample_Name))) %>%
 filter(!(Date == "2022-11-23" & grepl("EC_38", Sample_Name))) %>%
 filter(!(Date == "2022-09-28" & grepl("EC_21", Sample_Name))) %>%
 filter(!(Date == "2022-09-28" & grepl("EC_33", Sample_Name))) %>%
 filter(Sample_Name != "EC_72_INC-D5") %>% # don't have these replicates
 filter(Sample_Name != "EC_72_INC-W5")  %>%
  separate(Sample_Name, into = c("EC", "Site", "INC"), sep = "_", remove = FALSE)

# Add "0" to start of sample kit names that don't have it for EC samples
if (study.code == "EC") {

for (i in 1:nrow(all_moisture)){
  
  if (str_count(all_moisture$Site[i], "[0-9]") <= 2){
    
    all_moisture$Site[i] = paste0("0", all_moisture$Site[i])
    
  }
  
  else {
    
    all_moisture$Site[i] = all_moisture$Site[i]
  }
  
}}

if (study.code == "EC"){

merged <- all_moisture %>% 
  unite(Sample_ID, c("EC", "Site")) %>% 
  unite(Sample_Name, c("Sample_ID", "INC"), remove = FALSE) %>% 
  left_join(dry_wt_averages, by = "Sample_ID") 

} else if (study.code == "EL") {
    
  merged <- all_moisture %>% 
    separate(INC, c("INC", "Rep")) %>% 
    unite(Sample_ID, c("EC", "Site", "INC"), remove = F) %>% 
    left_join(dry_wt_averages, by = "Sample_ID") 
    
  }

## All Samples - wet and dry (all dates) ####

## Calculate wet sediment mass, dry sediment mass, added water, and total water for each date of 21 day incubation

location = c("-W1", "-W2", "-W3", "-W4", "-W5", "-D1", "-D2", "-D3", "-D4", "-D5")

wet_dry <- as.data.frame(matrix(NA, ncol = 9, nrow =1))

colnames(wet_dry) = c("Sample_Name","Date",  "Tare_weight_g", "Sample_weight_g","Sample_weight_Fill_g", "Dry_weight_sed_g", "Water_added_g", "Water_total_g", "Gravimetric_Moisture")

for (i in 1:length(location)){
  
  all_dates= as.data.frame(matrix(NA, ncol =9, nrow = length(unique(merged$Sample_Name))))
  
  colnames(all_dates) = c("Sample_Name","Date",  "Tare_weight_g", "Sample_weight_g","Sample_weight_Fill_g", "Dry_weight_sed_g","Water_added_g", "Water_total_g", "Gravimetric_Moisture")
  
  #Subset by replicate
  data_location_subset = merged[grep(location[i],merged$Sample_Name),]
  
  unique.incubations = unique(data_location_subset$Sample_Name)
  
  
  for (j in 1:length(unique.incubations)){
    
    #Subset by site and replicate (e.g. all dates of added water for EV_001_INC-W1)
    data_site_subset = subset(data_location_subset, data_location_subset$Sample_Name == unique.incubations[j])
    
    #put in right date order
    data_site_subset <- data_site_subset[with(data_site_subset, order(Date)),]
    
    #Fill Tare Weight down column and get dry sediment mass in each vial
    data_site_subset = data_site_subset %>%
      mutate(Tare_weight_g = first(Tare_weight_g)) %>% 
      mutate(Dry_weight_sed_g = (first(Sample_weight_g) - Tare_weight_g)*(1/average_wet_g_by_dry_g)) 
  
    merge_dates= as.data.frame(matrix(NA, ncol = 8, nrow = length(unique(all_moisture$Sample_Name))))
    
    colnames(merge_dates) = c("Sample_Name","Date",  "Tare_weight_g", "Sample_weight_g","Sample_weight_Fill_g", "Dry_weight_sed_g", "Water_added_g", "Water_total_g")
    
  for (k in 1:nrow(data_site_subset)){
    
    # Separate Wet Treatments
    if (data_site_subset$Sample_weight_Fill_g[k] > -9999) {
     
      merge_dates$Sample_Name[k] = as.character(data_site_subset$Sample_Name[k])
    
    merge_dates$Date[k] = as.character(data_site_subset$Date[k]) 
    
    merge_dates$Tare_weight_g[k] = as.numeric(data_site_subset$Tare_weight_g[k]) 
    
    # Get actual sample weight before adding water without tare weight
    merge_dates$Sample_weight_g[k] = as.numeric(data_site_subset$Sample_weight_g[k] - data_site_subset$Tare_weight_g[k])
    
    # Get actual sample weight after adding water without tare weight
    merge_dates$Sample_weight_Fill_g[k] = as.numeric(data_site_subset$Sample_weight_Fill_g[k] - data_site_subset$Tare_weight_g[k])
    
    merge_dates$Dry_weight_sed_g[k] = as.numeric(data_site_subset$Dry_weight_sed_g[k])
    
    # Get mass of water added per day
    merge_dates$Water_added_g[k] = as.numeric(data_site_subset$Sample_weight_Fill_g[k] - data_site_subset$Sample_weight_g[k])
    
    # Get total water mass
  merge_dates$Water_total_g[k] = merge_dates$Sample_weight_Fill_g[k] - merge_dates$Dry_weight_sed_g[k]
  
  # calculate gravimetric moisture
  merge_dates$Gravimetric_Moisture[k] = merge_dates$Water_total_g[k]/merge_dates$Dry_weight_sed_g[k]
  

    }
    
    # Separate Dry Treatments
    else {
      
      merge_dates$Sample_Name[k] = as.character(data_site_subset$Sample_Name[k])
      
      merge_dates$Date[k] = as.character(data_site_subset$Date[k]) 
      
      merge_dates$Tare_weight_g[k] = as.numeric(data_site_subset$Tare_weight_g[k]) 
      
      # Get actual sample weight without tare weight
      merge_dates$Sample_weight_g[k] = as.numeric(data_site_subset$Sample_weight_g[k] - data_site_subset$Tare_weight_g[k])
      
      # No water added to dry treatments
      merge_dates$Water_added_g[k] = as.numeric(0)
      
      merge_dates$Dry_weight_sed_g[k] = as.numeric(data_site_subset$Dry_weight_sed_g[k])
      
      merge_dates$Av_dry_weight_sed_g[k] = as.numeric(data_site_subset$Av_dry_weight_sed_g[k])
      
      # Get total water mass in sample 
      merge_dates$Water_total_g[k] = merge_dates$Sample_weight_g[k] - merge_dates$Dry_weight_sed_g[k]
      
      merge_dates$Av_water_total_g[k] = merge_dates$Sample_weight_Fill_g[k] - merge_dates$Av_dry_weight_sed_g[k]
      
      merge_dates$Gravimetric_Moisture[k] = merge_dates$Water_total_g[k]/merge_dates$Dry_weight_sed_g[k]
      
    }
  
    }
    
    #all_dates = all_dates %>%
     # mutate(Sample_weight_initial_g = first(Sample_weight_initial_g), round, 4)
    
   all_dates = rbind(all_dates, merge_dates)
    
  }
  
  
  
  all_dates = all_dates[!is.na(all_dates$Sample_Name),]
  
  wet_dry = rbind(wet_dry, all_dates)
  
}

wet_dry = wet_dry[-1,] 

## Column Descriptions:

# Wet_Sediment_Mass_g: mass of wet sediment intially added to vial on first day of incubations for W and D treatments. For subsequent days of wet treatment: mass of original wet sed + water before adding more to bring to constant weight. For dry treatments: mass of wet sediment in vial

#Wet_Sediment_Mass_Added_Water_g: mass of wet sediment + water for wet treatments after adding more to bring to constant weight. 

#Dry_Sediment_Mass_g: dry sediment mass calculated from average moisture content of sediment tins (being published with ICON) and first day Wet_Sediment_Mass_g

#Added_Water_Mass_g: For wet treatments, on the first day, amount of water added to get to 1 cm headspace. Subsequent days are amount of water added to maintain water losses

#Total_Water_Mass_g: total added water + water from wet sediment mass 

#Gravimetric Moisture: Total_Water_Mass_g / Dry_Sediment_Mass_g

final_wet_dry <- wet_dry %>% 
  dplyr::select(-c(Tare_weight_g)) %>% 
  rename(Wet_Sediment_Mass_g = Sample_weight_g) %>% 
  rename(Wet_Sediment_Mass_Added_Water_g = Sample_weight_Fill_g) %>% 
  rename(Dry_Sediment_Mass_g = Dry_weight_sed_g) %>% 
  rename(Total_Water_Mass_g = Water_total_g) %>% 
  rename(Added_Water_Mass_g = Water_added_g) %>% 
  distinct(Sample_Name, Date, .keep_all = TRUE) %>% 
  mutate(across(c("Dry_Sediment_Mass_g": "Gravimetric_Moisture"),round,2)) %>% 
  mutate(Wet_Sediment_Mass_Added_Water_g= if_else(is.na(Wet_Sediment_Mass_Added_Water_g), -9999, Wet_Sediment_Mass_Added_Water_g)) %>% 
  mutate(Added_Water_Mass_g= if_else(Added_Water_Mass_g == 0, -9999, Added_Water_Mass_g)) %>% 
  mutate(Dry_Sediment_Mass_g = if_else(is.na(Dry_Sediment_Mass_g), -9999, Dry_Sediment_Mass_g)) %>% 
  mutate(Total_Water_Mass_g = if_else(is.na(Total_Water_Mass_g), -9999, Total_Water_Mass_g)) 



## Export File (all dates)

write.csv(final_wet_dry,paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/INC/03_ProcessedData/",study.code,"_Drying_Masses_ReadyForBoye_",Sys.Date(),".csv"), row.names = F)  

ggplot(final_wet_dry, aes(x = Dry_Sediment_Mass_g))+
  geom_histogram()


## Start Summary File ####

## Start loop (similar to above) to only pull out Initial water mass/gravimetric moisture (first day of incubation) and final water mass/gravimetric moisture (last day of incubation)

location = c("-W1", "-W2", "-W3", "-W4", "-W5", "-D1", "-D2", "-D3", "-D4", "-D5")

water_mass <- as.data.frame(matrix(NA, ncol = 6, nrow =1))

colnames(water_mass) = c("Sample_Name", "Initial_Water_mass_g", "Final_Water_mass_g", "Dry_Sediment_Mass_g", "Initial_Gravimetric_Moisture", "Final_Gravimetric_Moisture")

for (i in 1:length(location)){
  
  summary_file = as.data.frame(matrix(NA, ncol = 6, nrow = length(unique(final_wet_dry$Sample_Name))))
  
  colnames(summary_file) = c("Sample_Name", "Initial_Water_mass_g", "Final_Water_mass_g", "Dry_Sediment_Mass_g", "Initial_Gravimetric_Moisture", "Final_Gravimetric_Moisture")
  
  # Subset by replicate
  summary_location_subset = final_wet_dry[grep(location[i],final_wet_dry$Sample_Name),]
  
  # Note unique samples
  unique.incubations = unique(summary_location_subset$Sample_Name)
  
  for (j in 1:length(unique.incubations)){
    
    #subset by site
    summary_site_subset = subset(summary_location_subset, summary_location_subset$Sample_Name == unique.incubations[j])
    
    #put in chronological order 
    summary_site_subset <- summary_site_subset[with(summary_site_subset, order(Date)),]
    
    summary_site_subset <- summary_site_subset %>% 
      mutate(Initial_Water_mass_g = first(Total_Water_Mass_g)) %>% # pull out initial water mass from first day of incubations
      mutate(Final_Water_mass_g = last(Total_Water_Mass_g)) %>% # pull out final water mass from last day of incubations
      mutate(Initial_Gravimetric_Moisture = first(Gravimetric_Moisture)) %>% # pull out initial gravimetric moisture from first day of incubations
      mutate(Final_Gravimetric_Moisture = last(Gravimetric_Moisture)) %>%  # pull out final gravimetric moisture from last day of incubations
      select(c(Sample_Name, Initial_Water_mass_g, Final_Water_mass_g, Dry_Sediment_Mass_g, Initial_Gravimetric_Moisture, Final_Gravimetric_Moisture)) %>% 
      filter(row_number()==1)
    
    summary_file <- rbind(summary_file, summary_site_subset)
    
  }
  
  water_mass = rbind(water_mass, summary_file)
  
}  
    
water_mass <- water_mass %>% 
  drop_na(Sample_Name)

# Export file

write.csv(water_mass,paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/INC/03_ProcessedData/",study.code,"_Drying_Masses_Summary_ReadyForBoye_on_",Sys.Date(),".csv"), row.names = F)
