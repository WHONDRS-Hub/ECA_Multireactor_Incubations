library(readxl);library(corrplot);library(tidyverse)

rm(list=ls());graphics.off()

pnnl.user = 'laan208'

## ECA Dry Weights ####
dry_wt <- read.csv(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/MOI/03_ProcessedData/EV_MOI_ReadyForBoye_04-18-2024.csv"))

dry_wt_averages <- dry_wt %>% 
 # filter(is.na(flag)) %>% 
  separate(Sample_Name, c("EC", "Site", "Rep"), sep = "_", remove = FALSE) %>% 
  unite(Sample_ID, c("EC", "Site")) %>% 
  group_by(Sample_ID) %>% 
  mutate(average_grav = mean(Gravimetric_Moisture)) %>% 
  mutate(cv = (sd(Gravimetric_Moisture)/mean(Gravimetric_Moisture))*100)%>% 
  distinct(Sample_ID, .keep_all = TRUE) %>% 
  select(c(Sample_ID, average_grav)) %>% mutate(average_wet_g_by_dry_g = 1 + (average_grav/100)) #%>% 
 # mutate(Sample_ID = str_replace(Sample_ID, "CM", "EC"))
  
moisture <- paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/INC/01_RawData/2023_Data_Raw_INC_EV.xlsx")
#2022_Data_Raw_INC_ECA_EC.xlsx
#2023_Data_Raw_INC_EV

all_moisture <- read_xlsx(moisture)

all_moisture$Date <- as.Date(all_moisture$Date)

all_moisture$Sample_weight_Fill_g <- as.numeric(all_moisture$Sample_weight_Fill_g)

# For ECA: 37,38 was incubated on 11/23 so needs to be removed from sheet, 72-5 were not incubated (not enough sediment), 12-D5 didn't have 20 g of sediment, 21 and 33 incubated on 9/28

all_moisture <- all_moisture %>% 
  #filter(`Jars or 40 mL vials` != "Jars") %>% 
  #filter(!(Date == "2022-11-23" & grepl("EC_37", Sample_Name))) %>% 
 # filter(!(Date == "2022-11-23" & grepl("EC_38", Sample_Name))) %>% 
  #filter(!(Date == "2022-09-28" & grepl("EC_21", Sample_Name))) %>% 
 # filter(!(Date == "2022-09-28" & grepl("EC_33", Sample_Name))) %>% 
 # filter(Sample_Name != "EC_72_INC-D5") %>% 
 # filter(Sample_Name != "EC_72_INC-W5")  %>% 
  separate(Sample_Name, into = c("EC", "Site", "INC"), sep = "_", remove = FALSE)

#add "0" to start of sample kit names that don't have it

for (i in 1:nrow(all_moisture)){
  
  if (str_count(all_moisture$Site[i], "[0-9]") <= 2){
    
    all_moisture$Site[i] = paste0("0", all_moisture$Site[i])
    
  }
  
  else {
    
    all_moisture$Site[i] = all_moisture$Site[i]
  }
  
}

all_moisture <- all_moisture %>% 
  unite(Sample_ID, c("EC", "Site")) %>% 
  unite(Sample_Name, c("Sample_ID", "INC"), sep = "_", remove = FALSE)


merged <- left_join(all_moisture, dry_wt_averages, by = "Sample_ID") %>% 
  dplyr::select(-c(Notes, ...7, ...8, ...9, ...10, ...11, ...12, ...13, ...14, ...15, ...16)) %>% 
  distinct(Sample_Name, Date, .keep_all = TRUE)
  # dplyr::select(-c(Notes, `Jars or 40 mL vials`, ...8, ...9, ...10, `Study Code`, percent_water_content_wet, percent_water_content_dry, true_dry_weight_g, wet_weight_g))


##All Samples - wet and dry (all dates) ####

location = c("-W1", "-W2", "-W3", "-W4", "-W5", "-D1", "-D2", "-D3", "-D4", "-D5")

wet_dry <- as.data.frame(matrix(NA, ncol = 8, nrow =1))

colnames(wet_dry) = c("Sample_Name","Date",  "Tare_weight_g", "Sample_weight_g","Sample_weight_Fill_g", "Dry_weight_sed_g", "Water_added_g", "Water_total_g")

all_dates = as.data.frame(matrix(NA, ncol = 8, nrow = length(unique(all_moisture$Sample_Name))))

colnames(all_dates) = c("Sample_Name","Date",  "Tare_weight_g", "Sample_weight_g","Sample_weight_Fill_g", "Dry_weight_sed_g","Water_added_g", "Water_total_g")

for (i in 1:length(location)){
  
  all_dates= as.data.frame(matrix(NA, ncol =8, nrow = length(unique(merged$Sample_Name))))
  
  colnames(all_dates) = c("Sample_Name","Date",  "Tare_weight_g", "Sample_weight_g","Sample_weight_Fill_g", "Dry_weight_sed_g","Water_added_g", "Water_total_g")
  

  data_location_subset = merged[grep(location[i],merged$Sample_Name),]
  
  unique.incubations = unique(data_location_subset$Sample_Name)
  
  
  for (j in 1:length(unique.incubations)){
    
    data_site_subset = subset(data_location_subset, data_location_subset$Sample_Name == unique.incubations[j])
    
    data_site_subset <- data_site_subset[with(data_site_subset, order(Date)),]
    
    data_site_subset = data_site_subset %>%
      mutate(Tare_weight_g = first(Tare_weight_g))
    
    data_site_subset = data_site_subset %>% 
      mutate(Dry_weight_sed_g = (first(Sample_weight_g) - Tare_weight_g)*(1/average_wet_g_by_dry_g)) 
  
  
    merge_dates= as.data.frame(matrix(NA, ncol = 8, nrow = length(unique(all_moisture$Sample_Name))))
    
    colnames(merge_dates) = c("Sample_Name","Date",  "Tare_weight_g", "Sample_weight_g","Sample_weight_Fill_g", "Dry_weight_sed_g", "Water_added_g", "Water_total_g")
    
  for (k in 1:nrow(data_site_subset)){
    
    if (data_site_subset$Sample_weight_Fill_g[k] > -9999) {
     
      merge_dates$Sample_Name[k] = as.character(data_site_subset$Sample_Name[k])
    
    merge_dates$Date[k] = as.character(data_site_subset$Date[k]) 
    
    merge_dates$Tare_weight_g[k] = as.numeric(data_site_subset$Tare_weight_g[k]) 
    
    merge_dates$Sample_weight_g[k] = as.numeric(data_site_subset$Sample_weight_g[k] - data_site_subset$Tare_weight_g[k])
    
    merge_dates$Sample_weight_Fill_g[k] = as.numeric(data_site_subset$Sample_weight_Fill_g[k] - data_site_subset$Tare_weight_g[k])
    
    merge_dates$Dry_weight_sed_g[k] = as.numeric(data_site_subset$Dry_weight_sed_g[k])
    

    merge_dates$Water_added_g[k] = as.numeric(data_site_subset$Sample_weight_Fill_g[k] - data_site_subset$Sample_weight_g[k])
    
  merge_dates$Water_total_g[k] = merge_dates$Sample_weight_Fill_g[k] - merge_dates$Dry_weight_sed_g[k]
  

    }
    
    else {
      
      merge_dates$Sample_Name[k] = as.character(data_site_subset$Sample_Name[k])
      
      merge_dates$Date[k] = as.character(data_site_subset$Date[k]) 
      
      merge_dates$Tare_weight_g[k] = as.numeric(data_site_subset$Tare_weight_g[k]) 
      
      merge_dates$Sample_weight_g[k] = as.numeric(data_site_subset$Sample_weight_g[k] - data_site_subset$Tare_weight_g[k])
      
      merge_dates$Water_added_g[k] = as.numeric(0)
      
      merge_dates$Dry_weight_sed_g[k] = as.numeric(data_site_subset$Dry_weight_sed_g[k])
      
      merge_dates$Av_dry_weight_sed_g[k] = as.numeric(data_site_subset$Av_dry_weight_sed_g[k])
      
      merge_dates$Water_total_g[k] = merge_dates$Sample_weight_g[k] - merge_dates$Dry_weight_sed_g[k]
      
      merge_dates$Av_water_total_g[k] = merge_dates$Sample_weight_Fill_g[k] - merge_dates$Av_dry_weight_sed_g[k]
      
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

final_wet_dry <- wet_dry %>% 
  dplyr::select(-c(Tare_weight_g)) %>% 
  rename(Wet_Sediment_Mass_g = Sample_weight_g) %>% 
  rename(Wet_Sediment_Mass_Added_Water_g = Sample_weight_Fill_g) %>% 
  rename(Dry_Sediment_Mass_g = Dry_weight_sed_g) %>% 
  rename(Total_Water_Mass_g = Water_total_g) %>% 
  rename(Added_Water_Mass_g = Water_added_g) %>% 
  distinct(Sample_Name, Date, .keep_all = TRUE) %>% 
  mutate(across(c("Dry_Sediment_Mass_g", "Total_Water_Mass_g"),round,2)) %>% 
  mutate(Wet_Sediment_Mass_Added_Water_g= if_else(is.na(Wet_Sediment_Mass_Added_Water_g), -9999, Wet_Sediment_Mass_Added_Water_g)) %>% 
  mutate(Added_Water_Mass_g= if_else(Added_Water_Mass_g == 0, -9999, Added_Water_Mass_g)) %>% 
  mutate(Dry_Sediment_Mass_g = if_else(is.na(Dry_Sediment_Mass_g), -9999, Dry_Sediment_Mass_g)) %>% 
  mutate(Total_Water_Mass_g = if_else(is.na(Total_Water_Mass_g), -9999, Total_Water_Mass_g)) 


write.csv(final_wet_dry,paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/INC/03_ProcessedData/EV_Drying_Masses_ReadyForBoye_",Sys.Date(),".csv"), row.names = F)  

ggplot(final_wet_dry, aes(x = Dry_Sediment_Mass_g))+
  geom_histogram()


location = c("-W1", "-W2", "-W3", "-W4", "-W5", "-D1", "-D2", "-D3", "-D4", "-D5")

water_mass <- as.data.frame(matrix(NA, ncol = 4, nrow =1))

colnames(water_mass) = c("Sample_Name", "Initial_Water_mass_g", "Final_Water_mass_g", "Dry_Sediment_Mass_g")

summary_file = as.data.frame(matrix(NA, ncol = 4, nrow = length(unique(final_wet_dry$Sample_Name))))

colnames(summary_file) = c("Sample_Name", "Initial_Water_mass_g", "Final_Water_mass_g", "Dry_Sediment_Mass_g")

for (i in 1:length(location)){
  
  summary_file = as.data.frame(matrix(NA, ncol = 4, nrow = length(unique(final_wet_dry$Sample_Name))))
  
  colnames(summary_file) = c("Sample_Name", "Initial_Water_mass_g", "Final_Water_mass_g", "Dry_Sediment_Mass_g")
  
  summary_location_subset = final_wet_dry[grep(location[i],final_wet_dry$Sample_Name),]
  
  unique.incubations = unique(summary_location_subset$Sample_Name)
  
  for (j in 1:length(unique.incubations)){
    
    summary_site_subset = subset(summary_location_subset, summary_location_subset$Sample_Name == unique.incubations[j])
    
    summary_site_subset <- summary_site_subset[with(summary_site_subset, order(Date)),]
    
    summary_site_subset <- summary_site_subset %>% 
      mutate(Initial_Water_mass_g = first(Total_Water_Mass_g))
    
    summary_site_subset <- summary_site_subset %>% 
      mutate(Final_Water_mass_g = last(Total_Water_Mass_g))
    
    summary_site_subset <- summary_site_subset %>% 
      select(c(Sample_Name, Initial_Water_mass_g, Final_Water_mass_g, Dry_Sediment_Mass_g)) %>% 
      filter(row_number()==1)
    
    summary_file <- rbind(summary_file, summary_site_subset)
    
  }
  
  water_mass = rbind(water_mass, summary_file)
  
}  
    
water_mass <- water_mass %>% 
  drop_na(Sample_Name) %>% 
  mutate(across(c("Initial_Water_mass_g", "Final_Water_mass_g", "Dry_Sediment_Mass_g"),round,2)) %>% 
  rename(Initial_Water_Mass_g = Initial_Water_mass_g) %>% 
  rename(Final_Water_Mass_g = Final_Water_mass_g)

write.csv(water_mass,paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/INC/03_ProcessedData/ECA_Drying_Masses_Summary_merged_by_",pnnl.user,"_on_",Sys.Date(),".csv"), row.names = F)

wet_wt <- paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/MOI/03_ProcessedData/EC_Moisture_Content_2022.csv")

corr <- read_csv(wet_wt)

mean_wet_wt <- corr %>% 
  separate(col = sample_name, into = c("Project", "kit", "analysis"), sep = "_") %>%  
  separate(col = analysis, into = c("Analysis", "Replicate"), sep = "-") %>% 
  group_by(kit) %>% 
  summarise(mean_wet_grav = (mean(percent_water_content_wet)/100))

pivot = wet_dry %>% 
  pivot_longer(cols = ends_with("_g"),
               names_to = "Sample_weight_type",
               values_to = "Sample_weight_g")

grav <- wet_dry %>% 
  group_by(Sample_Name) %>% 
  mutate(Sample_weight_initial_g = first(Sample_weight_g)) %>% 
  mutate(Sample_weight_final_g = last(Sample_weight_g)) %>% 
  ungroup() %>% 
  separate(col = Sample_Name, into = c("Project", "kit", "analysis"), sep = "_") %>%  
  separate(col = analysis, into = c("Analysis", "Replicate"), sep = "-")

final_grav <- merge(grav, mean_wet_wt, by = "kit")

all_grav <- final_grav %>% 
  mutate(mass_sed = (Sample_weight_initial_g - (Sample_weight_initial_g * mean_wet_grav))) %>% 
  mutate(mass_water_initial = (Sample_weight_initial_g * mean_wet_grav) + Water_added_initial_g) %>%  
  mutate(grav_dry_initial = mass_water_initial/mass_sed) %>% 
  mutate(mass_water_final_g = (Sample_weight_final_g - mass_sed)+Water_added_initial_g) %>% 
  mutate(grav_dry_final = mass_water_final_g/mass_sed) %>% 
  mutate(lost_grav_perc = grav_dry_initial - grav_dry_final)
  
  
all_grav_final <- all_grav %>% 
  mutate(Treat = case_when(
    endsWith(Replicate, "W1") ~ "Wet",
    endsWith(Replicate, "W2") ~ "Wet",
    endsWith(Replicate, "W3") ~ "Wet",
    endsWith(Replicate, "W4") ~ "Wet", 
    endsWith(Replicate, "W5") ~ "Wet",
    endsWith(Replicate, "D1") ~ "Dry",
    endsWith(Replicate, "D2") ~ "Dry",
    endsWith(Replicate, "D3") ~ "Dry",
    endsWith(Replicate, "D4") ~ "Dry",
    endsWith(Replicate, "D5") ~ "Dry"
  )) %>% 
  unite(kit_Treat, kit, Treat) %>% 
  select(-c(Project, Analysis, Date, Tare_weight_g, Sample_weight_g)) %>% 
  distinct(kit_Treat, Replicate, .keep_all = TRUE)

average_grav <- all_grav_final %>% 
  group_by(kit_Treat) %>% 
  mutate(average_grav_intial = mean(grav_dry_initial)) %>% 
  mutate(average_grav_final = mean(grav_dry_final)) %>% 
  mutate(average_grav_lost_subt = average_grav_intial - average_grav_final) %>% 
  mutate(average_grav_lost = mean(lost_grav_perc)) %>% 
  distinct(kit_Treat, .keep_all = TRUE) %>% 
  dplyr::select(c(kit_Treat,average_grav_intial,average_grav_final,average_grav_lost))

average_grav_lost <- average_grav %>% 
  separate(col = kit_Treat, into = c("kit", "Treat"), sep = "_") %>% 
  group_by(kit) %>% 
mutate(grav_final_diff = (average_grav_final[Treat == "Wet"] - average_grav_final[Treat == "Dry"])) %>% 
  dplyr::select(c(kit, grav_final_diff)) %>% 
  distinct(kit, .keep_all = TRUE)

## EV Dry Weights ####

dry_wt <- read.csv(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/MOI/03_ProcessedData/EV_Moisture_Content_2023.csv"))

moisture <- paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/INC/01_RawData/2023_Data_Raw_INC_EV.xlsx")

all_moisture <- read_xlsx(moisture)

all_moisture$Date <- as.Date(all_moisture$Date)

all_moisture$Sample_weight_Fill_g <- as.numeric(all_moisture$Sample_weight_Fill_g)

mean_dry_wt <- dry_wt %>% 
  mutate(wet_g_by_dry_g = wet_weight_g/true_dry_weight_g) %>% 
  separate(Sample_Name, sep = "-", c("Sample_Name", "replicate")) %>% 
  group_by(Sample_Name) %>% 
  mutate(average = mean(wet_g_by_dry_g)) %>% 
  mutate(cv = sd(wet_g_by_dry_g)/mean(wet_g_by_dry_g)) %>% 
  separate(Sample_Name, sep = "_", c("Study Code", "Site"))

all_moisture <- all_moisture %>% 
  separate(Sample_Name, into = c("EC", "Site", "INC"), sep = "_", remove = FALSE)

merged <- merge(all_moisture, mean_dry_wt, by = "Site")

merged_clean <- merged %>% 
  dplyr::select(-c(Notes, ...7 , ...8, ...9, ...10,...11, ...12, ...13, ...14, ...15, ...16, `Study Code`, percent_water_content_wet, percent_water_content_dry, true_dry_weight_g, wet_weight_g)) %>% 
  unite("Sample_Name", c("EC", "Site", "INC"), sep = "_")%>% 
  rename(Average_Wet_g_by_Dry_g = average) %>% 
  distinct(Sample_Name, Date, .keep_all = TRUE)

##All Samples - wet and dry (all dates) ####

location = c("-W1", "-W2", "-W3", "-W4", "-W5", "-D1", "-D2", "-D3", "-D4", "-D5")

wet_dry <- as.data.frame(matrix(NA, ncol = 8, nrow =1))

colnames(wet_dry) = c("Sample_Name","Date",  "Tare_weight_g", "Sample_weight_g","Sample_weight_Fill_g", "Dry_weight_sed_g", "Water_added_g", "Water_total_g")

all_dates = as.data.frame(matrix(NA, ncol = 8, nrow = length(unique(all_moisture$Sample_Name))))

colnames(all_dates) = c("Sample_Name","Date",  "Tare_weight_g", "Sample_weight_g","Sample_weight_Fill_g", "Dry_weight_sed_g","Water_added_g", "Water_total_g")

for (i in 1:length(location)){
  
  all_dates= as.data.frame(matrix(NA, ncol = 8, nrow = length(unique(merged_clean$Sample_Name))))
  
  colnames(all_dates) = c("Sample_Name","Date",  "Tare_weight_g", "Sample_weight_g","Sample_weight_Fill_g", "Dry_weight_sed_g","Water_added_g", "Water_total_g")
  
  data_location_subset = merged_clean[grep(location[i],merged_clean$Sample_Name),]
  
  unique.incubations = unique(data_location_subset$Sample_Name)
  
  
  for (j in 1:length(unique.incubations)){
    
    data_site_subset = subset(data_location_subset, data_location_subset$Sample_Name == unique.incubations[[j]])
    
    data_site_subset <- data_site_subset[with(data_site_subset, order(Date)),]
    
    data_site_subset = data_site_subset %>%
      mutate(Tare_weight_g = first(Tare_weight_g)) %>% 
      mutate(Dry_weight_sed_g = (first(Sample_weight_g) - Tare_weight_g)*(1/Average_Wet_g_by_Dry_g))
    
    merge_dates= as.data.frame(matrix(NA, ncol = 8, nrow = length(unique(all_moisture$Sample_Name))))
    
    colnames(merge_dates) = c("Sample_Name","Date",  "Tare_weight_g", "Sample_weight_g","Sample_weight_Fill_g", "Dry_weight_sed_g", "Water_added_g", "Water_total_g")
    
    for (k in 1:nrow(data_site_subset)){
      
      if (data_site_subset$Sample_weight_Fill_g[k] > -9999) {
        
        merge_dates$Sample_Name[k] = as.character(data_site_subset$Sample_Name[k])
        
        merge_dates$Date[k] = as.character(data_site_subset$Date[k]) 
        
        merge_dates$Tare_weight_g[k] = as.numeric(data_site_subset$Tare_weight_g[k]) 
        
        merge_dates$Sample_weight_g[k] = as.numeric(data_site_subset$Sample_weight_g[k] - data_site_subset$Tare_weight_g[k])
        
        merge_dates$Sample_weight_Fill_g[k] = as.numeric(data_site_subset$Sample_weight_Fill_g[k] - data_site_subset$Tare_weight_g[k])
        
        merge_dates$Dry_weight_sed_g[k] = round(as.numeric(data_site_subset$Dry_weight_sed_g[k]), 2)
        
        merge_dates$Water_added_g[k] = as.numeric(data_site_subset$Sample_weight_Fill_g[k] - data_site_subset$Sample_weight_g[k])
        
        merge_dates$Water_total_g[k] = round((merge_dates$Sample_weight_Fill_g[k] - merge_dates$Dry_weight_sed_g[k]), 2)
        
         }
      
      else {
        
        merge_dates$Sample_Name[k] = as.character(data_site_subset$Sample_Name[k])
        
        merge_dates$Date[k] = as.character(data_site_subset$Date[k]) 
        
        merge_dates$Tare_weight_g[k] = as.numeric(data_site_subset$Tare_weight_g[k]) 
        
        merge_dates$Sample_weight_g[k] = as.numeric(data_site_subset$Sample_weight_g[k] - data_site_subset$Tare_weight_g[k])
        
        merge_dates$Water_added_g[k] = as.numeric(0)
        
        merge_dates$Dry_weight_sed_g[k] = as.numeric(data_site_subset$Dry_weight_sed_g[k])
        
        merge_dates$Water_total_g[k] = merge_dates$Sample_weight_g[k] - merge_dates$Dry_weight_sed_g[k]
      
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

final_wet_dry <- wet_dry %>% 
  dplyr::select(-c(Tare_weight_g)) %>% 
  rename(Wet_Sediment_Mass_g = Sample_weight_g) %>% 
  rename(Wet_Sediment_Mass_Added_Water_g = Sample_weight_Fill_g) %>% 
  rename(Dry_Sediment_Mass_g = Dry_weight_sed_g) %>% 
  rename(Total_Water_Mass_g = Water_total_g) %>% 
  rename(Added_Water_Mass_g = Water_added_g) 

write.csv(final_wet_dry,paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/INC/03_ProcessedData/EV_Drying_Masses_merged_by_",pnnl.user,"_on_",Sys.Date(),".csv"), row.names = F)  

location = c("-W1", "-W2", "-W3", "-W4", "-W5", "-D1", "-D2", "-D3", "-D4", "-D5")

water_mass <- as.data.frame(matrix(NA, ncol = 4, nrow =1))

colnames(water_mass) = c("Sample_Name", "Initial_Water_Mass_g", "Final_Water_Mass_g", "Dry_Sediment_Mass_g")

summary_file = as.data.frame(matrix(NA, ncol = 4, nrow = length(unique(final_wet_dry$Sample_Name))))

colnames(summary_file) = c("Sample_Name", "Initial_Water_Mass_g", "Final_Water_Mass_g", "Dry_Sediment_Mass_g")

for (i in 1:length(location)){
  
  summary_file = as.data.frame(matrix(NA, ncol = 4, nrow = length(unique(final_wet_dry$Sample_Name))))
  
  colnames(summary_file) = c("Sample_Name", "Initial_Water_Mass_g", "Final_Water_Mass_g", "Dry_Sediment_Mass_g")
  
  summary_location_subset = final_wet_dry[grep(location[i],final_wet_dry$Sample_Name),]
  
  unique.incubations = unique(summary_location_subset$Sample_Name)
  
  for (j in 1:length(unique.incubations)){
    
    summary_site_subset = subset(summary_location_subset, summary_location_subset$Sample_Name == unique.incubations[j])
    
    summary_site_subset <- summary_site_subset[with(summary_site_subset, order(Date)),]
    
    summary_site_subset <- summary_site_subset %>% 
      mutate(Initial_Water_Mass_g = round(first(Total_Water_Mass_g), 2))
    
    summary_site_subset <- summary_site_subset %>% 
      mutate(Final_Water_Mass_g = round(last(Total_Water_Mass_g), 2))
    
    summary_site_subset <- summary_site_subset %>% 
      select(c(Sample_Name, Initial_Water_Mass_g, Final_Water_Mass_g, Dry_Sediment_Mass_g)) %>% 
      filter(row_number()==1)
    
    summary_file <- rbind(summary_file, summary_site_subset)
    
  }
  
  water_mass = rbind(water_mass, summary_file)
  
}  

water_mass <- water_mass %>% 
  drop_na(Sample_Name) %>% 
  filter(Dry_Sediment_Mass_g > 0) %>% 
  mutate(Dry_Sediment_Mass_g = round(Dry_Sediment_Mass_g, 2))

write.csv(water_mass,paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/INC/03_ProcessedData/EV_Drying_Masses_Summary_merged_by_",pnnl.user,"_on_",Sys.Date(),".csv"), row.names = F)
