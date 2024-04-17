###### Load Library ######

library(dplyr); library(ggplot2);library(ggsignif)
library(ggpubr);library(reshape2);library(ggpmisc)
library(segmented);library(broom);library(lmtest);library(car)
library(ggpmisc);library(lubridate); library(readxl);
library(tidyverse);library(patchwork)
library(readr)


##### Load data ######
rm(list=ls());graphics.off()

# Set working directory to data file
#Example:
pnnl.user = 'laan208'

input.path = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation")

map.path = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/")

output.path = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/INC/03_ProcessedData/Raw_DO_by_Min/")

#setwd(input.path)

#path for reading in 100% saturation values for each kit based on pressure/temperature during disk calibration
fast.sat <- paste0("C:/Users/", pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/rates/ev_fast_rate_calculations.xlsx")

fast.rates <- read_excel(fast.sat)  

fast.rates.kits <- fast.rates %>% 
  rename("DO_mg_L" = "DO_sat_mg_L") 

#read in respiration data and clean
import_data = function(input.path){
  
  # pull a list of file names in the target folder with the target pattern
  # then read all files and combine
  
  filePaths <- list.files(path = input.path, recursive = T, pattern = "\\.csv$", full.names = TRUE)
  
  filePaths <- filePaths[grepl("EV", filePaths)]
  
  filePaths <- filePaths[grepl("results", filePaths)]
  
   # dat <- 
  do.call(rbind, lapply(filePaths, function(input.path){
    # then add a new column `source` to denote the file name
    df <- read.csv(input.path, skip = 4)
    df[["source_file"]] <- rep(input.path, nrow(df)) # add a column for the source file
    
    df %>%
      na.omit () %>% 
      as_tibble(row.names = 1:nrow(df)) %>% 
      # since all rows are in 2-min increments, just multiply the row number by 2
      tibble::rownames_to_column() %>%
      mutate(elapsed_min = as.numeric(rowname)*2) %>%
      dplyr::select(-rowname) %>% 
      # make longer, so all the data are in a single column
      pivot_longer(-c(elapsed_min,source_file), names_to = "disc_number", values_to = "DO_mg_L", values_transform = as.numeric) %>% 
      # remove unnecessary strings from the source_file name
      mutate(source_file = str_remove_all(source_file, paste0(input.path, "/")))
  }
  ))
}
data = import_data(input.path)

##### Clean Data ####

data_long = 
  data %>% 
  mutate(disc_number = str_remove_all(disc_number, "X")) %>%
  mutate(source_file = str_remove_all(source_file, ".*Optode_multi_reactor_incubation/"),source_file = str_remove_all(source_file, "/results/results.csv")) %>%  
  filter(elapsed_min > 0)

import_data = function(map.path){ 
  
  #pull a list of files in target folder with correct pattern
  #read all files and combine
  
  map.file <-  list.files(map.path, recursive = T, pattern = "\\.xlsx$",full.names = T)
  
  map.file <- map.file[grepl("EV", map.file)]
  
  mapping <- lapply(map.file, read_xlsx)
  
  for (i in 1:length(mapping)){mapping[[i]] <- cbind(mapping[[i]], map.file[i])}
  
  all.map <- 
    do.call(rbind,mapping)
}

all.map = import_data(map.path)

all.map = all.map %>% 
  rename("source_file" = "map.file[i]") %>% 
  rename("disc_number" = "Disk_ID") %>% 
  mutate(source_file = str_remove_all(source_file, ".*//")) %>% 
  separate(source_file, sep = "/", c("source_file", "file")) %>% 
  mutate(source_file = str_replace(source_file, "EC", "results")) %>% 
  dplyr::select(-`Disk Calibration date`, -Location, -`Time on`, -`Time off`, -SpC, -pH, -Temp, -Notes, -file)

all.samples <- merge(data_long, all.map)


bind <- merge(all.samples, fast.rates.kits, all = TRUE)

bind <- bind %>% 
  filter(!grepl("Blank", Sample_ID)) %>% 
  separate(Sample_ID, into = c("EC", "Kit", "Rep"), sep = "_")

cleaned_data <- bind %>% 
  unite(Sample_Name, c("EC", "Kit", "Rep"), sep = "_") %>% 
  relocate(Sample_Name, .before = elapsed_min) %>% 
  relocate(DO_mg_L, .before = elapsed_min) %>% 
  mutate(DO_mg_per_L = DO_mg_L) %>% 
  dplyr::select(c(Sample_Name, DO_mg_per_L, elapsed_min)) %>% 
  arrange(Sample_Name) %>% 
  mutate(Methods_Deviation = "N/A") %>% 
  mutate(Methods_Deviation = ifelse(elapsed_min == 0, "RATE_004", "N/A")) %>% 
  rename(Elapsed_Minute = elapsed_min)

write.csv(cleaned_data, file.path(output.path,"EV_INC_Raw_DO_By_Min_ReadyForBoye_04-17-2024.csv"), quote = F, row.names = F) 

splits <- all.samples %>% 
  #group_by(source_file) %>% 
  group_split(source_file) %>% 
  as.list() 

for (i in 1:length(splits)) {
  
  results <- splits[[i]]
  
  results <- results %>% 
    separate(Sample_ID, sep = "_", c("EC", "Site", "Rep"), remove = FALSE)
  
  unique.sites = unique(results$Site)
  unique.sites = na.omit(unique.sites) 
  
  
  for (j in 1:length(unique.sites)) {
    
    site_results <- subset(results, results$Site == unique.sites[j])
    
    for (k in 1:nrow(site_results)){
      
      if (str_count(site_results$Site[k], "[0-9]") <= 2){
        
        site_results$Site[k] = paste0("0", site_results$Site[k])
        
      }
      
      else {
        
        site_results$Site[k] = site_results$Site[k]
      }
      
    }
    
    site_no <- site_results$Site[1] 
    
    site_results <- site_results %>% 
      unite("Sample_ID",  c("EC", "Site", "Rep"), sep = "_") 
    
    site_results = site_results[order(site_results$elapsed_min, decreasing = FALSE),]
    site_results$elapsed_min = as.numeric(site_results$elapsed_min)
    
    site_results <- site_results %>% 
      dplyr::select(-c(disc_number, source_file)) %>% 
      pivot_wider(names_from = Sample_ID, values_from = c(DO_mg_L))
    
    order = sort(colnames(site_results))
    
    site_results <- site_results[, order]
    
    site_results <- site_results[,c(ncol(site_results), 1:(ncol(site_results)-1))]
    
    
    write.csv(site_results, paste0(output.path,"Raw_DO_by_Minute_EC_",site_no,".csv"), row.names = FALSE)
    
  }
  
}


