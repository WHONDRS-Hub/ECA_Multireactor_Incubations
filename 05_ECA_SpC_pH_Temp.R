## Pull out SpC, pH, and Temp values from ECA Mapping Files
library(readxl)
library(tidyverse)
library(dplyr)

rm(list=ls());graphics.off()

# Set user inputs
pnnl.user = 'laan208'
study.code = 'EL'

# Set working directory to data file
setwd(paste0("C:/Users/",pnnl.user,"/OneDrive - PNNL/Shared Documents - Core Richland and Sequim Lab-Field Team/Data Generation and Files/ECA/"))

#All incubation pH, SpC, temp
chemistry <- paste0("Optode multi reactor/Optode_multi_reactor_incubation/")

import_data = function(chemistry){ 
  
  #pull a list of files in target folder with correct pattern
  #read all files and combine
  
  map.file <-  list.files(chemistry, recursive = T, pattern = "Mapping.xlsx$",full.names = T)
  
  map.file <- map.file[grepl(study.code, map.file)]
  
  mapping <- lapply(map.file, read_xlsx)
  
    for (i in 1:length(mapping)){mapping[[i]] <- cbind(mapping[[i]], map.file[i])}
  
  all.map <- 
    do.call(rbind,mapping)
}

map = import_data(chemistry)


if (study.code == "EC") {
# Add extra 0 to site name in mapping files and remove blanks
all_chem <- map %>% 
  dplyr::select(c(Sample_ID, SpC, Temp, pH, Notes)) %>% 
  filter(!grepl("Blank", Sample_ID)) %>% 
  tidyr::separate(Sample_ID, into = c("EC", "kit", "rep"), sep = "_", remove = FALSE)

for (i in 1:nrow(all_chem)){
  
  if (str_count(all_chem$kit[i], "[0-9]") <= 2){
    
    all_chem$kit[i] = paste0("0", all_chem$kit[i])
    
  }
  
  else {
    
    all_chem$kit[i] = all_chem$kit[i]
  }
  
}
}


# Add method deviations and material columns, look at CV of SpC, pH, and Temperatures
all_chem_corr <- all_chem %>% 
  unite(Sample_ID, c("EC", "kit", "rep"), sep = "_") %>% 
  relocate(Sample_ID, .before = SpC) %>% 
  rename(Sample_Name = Sample_ID) %>% 
  mutate(Methods_Deviation = "N/A") %>% 
  separate(Sample_Name, into = c("EC", "kit", "rep"), remove = FALSE, sep = "_") %>% 
  separate(rep, into = c("INC", "rep"), remove = FALSE, sep = "-") %>% 
  mutate(Treat = case_when(grepl("W",rep)~"Wet",
                           grepl("D", rep) ~"Dry")) %>% 
  relocate(Treat, .after = rep) %>%
  group_by(kit, Treat) %>% 
  mutate(CV_SpC = (sd(SpC)/mean(SpC))*100) %>% 
  mutate(CV_pH = (sd(pH)/mean(pH))*100) %>% 
  mutate(CV_Temp = (sd(Temp)/mean(Temp))*100) %>% 
  mutate(Material = "Sediment") %>% 
  ungroup()

chem_final <- all_chem_corr %>% 
  dplyr::select(c(Sample_Name, Material, SpC, Temp, pH, Methods_Deviation))

# Export final file
write.csv(chem_final, file.path(paste0(chemistry,"SpC_pH_Temp_Processed_Data/",study.code,"_SpC_pH_Temp_ReadyForBoye_",Sys.Date(),".csv")), row.names = F)


