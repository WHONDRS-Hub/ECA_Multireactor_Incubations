# library(lubridate);library(writexl);library(raster);library(devtools)
# library(readxl)
# library(corrplot)
# library(corrr)
# library(vegan)
# library(FactoMineR)
# library(factoextra)
library(readxl)
library(tidyverse)
library(dplyr)

rm(list=ls());graphics.off()

# Set working directory to data file
#Example:
pnnl.user = 'laan208'

#Read in all data
setwd(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/"))

#All incubation pH, SpC, temp
chemistry <- paste0("Optode multi reactor/Optode_multi_reactor_incubation/")

import_data = function(chemistry){ 
  
  #pull a list of files in target folder with correct pattern
  #read all files and combine
  
  map.file <-  list.files(chemistry, recursive = T, pattern = "\\.xlsx$",full.names = T)
  
  map.file <- map.file[!grepl("red|Red|EC_01_|EC_02_|EC_03_|EC_06_07|EC_10_15|EC_04_08|fast|combined|SPC|IC|QA", map.file)]
  
  mapping <- lapply(map.file, read_xlsx)
  
  for (i in 1:length(mapping)){mapping[[i]] <- cbind(mapping[[i]], map.file[i])}
  
  all.map <- 
    do.call(rbind,mapping)
}

map = import_data(chemistry)

all_chem <- map %>% 
  select(c(Sample_ID, SpC, Temp, pH, Notes)) %>% 
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
  mutate(CV_Temp = (sd(Temp)/mean(Temp))*100)
