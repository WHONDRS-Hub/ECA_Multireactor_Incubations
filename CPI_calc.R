library(lubridate);library(writexl);library(raster);library(tidyverse);library(devtools)

#### Load data #####
rm(list=ls());graphics.off()

# Set working directory to data file

pnnl.user = 'laan208'

input.path <- paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/rates/Plots/All_Respiration_Rates")

setwd(input.path)

#change date to most recent respiration rate csv

date = '2023-05-25'

resp <- read.csv(paste0(input.path,"/removed_respiration_merged_by_",pnnl.user,"_on_",date,".csv"))

wet <- resp %>% 
  filter(Treat == "Wet")

wet$rate_tot = sum(as.numeric(wet$rate_mg_per_L_per_min), na.rm = TRUE)

wet$rate_median = median(as.numeric(wet$rate_mg_per_L_per_min), na.rm = TRUE)

wet$great_median_tot = sum(wet$rate_mg_per_L_per_min[wet$rate_mg_per_L_per_min > wet$rate_median])

wet$cpi = wet$great_median_tot/wet$rate_tot

dry <- resp %>% 
  filter(Treat == "Dry")

dry$rate_tot = sum(as.numeric(dry$rate_mg_per_L_per_min), na.rm = TRUE)

dry$rate_median = median(as.numeric(dry$rate_mg_per_L_per_min), na.rm = TRUE)

dry$great_median_tot = sum(dry$rate_mg_per_L_per_min[dry$rate_mg_per_L_per_min > dry$rate_median])

dry$cpi = dry$great_median_tot/dry$rate_tot
                 