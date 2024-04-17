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

output.path = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/rates/Raw_DO_by_Min/")

#setwd(input.path)

#path for reading in 100% saturation values for each kit based on pressure/temperature during disk calibration
fast.sat <- paste0("C:/Users/", pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/rates/fast_rate_calculations.xlsx")

fast.rates <- read_excel(fast.sat)  

fast.rates.kits <- fast.rates %>% 
  rename("DO_mg_L" = "DO_sat_mg_L") 


df = df %>% mutate(source = basename(raw.data)) %>% 
  select(-c(software_version)) %>% 
  filter(!is.na(x3)) %>% 
  filter(!is.na(x3_05_11)) %>% 
  filter(!is.na(x15))