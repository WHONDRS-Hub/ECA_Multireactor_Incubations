# ==============================================================================
#
# Make a summary file of means for each analyte file going into a data package
#
# Status: In progress
#
# known issue: putting NA in detection limit and precision row 
# 
# ==============================================================================
#
# Author: Brieanne Forbes, brieanne.forbes@pnnl.gov
# 30 Sept 2022
#
# ==============================================================================

library(tidyverse)
library(janitor)
library(stringr)
rm(list=ls(all=T))

# ================================= User inputs ================================

dir <- 'C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/'

RC <- 'ECA'

study_code <- 'EC'

material <- 'Sediment'

# ====================================== Build dir name ========================

boye_dir <- paste0(dir, RC, '/')

# ================================= Wrangle data and summarize =================

analyte_files <- list.files(boye_dir, 'ReadyForBoye', full.names = T, recursive = T)

analyte_files <- analyte_files[grepl('Summary',analyte_files)]

analyte_files <- analyte_files[!grepl('Archive',analyte_files)]

# ====================== create combined data frames ===========================

atp_file = analyte_files[grepl("ATP", analyte_files)]
  
atp <- read_csv(atp_file) %>%
  mutate(Sample_Name = str_remove(Sample_Name, "_ATP")) 

sfe <- read_csv(analyte_files[2]) %>%
  mutate(Sample_Name = str_remove(Sample_Name, "_SFE")) 

inc <- read_csv(analyte_files[4]) %>%
  mutate(Sample_Name = str_remove(Sample_Name, "_INC"))   


