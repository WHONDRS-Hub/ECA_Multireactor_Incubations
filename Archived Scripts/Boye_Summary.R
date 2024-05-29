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

# ================================= Wrangle Not Summarized Data =================

NPOC_TN_file <- list.files(boye_dir, 'NPOC_TN_Boye', full.names = T, recursive = T)

NPOC_TN_file <- NPOC_TN_file[!grepl('Archive',NPOC_TN_file)]

qaqc_files <- list.files(boye_dir, 'CombinedQAQC', full.names = T, recursive = T)

qaqc_files <- qaqc_files[!grepl('Archive',qaqc_files)]


## NPOC/TN ####

if (length(NPOC_TN_file) > 0) {
  
  NPOC_TN_boye_headers <- read_csv(NPOC_TN_file, n_max = 11, skip = 2)%>%
    select(-'Methods_Deviation')
 
  npoc_tn <- read_csv(NPOC_TN_file, skip = 2)
  
 # NPOC_TN_NPOC_qaqc <- qaqc_files[grepl("NPOC_TN", qaqc_files)] %>%
  #  read_csv() %>%
  #  filter(NPOC_Outlier == T)
  
 # NPOC_TN_TN_qaqc <- qaqc_files[grepl("NPOC_TN", qaqc_files)] %>%
  #  read_csv() %>%
  #  filter(TN_Outlier == T)
  
  NPOC_TN_data <- read_csv(NPOC_TN_file,  na = '-9999') %>%
    filter(!Sample_ID %in% c('N/A', '-9999', NA),
           Field_Name != '#End_Data') %>%
    mutate(Field_Name = 'N/A',
           `00681_NPOC_mg_per_L_as_C` = ifelse(Sample_Name %in% NPOC_TN_NPOC_qaqc$Sample_ID, NA, as.numeric(`00681_NPOC_mg_per_L_as_C`)),
           `00602_TN_mg_per_L_as_N` = ifelse(Sample_Name %in% NPOC_TN_TN_qaqc$Sample_ID, NA, as.numeric(`00602_TN_mg_per_L_as_N`)),
           Sample_Name = str_remove(Sample_Name, '-1'),
           Sample_Name = str_remove(Sample_Name, '-2'),
           Sample_Name = str_remove(Sample_Name, '-3'),
           Sample_Name = str_remove(Sample_Name, '_OCN')
    )
  
  NPOC_TN_summary <- NPOC_TN_data %>%
    group_by(Sample_Name) %>%
    mutate(count_NPOC = sum(!is.na(`00681_NPOC_mg_per_L_as_C`)),
           count_TN = sum(!is.na(`00602_TN_mg_per_L_as_N`))) %>%
    summarize(
      Field_Name = NA,
      Material = unique(Material),
      `Mean_00681_NPOC_mg_per_L_as_C` = mean(`00681_NPOC_mg_per_L_as_C`, na.rm = T),
      `Mean_00602_TN_mg_per_L_as_N` = mean(`00602_TN_mg_per_L_as_N`, na.rm = T),
      Mean_Missing_Reps = ifelse(count_NPOC<3, TRUE, FALSE),
      Mean_Missing_Reps = ifelse(count_TN<3, TRUE, Mean_Missing_Reps),
      count_NPOC = unique(count_NPOC),
      count_TN = unique(count_TN)
    ) %>%
    filter(!is.na(Sample_Name)) %>%
    select(Field_Name, Sample_Name, Material, `Mean_00681_NPOC_mg_per_L_as_C`, `Mean_00602_TN_mg_per_L_as_N`, Mean_Missing_Reps)%>%
    distinct()
  
  combine <- combine %>%
    full_join(NPOC_TN_summary, by = c("Field_Name", "Sample_Name", "Material")) %>%
    arrange(Sample_Name) %>%
    unite(Mean_Missing_Reps, Mean_Missing_Reps.x, Mean_Missing_Reps.y, remove = T, na.rm = T)%>%
    mutate(Mean_Missing_Reps = ifelse(str_detect(Mean_Missing_Reps, 'TRUE'), TRUE, FALSE))%>%
    filter(!is.na(Sample_Name))
  
  combine_headers <- combine_headers %>%
    left_join(NPOC_TN_boye_headers)
  
}


# ================================= Wrangle Summarized Data =================

analyte_files <- list.files(boye_dir, 'ReadyForBoye', full.names = T, recursive = T)

analyte_files <- analyte_files[grepl('Summary',analyte_files)]



# ====================== create combined data frames ===========================

atp_file = analyte_files[grepl("ATP", analyte_files)]
  
atp <- read_csv(atp_file) %>%
  mutate(Sample_Name = str_remove(Sample_Name, "_ATP")) 

sfe_file = analyte_files[grepl("SFE", analyte_files)]

sfe <- read_csv(sfe_file) %>%
  mutate(Sample_Name = str_remove(Sample_Name, "_SFE")) 

inc_file = analyte_files[grepl("Respiration", analyte_files)]

inc <- read_csv(inc_file) %>%
  mutate(Sample_Name = str_remove(Sample_Name, "_INC"))   

combine = left_join(inc, atp, by = c("Sample_Name", "Material")) %>% 
  left_join(sfe, by = c("Sample_Name", "Material")) %>% 
  mutate(Missing_Reps = ifelse(!is.na(count.x ), TRUE, ifelse(!is.na(count.y), TRUE, FALSE)))

