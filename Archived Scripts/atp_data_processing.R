# atp_data_processing.R -----------------------------------------------

# This script is used to process multiple days of ATP runs 

# 1. Fill out user inputs
# 2. User inputs are used to read in specified files
# 3. The script calculates a calibration curve and LOD for each run
# 4. The script corrects sample RLU values to concentrations and accounts for dilutions
# 5. The script corrects nanomole/L to picomole/g values

# Read in libraries 
library(tidyverse);library(readxl);library(ggpmisc)

rm(list=ls());graphics.off()



# Read in Data ------------------------------------------------------------

# Set user inputs
print.files = FALSE # switch to TRUE if you would like to save csv's to respective folders
pnnl.user = 'laan208'
project = 'ECA'
study.code = 'EC'
cv.flag = 30 # Set cut off for coefficient of variation flags in replicates
# used 30% for ICON, none for ECA

# Set Data Paths
eca.data.path = ("ECA/ATP/")
cm.data.path = ("ICON_ModEx_SSS/11_ATP/")

# Define path and working directory

if (project == 'ECA'){
  
  input.path = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/", eca.data.path,"")
  
} else if (project == 'RC4' & study.code == 'CM'){
  
  input.path = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/", cm.data.path,"")
  
} 

setwd(input.path)

raw.data = paste0("01_RawData/")

processed.data = paste0("03_ProcessedData/")

# Import data
# List raw data file paths
filePaths_atp<- list.files(path = raw.data, pattern = "Data_Raw_ATP", full.names = TRUE, recursive = TRUE) 

# Bind all raw data files together
atp_data <- do.call(bind_rows, lapply(filePaths_atp, function(raw.data) {
  df <- read_xlsx(raw.data, skip = 1)%>% mutate_all(as.character) %>% janitor::clean_names()
  df = df %>% mutate(source = basename(raw.data)) # set filename as source column
  df}))

# Shorten source file to get date of run
atp_data$source <- substr(atp_data$source, start = 1, stop = 8)

## Get RLU ranges for each standard curve per day
# Set numeric columns
atp <- atp_data %>% 
  mutate(rlu = as.numeric(rlu)) %>% 
  mutate(standard_curve_mg_l = as.numeric(standard_curve_mg_l)) %>% 
  group_by(source) %>% # group by date of run
  mutate(low_beg_range = rlu[standard_curve_mg_l == 0.00][1]) %>%  # pull out 0 nM RLU from each date (start of low and medium standard curves)
  mutate(low_end_range = rlu[standard_curve_mg_l == 5][1]) %>%  # pull out 5 nM RLU from each date (end of low standard curve)
  mutate(med_end_range = rlu[standard_curve_mg_l == 100][1]) %>%  # pull out 100 nM RLU from each date (end of medium standard curve)
  mutate(dilution_factor = as.numeric(dilution_factor)) %>%   
  mutate(dried_sediment_mass_g = as.numeric(dried_sediment_15_m_l_mass_g) - as.numeric(x15_m_l_mass_g)) %>% # calculate dried sediment mass by substracting 15 mL tube mass
  mutate(edta_est_g = as.numeric(x15_m_l_sediment_edta_wet_mass_g) - as.numeric(x15_m_l_sediment_wet_mass_g_after_decanting_reagent)) %>% #estimate mass of EDTA in sediment by subtracting total sample mass (sediment + EDTA) - (sediment after decanting EDTA), used for ECA when sample weights weren't recorded (before 12/1/2022)
  mutate(water_mass_g = as.numeric(x15_m_l_sediment_edta_wet_mass_g) - as.numeric(dried_sediment_15_m_l_mass_g) - as.numeric(edta_mass_g)) %>% #estimate mass of liquid in sediment by subtracting (wet sediment + edta) - (dried sediment) - (edta mass) 
  mutate(water_est_g = as.numeric(x15_m_l_sediment_wet_mass_g_after_decanting_reagent) - as.numeric(dried_sediment_15_m_l_mass_g)) %>% # estimate mass of liquid in sediment by subtracting (wet sediment after decanting edta) - (dried sediment), used for ECA when sample weights weren't recorded (before 12/1/2022)
  mutate(volume_edta_water_g = as.numeric(water_mass_g) + as.numeric(edta_mass_g)) %>% # EDTA + liquid in sediment
  mutate(volume_edta_water_est_g = water_est_g + edta_est_g) %>% # used for ECA when sample weights weren't recorded (before 12/1/2022)
  mutate(total_volume_g = volume_edta_water_g + 2.5) %>% #  EDTA + liquid in sediment + reagent (2.5 mL)
  mutate(total_volume_est_g = volume_edta_water_est_g + 2.5) %>% # used for ECA when sample weights weren't recorded (before 12/1/2022) + 2.5 mL reagent
  select(-c(x15_m_l_mass_g, x15_m_l_sediment_edta_wet_mass_g, x15_m_l_sediment_wet_mass_g_after_decanting_reagent, eppi_sediment_edta_wet_mass_g, eppi_mass_g, dried_sediment_15_m_l_mass_g, method_status, lab_temperature_deg_c))

## This function calculates standard curves for low (0 - 5 nM) medium (0 - 100 nM), and high (100 - 500 nM) ranges ####

## Run lines in the function if you want to check calculations

## Also has lines to look at figures of standard curves

calibrate_atp_data = function(data_formatted){
  
  #Pull out standards from each day
  standards = 
    atp %>% 
    filter(grepl("standard", sample_id)) %>% 
    dplyr::select(randomized_id, rlu, standard_curve_mg_l, source) %>% 
    drop_na(standard_curve_mg_l) %>% 
    drop_na(rlu)
  
  ## Standard Images #### 
  #use to look at figures of standard curves
  
  # All standards
  # standards %>%
  #   ggplot(aes(y = standard_curve_mg_l, x = rlu, color = as.character(source)))+
  #   geom_point()+
  #   stat_poly_line()+
  #   stat_poly_eq(use_label(c("eq")), label.y = 0.9)+
  #   stat_poly_eq(label.y = 0.85)+
  #   #geom_smooth(method = "lm", se = F)+
  #   facet_wrap(~source) +
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  # 
  
  # Low standard curve (0 - 5 nM)
  # standards %>%
  #   filter(standard_curve_mg_l <= 5) %>%
  #   ggplot(aes(y = standard_curve_mg_l, x = rlu, color = as.character(source)))+
  #   geom_point()+
  #   facet_wrap(~source) +
  #   #stat_poly_line(method = 'lm',
  #                  #formula = y ~ poly(x,2))+
  #   stat_poly_line()+
  #   stat_poly_eq(#method = 'lm',
  #                #formula = y ~ poly(x,2),
  #                use_label(c("eq")), label.y = 0.9)+
  #   stat_poly_eq(#method = 'lm',
  #                #formula = y ~ poly(x,2),
  #                label.y = 0.85, rr.digits = 3)+
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  
  # Medium standard curve (0 - 100 nM)
  # standards %>%
  #   filter(between(standard_curve_mg_l, 0, 100)) %>%
  #   ggplot(aes(y = standard_curve_mg_l, x = rlu, color = as.character(source)))+
  #   geom_point()+
  #   #stat_poly_line(method = 'lm',
  #   #formula = y ~ poly(x,2))+
  #   stat_poly_line()+
  #   stat_poly_eq(#method = 'lm',
  #     #formula = y ~ poly(x,2),
  #     use_label(c("eq")), label.y = 0.9)+
  #   stat_poly_eq(#method = 'lm',
  #     #formula = y ~ poly(x,2),
  #     label.y = 0.85, rr.digits = 3)+
  #   facet_wrap(~source) +
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  # 
  
  # High standard curve (100 - 500 nM)
  # standards %>%
  #   filter(standard_curve_mg_l >= 100) %>%
  #   ggplot(aes(y = standard_curve_mg_l, x = rlu, color = as.character(source)))+
  #   geom_point()+
  #   #stat_poly_line(method = 'lm',
  #   #formula = y ~ poly(x,2))+
  #   stat_poly_line()+
  #   stat_poly_eq(#method = 'lm',
  #     #formula = y ~ poly(x,2),
  #     use_label(c("eq")), label.y = 0.9)+
  #   stat_poly_eq(#method = 'lm',
  #     #formula = y ~ poly(x,2),
  #     label.y = 0.85, rr.digits = 3)+
  #   facet_wrap(~source) +
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  
  #####
  
  # calculate slope, r2, and intercept of each calibration curve
  calibration_coef = 
    standards %>% 
    dplyr::group_by(source) %>% # group by dates
    dplyr::summarize(
      low_slope = lm(rlu ~ standard_curve_mg_l, 
                     subset = standard_curve_mg_l <= 5)$coefficients["standard_curve_mg_l"], # pull out slope of low standard curve (0 - 5 nM) 
      low_intercept = lm(rlu ~ standard_curve_mg_l, 
                         subset = standard_curve_mg_l <= 5)$coefficients["(Intercept)"], # get low standard curve intercept
      low_r2 = summary(lm(rlu ~ standard_curve_mg_l, 
                          subset = standard_curve_mg_l <= 5))$r.squared, # get low standard curve R2
      low_count = sum(standard_curve_mg_l <= 5), #count number of standards in low standard curve
      med_slope = lm(rlu ~ standard_curve_mg_l, 
                     subset = between(standard_curve_mg_l, 0, 100))$coefficients["standard_curve_mg_l"], # pull out slope of medium standard curve (0 - 100 nM) 
      med_intercept = lm(rlu ~ standard_curve_mg_l, 
                         subset = between(standard_curve_mg_l, 0, 100))$coefficients["(Intercept)"], # get medium standard curve intercept
      med_r2 = summary(lm(rlu ~ standard_curve_mg_l, 
                          subset = between(standard_curve_mg_l, 0, 100)))$r.squared, # get medium standard curve R2
      high_slope = lm(rlu ~ standard_curve_mg_l, 
                      subset = standard_curve_mg_l >= 100)$coefficients["standard_curve_mg_l"], # pull out slope of high standard curve (0 - 100 nM) 
      high_intercept = lm(rlu ~ standard_curve_mg_l, 
                          subset = standard_curve_mg_l >= 100)$coefficients["(Intercept)"], # pull out intercept of high standard curve (0 - 100 nM) 
      high_r2 = summary(lm(rlu ~ standard_curve_mg_l, 
                           subset = standard_curve_mg_l >= 100))$r.squared
    ) # pull out R2 of high standard curve (0 - 100 nM)
  
  # y = mx + c
  # abs = m*ppm + c
  # ppm = abs-c/m
  
  # Calculate nM value for low, medium, and high range of standards
  data_formatted = atp %>% 
    left_join(calibration_coef) %>% 
    mutate(low_ppm_calculate = (rlu - low_intercept) / low_slope) %>% 
    mutate(med_ppm_calculate = (rlu - med_intercept) / med_slope ) %>% 
    mutate(high_ppm_calculate = (rlu - high_intercept) / high_slope ) 
  
}

## Check QC values from throughout runs
qa = 
  calibrate_atp_data(data_formatted) %>% 
  filter(grepl("standard", sample_id)) %>% # pull out only standards
  filter(is.na(standard_curve_mg_l)) %>% 
  mutate(actual_nM = if_else(rlu <= low_end_range, low_ppm_calculate * dilution_factor,
                             if_else(rlu > low_end_range & rlu <= med_end_range, med_ppm_calculate * dilution_factor, high_ppm_calculate * dilution_factor))) %>% #assign correct standard curve based on RLU ranges set earlier in script
  relocate(actual_nM, .after = rlu) %>% 
  relocate(source, .after = actual_nM) %>% 
  drop_na(actual_nM) %>% 
  select(c(sample_id, randomized_id, rlu, actual_nM, source)) %>% 
  group_by(randomized_id) %>% # group by standard value (ie 5 nM)
  mutate(cv_all = (sd(actual_nM)/mean(actual_nM))*100) %>% #Check CV of standards across all runs
  ungroup() %>% 
  group_by(source, randomized_id) %>% 
  mutate(cv_day = (sd(actual_nM)/mean(actual_nM))*100) #Check CV of standards for each day

## Calculate LOD from each day using 0 nM QC values
#range from +/- 0.04
blanks <- qa %>%
  filter(randomized_id == "0_std") %>% 
  group_by(source) %>% 
  mutate(day_lod_sd = mean(actual_nM) + 3*sd(actual_nM)) %>% # calculate  LOD as mean of blanks in run + 3*sd of nM calculation of blanks in run
  mutate(range = max(actual_nM) - min(actual_nM)) %>% # calculate range of blanks
  ungroup()

# Pull out LOD for each day to assign strings
lod <-  blanks %>% 
  select(c(source, day_lod_sd)) %>% 
  distinct(source, .keep_all = TRUE) %>% 
  mutate(day_lod_sd = round(day_lod_sd, 3))

## Check reference sediments throughout runs
# For CM and ECA, CV of reference sediments across all runs ~14%

reference = 
  calibrate_atp_data(data_formatted) %>% # calculate nM values of reference sediments
  ungroup() %>% 
  filter(grepl("PE", sample_id)) %>% # pull out reference sediments
  mutate(pmol_per_g = (med_ppm_calculate * (((total_volume_g/1000))/dried_sediment_mass_g))*1000) %>% #calculate pmol/g value by multiplying nmol * (total volume of liquid in sample (EDTA + water + 2.5 mL reagent) /dried sediment mass
  mutate(cv_all = (sd(pmol_per_g)/mean(pmol_per_g))*100) %>% # calculate CV of all reference sediments across days
  group_by(source) %>% 
  mutate(cv_day = (sd(pmol_per_g)/mean(pmol_per_g))*100) %>% # calculate CV of reference sediments each day 
  ungroup() %>% 
  select(c(sample_id, dilution_factor, rlu, med_ppm_calculate, pmol_per_g, cv_all, cv_day, source))

# boxplots of reference sediments

# mean_ref = mean(reference$pmol_per_g)

# ggplot(reference, aes(y = pmol_per_g)) +
#   geom_boxplot()+
#   facet_grid(~source) +
#   geom_hline(aes(yintercept = mean_ref))


## Use calibration coefficients to calculate sample nanomole/L value and correct to picomole/g value
samples_nm = 
  calibrate_atp_data(data_formatted) %>% 
  ungroup() %>% 
  select(-c(randomized_id, standard_curve_mg_l, n_m_calculation)) %>% 
  filter(grepl(study.code, sample_id)) %>% 
  mutate(actual_nM = if_else(rlu <= low_end_range, low_ppm_calculate * dilution_factor, if_else(rlu > low_end_range & rlu <= med_end_range, med_ppm_calculate * dilution_factor, high_ppm_calculate * dilution_factor))) %>% # if sample RLU < 5 nM standard (low curve), multiply calculated low nM * dilution factor, if sample RLU > 5 nM and < 100 nM (medium curve), multiply calculated medium nM * dilution factor, if sample RLU >100 nM (high curve), multiply calculated high nM * dilution factor 
  mutate(sc_used = if_else(rlu <= low_end_range, "low", if_else(rlu > low_end_range & rlu <= med_end_range, "med", "high"))) %>% # notate which standard curve was used
  mutate(pmol_g = (actual_nM * (((total_volume_g/1000))/dried_sediment_mass_g))*1000) %>% #calculate pmol/g ATP values by correcting for solid:solution ratio of tubes
  mutate(pmol_g_est = if_else(is.na(pmol_g), (actual_nM * (((total_volume_est_g/1000))/dried_sediment_mass_g))*1000, pmol_g)) %>%  # used for ECA when sample weights weren't recorded (before 12/1/2022)
  select(c(sample_id, dilution_factor, actual_nM, pmol_g, pmol_g_est, sc_used, source, rlu)) 

#check that samples calculated with low and medium standard curves don't overlap
# samples_nm %>%
#    filter(dilution_factor == 1) %>%
# ggplot(samples_nm, mapping = aes(y = actual_nM, x = rlu, color = sc_used))+
#   geom_point()+
#   facet_wrap(~source) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

#boxplots of sample nanomol/L value

# ggplot(samples_nm, aes(y = actual_nM, x = Sample)) +
#   geom_boxplot()

#boxplots of sample picomol/g value

# ggplot(samples_nmol, aes(y = pmol_g, x = Sample)) + 
#   geom_boxplot()

# Flag samples below LOD 
samples_lod <- left_join(samples_nm,lod, by = "source") %>%  
  mutate(actual_nM = round(actual_nM, 2)) %>% 
  mutate(pmol_g = round(pmol_g, 2)) %>% 
  mutate(actual_nM_lod = if_else(actual_nM <= day_lod_sd, paste0("ATP_Below_",day_lod_sd,"_nanomol_per_L_LOD|", actual_nM, "_nanomol_per_L_Raw_Not_Corrected|", actual_nM,"_nanomol_per_L_Final_Corrected"), as.character(actual_nM))) %>% 
  mutate(pmol_lod = if_else(actual_nM <= day_lod_sd, paste0("ATP_Below_",day_lod_sd,"_nanomol_per_L_LOD|", actual_nM, "_nanomol_per_L_Raw_Not_Corrected|", pmol_g,"_picomol_per_g_Final_Corrected"), as.character(pmol_g)))

samples_final <- samples_lod %>% 
  rename(ATP_nanomol_per_L = actual_nM_lod) %>% 
  rename(ATP_picomol_per_g = pmol_lod) %>% 
  mutate(ATP_nanomol_per_L = if_else(rlu == -9999, "-9999", ATP_nanomol_per_L)) %>% 
  mutate(ATP_nanomol_per_L = if_else(is.na(ATP_nanomol_per_L), "-9999", ATP_nanomol_per_L)) %>% 
  mutate(ATP_picomol_per_g = if_else(rlu == -9999, "-9999", ATP_picomol_per_g)) %>% 
  mutate(ATP_picomol_per_g = if_else(is.na(ATP_picomol_per_g), "-9999", ATP_picomol_per_g)) %>%
  mutate(Methods_Deviation = if_else(dilution_factor == 2, "ATP_001", "N/A")) %>% # Add deviation if samples were diluted
  rename(Sample_Name = sample_id) %>% 
  dplyr::select(c(Sample_Name, ATP_nanomol_per_L, ATP_picomol_per_g, Methods_Deviation))


# Add method deviation and flag for replicates with high CVs - 30% for ICON, did not flag for ECA

# Make function to calculate coefficient of variation
cv <- function(x) {
  cv_value <- (sd(x) / mean(x)) * 100
  return(cv_value)
}

# Count number of missing samples
missing = samples_final %>% 
  filter(ATP_nanomol_per_L == "-9999") %>% 
  separate(Sample_Name, c("Sample_Name", "Rep"), sep = -1) %>% 
  group_by(Sample_Name) %>% 
  summarise(Sample_Name, count = n())

# Calculate CV of samples in picomoles/g
atp_outliers = samples_final %>% 
  separate(Sample_Name, c("Sample_ID", "Rep"), sep = -1) %>% 
  filter(ATP_nanomol_per_L != -9999) %>% 
  select(c(Sample_ID, Rep, ATP_nanomol_per_L, ATP_picomol_per_g)) %>%   mutate(ATP_nanomol_per_L = as.numeric(ATP_nanomol_per_L)) %>% 
  mutate(ATP_picomol_per_g = as.numeric(ATP_picomol_per_g)) %>% 
  group_by(Sample_ID) %>% 
  mutate(ATP_picomol_per_g_cv = cv(ATP_picomol_per_g)) 

# Make new dataframe for samples being removed
atp_removed_outliers = as.data.frame(matrix(NA, ncol = 7, nrow = 1))

# Assign column names
colnames(atp_removed_outliers) = c("Sample_ID", "Rep", "ATP_nanomol_per_L", "ATP_picomol_per_g", "ATP_picomol_per_g_cv", "CV_after", "flag")

# Determine unique samples (ie CM_041)
unique.samples = unique(atp_outliers$Sample_ID)

for (i in 1:length(unique.samples)) {
  
  # Subset by sample ID
  data_subset = subset(atp_outliers, atp_outliers$Sample_ID == unique.samples[i])
  
  #matrix of pmol/g values from all analytical replicate
  conc.temp = as.numeric(data_subset$ATP_picomol_per_g)
  
  #Get st.dev., mean, CV of replicates
  conc.temp.sd <- sd(conc.temp)
  conc.temp.mean <- mean(conc.temp)
  CV = (conc.temp.sd/conc.temp.mean)*100
  
  #looping to get 2 out of 3 best samples
  for (sample.reduction in 1:3)  {
    
    if (conc.temp.mean == 0) {
      #if mean and range is 0, don't do anything
      CV = 0
      
    }
    
    else if (length(conc.temp) > 2 & CV >= cv.flag) {
      # If more than 2 analytical reps with CV > 30%, start removing samples 
      
      #Make distance matrix of an.reps.
      dist.temp = as.matrix(abs(dist(conc.temp)))
      dist.comp = numeric()
      
      for(conc.now in 1:ncol(dist.temp)) {
        # add all rows in each column together 
        dist.comp = rbind(dist.comp,c(conc.now,sum(dist.temp[,conc.now])))
        
      }
      # make sure 2nd column is numeric
      dist.comp[,2] = as.numeric(dist.comp[,2])
      conc.temp = conc.temp[-which.max(dist.comp[,2])]
      
      # recalculate st.dev., mean, cv
      conc.temp.sd <- sd(conc.temp)
      conc.temp.mean <- mean(conc.temp)
      conc.temp.cv <- (conc.temp.sd/conc.temp.mean)*100
      CV = conc.temp.cv
      
    }
    
  }
  
  if(length(conc.temp) >= 2) {
    
    if(CV > cv.flag){
      # If removals don't help meet CV and range thresholds, flag samples 
      atp_removed_outliers$ATP_picomol_per_g[which(atp_removed_outliers$Sample_Name == unique.samples[i])] = "Samples too Variable"
      
    }
    
    else {
      #make new data frame of final replicates
      samples_removed_combined = as.data.frame(conc.temp)
      # rejoin with sample names by remaining samples
      samples_remove = merge(samples_removed_combined, data_subset, by.x = "conc.temp", by.y = "ATP_picomol_per_g", all.x = TRUE)
      
      samples_remove_omit = merge(samples_remove, data_subset, by = c("Sample_ID", "Rep"), all = TRUE)
      
      # calculate new CV
      samples_remove_omit$CV_after = as.numeric(abs((sd(conc.temp)/mean(conc.temp))*100))
      
      # flag removed samples with "OMIT"
      samples_remove_omit$flag[is.na(samples_remove_omit$conc.temp)] = "OMIT"
      
      # clean up data frame
      samples_remove_omit <- samples_remove_omit %>% 
        dplyr::select(-c(ATP_picomol_per_g_cv.x, ATP_nanomol_per_L.x, conc.temp)) %>% 
        rename(ATP_picomol_per_g_cv = ATP_picomol_per_g_cv.y) %>% 
        rename(ATP_nanomol_per_L = ATP_nanomol_per_L.y)
      
    }
    
  }
  
  # create final data frame
  atp_removed_outliers = rbind(samples_remove_omit, atp_removed_outliers)
  
  rm('conc.temp')
}


samples_final_outliers = atp_removed_outliers %>% 
  unite(Sample_Name, c("Sample_ID", "Rep"), sep = "", remove = FALSE) %>% # add replicate back to name
  left_join(samples_final, by = c("Sample_Name")) %>% # join with samples_final data frame to add Methods_Deviations
  select(Sample_Name, Sample_ID, ATP_nanomol_per_L.y, ATP_picomol_per_g.y, flag, Methods_Deviation) %>% 
  rename(ATP_nanomol_per_L = ATP_nanomol_per_L.y) %>% 
  rename(ATP_picomol_per_g = ATP_picomol_per_g.y) %>% 
  group_by(Sample_ID) %>% 
  mutate(Methods_Deviation.flag = ifelse(flag == "OMIT", "ATP_CV_030", "N/A")) %>% # Flag samples that have high CV with CV deviation
  fill(Methods_Deviation.flag, .direction = 'downup') %>% 
  mutate(Methods_Deviation_new = ifelse(is.na(Methods_Deviation.flag), Methods_Deviation, ifelse(Methods_Deviation != "N/A", paste(Methods_Deviation, Methods_Deviation.flag, sep = ";"), Methods_Deviation.flag))) %>% 
  ungroup() %>% 
  drop_na(Sample_ID) %>% 
  select(c(Sample_Name, ATP_nanomol_per_L, ATP_picomol_per_g, Methods_Deviation_new, flag)) %>% 
  rename(Methods_Deviation = Methods_Deviation_new)

if (print.files == TRUE) { 
  
  write.csv(samples_final_outliers, paste0(processed.data,"",study.code,"_ATP_ReadyForBoye_", Sys.Date(),".csv"), row.names = F) 
  
}