library(tidyverse)
library(readxl)
library(ggpmisc)

rm(list=ls());graphics.off()

pnnl.user = 'laan208'

raw.data = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ICON_ModEx_SSS/11_ATP/01_Raw_Data/")

processed.data = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ICON_ModEx_SSS/11_ATP/03_ProcessedData/")


filePaths_atp<- list.files(path = raw.data, pattern = "Data_Raw_ATP", full.names = TRUE, recursive = TRUE) 

atp_data <- do.call(bind_rows, lapply(filePaths_atp, function(raw.data) {
  df <- read_xlsx(raw.data, skip = 1)%>% mutate_all(as.character) %>% janitor::clean_names()
  df = df %>% mutate(source = basename(raw.data))
  df}))

atp_data$source <- substr(atp_data$source, start = 1, stop = 8)

## Get RLU ranges for each standard curve per day
atp <- atp_data %>% 
  mutate(rlu = as.numeric(rlu)) %>% 
  mutate(standard_curve_mg_l = as.numeric(standard_curve_mg_l)) %>% 
  group_by(source) %>% 
  mutate(low_beg_range = rlu[standard_curve_mg_l == 0.00][1]) %>% 
  mutate(low_end_range = rlu[standard_curve_mg_l == 5][1]) %>% 
  mutate(med_end_range = rlu[standard_curve_mg_l == 100][1]) %>% 
   mutate(dilution_factor = as.numeric(dilution_factor))

## Calculate standard curves for low (0 - 5 nM) medium (0 - 100 nM), and high (100 - 500 nM) ranges ####

calibrate_atp_data = function(data_formatted){
  
  standards = 
    atp %>% 
    filter(grepl("standard", sample_id)) %>% 
    dplyr::select(randomized_id, rlu, standard_curve_mg_l, source) %>% 
    drop_na(standard_curve_mg_l) %>% 
    drop_na(rlu)
   
## Standard Images ####  
  # png(file = paste0(processed.data,"all_standards_linear.png"), width = 15, height = 15, units = "in", res = 300) 
  # 
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
  # dev.off()
  # 
  # png(file = paste0(processed.data,"all_standards_nonlinear.png"), width = 15, height = 15, units = "in", res = 300) 
  # 
  # standards %>% 
  #   ggplot(aes(y = standard_curve_mg_l, x = rlu, color = as.character(source)))+
  #   geom_point()+
  #   stat_poly_line(method = 'lm', 
  #               formula = y ~ poly(x,2))+
  #   #stat_poly_line()+
  #   stat_poly_eq(method = 'lm', 
  #                formula = y ~ poly(x,2),
  #                use_label(c("eq")), label.y = 0.9)+
  #   stat_poly_eq(method = 'lm', 
  #                formula = y ~ poly(x,2),
  #                label.y = 0.85)+
  #   facet_wrap(~source) + 
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  # 
  # dev.off()
  # 
  
  # png(file = paste0(processed.data,"low_standards_linear.png"), width = 15, height = 15, units = "in", res = 300)
  # 
  # standards %>%
  #   filter(standard_curve_mg_l <= 1) %>%
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
  # 
  # dev.off()
  # 
  # 
  #  png(file = paste0(processed.data,"med_standards_linear.png"), width = 15, height = 15, units = "in", res = 300)   
  # 
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
  # dev.off()
  # 
  # png(file = paste0(processed.data,"high_standards_linear.png"), width = 15, height = 15, units = "in", res = 300)
  # 
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
  # 
  # dev.off()
#####
  
  calibration_coef = 
    standards %>% 
    dplyr::group_by(source) %>% 
    dplyr::summarize(
     # all_poly_x = lm(standard_curve_mg_l~ poly(rlu,2))$coefficients["poly(rlu, 2)1"],
     # all_poly_x2 = lm(standard_curve_mg_l~ poly(rlu,2))$coefficients["poly(rlu, 2)2"],
    # all_poly_intercept = lm(standard_curve_mg_l ~ poly(rlu,2))$coefficients["(Intercept)"], 
    # poly_r2 = summary(lm(standard_curve_mg_l ~ poly(rlu,2)))$r.squared, 
      #stat_poly_line(method = 'lm', 
                    # formula = y ~ poly(x,2))
    # all_slope = lm(standard_curve_mg_l~rlu)$coefficients["rlu"], 
    # all_intercept = lm(standard_curve_mg_l~rlu)$coefficients["(Intercept)"], 
    # all_r2 = summary(lm(standard_curve_mg_l~rlu))$r.squared,
      low_slope = lm(rlu ~ standard_curve_mg_l, 
                     subset = standard_curve_mg_l <= 5)$coefficients["standard_curve_mg_l"], 
      low_intercept = lm(rlu ~ standard_curve_mg_l, 
                     subset = standard_curve_mg_l <= 5)$coefficients["(Intercept)"], 
      low_r2 = summary(lm(rlu ~ standard_curve_mg_l, 
                  subset = standard_curve_mg_l <= 5))$r.squared, 
      low_count = sum(standard_curve_mg_l <= 5),
      med_slope = lm(rlu ~ standard_curve_mg_l, 
                     subset = between(standard_curve_mg_l, 0, 100))$coefficients["standard_curve_mg_l"], 
      med_intercept = lm(rlu ~ standard_curve_mg_l, 
                         subset = between(standard_curve_mg_l, 0, 100))$coefficients["(Intercept)"], 
      med_r2 = summary(lm(rlu ~ standard_curve_mg_l, 
                  subset = between(standard_curve_mg_l, 0, 100)))$r.squared, 
      high_slope = lm(rlu ~ standard_curve_mg_l, 
                     subset = standard_curve_mg_l >= 100)$coefficients["standard_curve_mg_l"], 
      high_intercept = lm(rlu ~ standard_curve_mg_l, 
                         subset = standard_curve_mg_l >= 100)$coefficients["(Intercept)"], 
      high_r2 = summary(lm(rlu ~ standard_curve_mg_l, 
                          subset = standard_curve_mg_l >= 100))$r.squared
      )
      
  # y = mx + c
  # abs = m*ppm + c
  # ppm = abs-c/m
  
  # Calculate nM value for low, medium, and high range of standards
  data_formatted = atp %>% 
    left_join(calibration_coef) %>% 
    mutate(low_ppm_calculate = (rlu - low_intercept) / low_slope) %>% 
      mutate(med_ppm_calculate = (rlu - med_intercept) / med_slope ) %>% 
     mutate(high_ppm_calculate = (rlu - high_intercept) / high_slope ) 
  # mutate(low_ppm_calculate = (rlu*low_slope) + low_intercept ) %>% 
  #   mutate(med_ppm_calculate = (rlu*med_slope) + med_intercept ) %>% 
  #   mutate(high_ppm_calculate = (rlu*high_slope) + high_intercept ) 
  }

## Check QC values from throughout runs
qa = 
  calibrate_atp_data(data_formatted) %>% 
  filter(grepl("standard", sample_id)) %>% 
  filter(is.na(standard_curve_mg_l)) %>% 
  mutate(actual_nM = if_else(rlu <= low_end_range, low_ppm_calculate * dilution_factor,
                     if_else(rlu > low_end_range & rlu <= med_end_range, med_ppm_calculate * dilution_factor, high_ppm_calculate * dilution_factor))) %>% 
  relocate(actual_nM, .after = rlu) %>% 
  relocate(source, .after = actual_nM) %>% 
  drop_na(actual_nM) %>% 
  select(c(sample_id, randomized_id, rlu, actual_nM, source)) %>% 
  group_by(randomized_id) %>% 
  mutate(cv_all = (sd(actual_nM)/mean(actual_nM))*100) %>% 
  ungroup() %>% 
  group_by(source, randomized_id) %>% 
  mutate(cv_day = (sd(actual_nM)/mean(actual_nM))*100)

## Calculate LOD from each day using 0 nM QC values
#range from +/- 0.04
blanks <- qa %>%
  filter(randomized_id == "0_std") %>% 
  group_by(source) %>% 
  mutate(day_lod = mean(actual_nM)) %>% 
  mutate(day_lod_sd = mean(actual_nM) + 3*sd(actual_nM)) %>% 
  mutate(range = max(actual_nM) - min(actual_nM)) %>% 
  ungroup() %>% 
  mutate(all_lod = mean(actual_nM)) 

lod <-  blanks %>% 
  select(c(source,day_lod, day_lod_sd)) %>% 
  distinct(source, .keep_all = TRUE) %>% 
  mutate(day_lod = round(day_lod, 3))

#lod = calibrate_atp_data(data_formatted) %>% 
 # filter(standard_curve_mg_l == 0) %>% 
 # mutate(std_err = sqrt(((0 - as.numeric(n_m_calculation))^2)/(low_count - 2))) %>% 
 # mutate(tray_lod = 3.3*(std_err/low_slope))


## Check reference sediments throughout runs
reference = 
  calibrate_atp_data(data_formatted) %>%
  ungroup() %>% 
  filter(grepl("PE", sample_id)) %>% 
  mutate(dried_sediment_mass_g = as.numeric(dried_sediment_15_m_l_mass_g) - as.numeric(x15_m_l_mass_g)) %>% 
  mutate(edta_est_g = as.numeric(x15_m_l_sediment_edta_wet_mass_g) - as.numeric(x15_m_l_sediment_wet_mass_g_after_decanting_reagent)) %>% 
  mutate(water_mass_g = as.numeric(x15_m_l_sediment_edta_wet_mass_g) - as.numeric(dried_sediment_15_m_l_mass_g) - as.numeric(edta_mass_g)) %>% 
  mutate(water_est_g = as.numeric(x15_m_l_sediment_wet_mass_g_after_decanting_reagent) - as.numeric(dried_sediment_15_m_l_mass_g)) %>% 
  mutate(total_volume_g = as.numeric(water_mass_g) + as.numeric(edta_mass_g)) %>% 
  mutate(total_volume_est = as.numeric(water_est_g) + as.numeric(edta_est_g)) %>%
  mutate(pmol_per_g = (med_ppm_calculate * ((0.0025 + (total_volume_g/1000))/dried_sediment_mass_g))*1000) %>% 
  mutate(cv_all = (sd(pmol_per_g)/mean(pmol_per_g))*100) %>%
  group_by(source) %>% 
  mutate(cv_day = (sd(pmol_per_g)/mean(pmol_per_g))*100) %>% 
  ungroup() %>% 
  select(c(sample_id, dilution_factor, rlu, med_ppm_calculate, pmol_per_g, cv_all, cv_day, source))

mean_ref = mean(reference$pmol_per_g)

# png(file = paste0(processed.data,"ref_standards_med_g.png"), width = 15, height = 15, units = "in", res = 300) 
# 
# ggplot(reference, aes(y = pmol_per_g)) + 
#   geom_boxplot()+
#   facet_grid(~source) + 
#   geom_hline(aes(yintercept = mean))
# 
# dev.off()

## Samples
samples_nm = 
  calibrate_atp_data(data_formatted) %>% 
  ungroup() %>% 
  select(-c(randomized_id, standard_curve_mg_l, n_m_calculation, lab_temperature_deg_c, method_status)) %>% 
  filter(grepl("CM", sample_id)) %>% 
  separate(sample_id, into = c("Sample_Name", "Treat"), remove = FALSE, sep = "_ATP") %>% 
  mutate(actual_nM = if_else(rlu <= low_end_range, low_ppm_calculate * dilution_factor, if_else(rlu > low_end_range & rlu <= med_end_range, med_ppm_calculate * dilution_factor, high_ppm_calculate * dilution_factor))) %>% 
mutate(sc_used = if_else(rlu <= low_end_range, "low", if_else(rlu > low_end_range & rlu <= med_end_range, "med", "high"))) %>% 
  separate(sample_id, c("Sample", "Rep"), sep = -1) %>% 
  group_by(Sample) %>% 
  mutate(cv_samp = (sd(actual_nM)/mean(actual_nM))*100) %>% 
  ungroup() %>% 
  mutate(dried_sed_mass_g = as.numeric(dried_sediment_15_m_l_mass_g) - as.numeric(x15_m_l_mass_g)) %>% 
  mutate(edta_est_g = as.numeric(x15_m_l_sediment_edta_wet_mass_g) - as.numeric(x15_m_l_sediment_wet_mass_g_after_decanting_reagent)) %>% 
  mutate(water_mass_g = (as.numeric(x15_m_l_sediment_edta_wet_mass_g) - as.numeric(edta_mass_g) - as.numeric(x15_m_l_mass_g)) - dried_sed_mass_g) %>% 
  mutate(water_est_g = as.numeric(x15_m_l_sediment_wet_mass_g_after_decanting_reagent) - as.numeric(dried_sediment_15_m_l_mass_g)) %>% 
  mutate(total_volume_g = as.numeric(water_mass_g) + as.numeric(edta_mass_g)) %>% 
  mutate(total_volume_est = water_est_g + edta_est_g) %>% 
  select(c(Sample, Sample_Name, Rep, dilution_factor, rlu, actual_nM, sc_used, source, low_ppm_calculate, med_ppm_calculate, high_ppm_calculate,  dried_sed_mass_g, water_mass_g, water_est_g, edta_mass_g, edta_est_g, total_volume_g, total_volume_est)) 

## Calculate pmol ATP/g sediment using water mass estimated from tube masses (353 out of 558) and water mass estimated from moisture tins (205 out of 558)

samples_pmol <- samples_nm %>% 
  mutate(pmol_g_na = (actual_nM * ((0.005)/dried_sed_mass_g))*1000) %>% 
  mutate(pmol_g_est = (actual_nM * ((0.0025 + (total_volume_est/1000))/dried_sed_mass_g))*1000) %>% 
  mutate(pmol_g_act = (actual_nM * ((0.0025 + (total_volume_g/1000))/dried_sed_mass_g))*1000) %>% 
  mutate(perc_diff = ((pmol_g_est - pmol_g_na)/pmol_g_est)*100) %>% 
  group_by(Sample) %>% 
  mutate(cv_samp_est = (sd(pmol_g_est)/mean(pmol_g_est))*100) %>% 
  ungroup() %>% 
  group_by(Sample_Name) %>% 
  mutate(cv_kit_est = (sd(pmol_g_est)/mean(pmol_g_est))*100)

# png(file = paste0(processed.data,"samples_5_nM_cutoff.png"), width = 15, height = 15, units = "in", res = 300)   
# 
# samples %>% 
#    filter(dilution_factor == 1) %>% 
# ggplot(samples, mapping = aes(y = actual_nM, x = rlu, color = sc_used))+
#   geom_point()+
#   facet_wrap(~source) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
# 
# dev.off()


png(file = paste0(processed.data,"samples_nM.png"), width = 15, height = 15, units = "in", res = 300) 

ggplot(samples_nm, aes(y = actual_nM, x = Sample)) + 
  geom_boxplot()

dev.off()

png(file = paste0(processed.data,"samples_pmol_g.png"), width = 15, height = 15, units = "in", res = 300) 

ggplot(samples_pmol, aes(y = pmol_g_est, x = Sample)) + 
  geom_boxplot()

dev.off()

#flag samples below LOD

samples_lod <- left_join(samples_pmol,lod, by = "source") %>%  
  mutate(actual_nM = round(actual_nM, 2)) %>% 
  mutate(pmol_g_est = round(pmol_g_est, 2)) %>% 
  mutate(actual_nM_lod = if_else(actual_nM <= day_lod, paste0("ATP_Below_",day_lod,"_nanomol_per_L_LOD|", actual_nM, "_nanomol_per_L_Raw_Not_Corrected|", actual_nM,"_nanomol_per_L_Final_Corrected"), as.character(actual_nM))) %>% 
mutate(pmol_lod = if_else(actual_nM <= day_lod, paste0("ATP_Below_",day_lod,"_nanomol_per_L_LOD|", actual_nM, "_nanomol_per_L_Raw_Not_Corrected|", pmol_g_est,"_picomol_per_g_Final_Corrected"), as.character(pmol_g_est)))

samples_final <- samples_lod %>% 
  unite(Sample_Name, Sample, Rep, sep = "") %>% 
  rename(ATP_nanomol_per_L = actual_nM_lod) %>% 
  rename(ATP_picomol_per_g = pmol_lod) %>% 
  mutate(ATP_nanomol_per_L = if_else(rlu == -9999, "-9999", ATP_nanomol_per_L)) %>% 
  mutate(ATP_nanomol_per_L = if_else(is.na(ATP_nanomol_per_L), "-9999", ATP_nanomol_per_L)) %>% 
  mutate(ATP_picomol_per_g = if_else(rlu == -9999, "-9999", ATP_picomol_per_g)) %>% 
  mutate(ATP_picomol_per_g = if_else(is.na(ATP_picomol_per_g), "-9999", ATP_picomol_per_g)) %>%
  mutate(Methods_Deviation = if_else(dilution_factor == 2, "ATP_001", "N/A")) %>%
  dplyr::select(c(Sample_Name, ATP_nanomol_per_L, ATP_picomol_per_g, Methods_Deviation))


### ATP Summary File ####

cv <- function(x) {
  cv_value <- (sd(x) / mean(x)) * 100
  return(cv_value)
}

missing = samples_final %>% 
  filter(ATP_nanomol_per_L == -9999) %>% 
  separate(Sample_Name, c("Sample_Name", "Rep"), sep = -1) %>% 
  group_by(Sample_Name) %>% 
  summarise(Sample_Name, count = n())

atp_outliers = samples_final %>% 
  separate(Sample_Name, c("Sample_ID", "Rep"), sep = -1) %>% 
  filter(ATP_nanomol_per_L != -9999) %>% 
  select(c(Sample_ID, Rep, ATP_nanomol_per_L, ATP_picomol_per_g)) %>%   mutate(ATP_nanomol_per_L = as.numeric(ATP_nanomol_per_L)) %>% 
  mutate(ATP_picomol_per_g = as.numeric(ATP_picomol_per_g)) %>% 
  group_by(Sample_ID) %>% 
  mutate(ATP_picomol_per_g_cv = cv(ATP_picomol_per_g))

atp_removed_outliers = as.data.frame(matrix(NA, ncol = 7, nrow = 1))

colnames(atp_removed_outliers) = c("Sample_ID", "Rep", "ATP_nanomol_per_L", "ATP_picomol_per_g", "ATP_picomol_per_g_cv", "CV_after", "flag")

unique.samples = unique(atp_outliers$Sample_ID)

for (i in 1:length(unique.samples)) {
  
  data_subset = subset(atp_outliers, atp_outliers$Sample_ID == unique.samples[i])
  
  conc.temp = as.numeric(data_subset$ATP_picomol_per_g)
  
  conc.temp.sd <- sd(conc.temp)
  conc.temp.mean <- mean(conc.temp)
  CV = (conc.temp.sd/conc.temp.mean)*100
  
  #looping to get 2 out of 3 best samples
  for (sample.reduction in 1:3)  {
    
    if (conc.temp.mean == 0) {
      
      CV = 0
      
    }
    
    else if (length(conc.temp) > 2 & CV >= 30) {
      
      dist.temp = as.matrix(abs(dist(conc.temp)))
      dist.comp = numeric()
      
      for(conc.now in 1:ncol(dist.temp)) {
        
        dist.comp = rbind(dist.comp,c(conc.now,sum(dist.temp[,conc.now])))
        
      }
      
      dist.comp[,2] = as.numeric(dist.comp[,2])
      conc.temp = conc.temp[-which.max(dist.comp[,2])]
      
      conc.temp.sd <- sd(conc.temp)
      conc.temp.mean <- mean(conc.temp)
      conc.temp.cv <- (conc.temp.sd/conc.temp.mean)*100
      CV = conc.temp.cv
      
    }
    
  }
  
  if(length(conc.temp) >= 2) {
    
    if(CV > 30){
      
      samples_removed_final$Concentration[which(samples_removed_final$Sample_ID == unique.samples[i])] = "Samples too Variable"
      
    }
    
    else {
      
      samples_removed_combined = as.data.frame(conc.temp)
      
      samples_remove = merge(samples_removed_combined, data_subset, by.x = "conc.temp", by.y = "ATP_picomol_per_g", all.x = TRUE)
      
      #samples_remove = samples_remove[!duplicated(samples_remove$Sample_ID), ]
      
      samples_remove_omit = merge(samples_remove, data_subset, by = c("Sample_ID", "Rep"), all = TRUE)
      
      samples_remove_omit$CV_after = as.numeric(abs((sd(conc.temp)/mean(conc.temp))*100))
      
      samples_remove_omit$flag[is.na(samples_remove_omit$conc.temp)] = "OMIT"
      
      samples_remove_omit <- samples_remove_omit %>% 
        dplyr::select(-c(ATP_picomol_per_g_cv.x, ATP_nanomol_per_L.x, conc.temp)) %>% 
        rename(ATP_picomol_per_g_cv = ATP_picomol_per_g_cv.y) %>% 
        rename(ATP_nanomol_per_L = ATP_nanomol_per_L.y)
      
    }
    
  }
  
  atp_removed_outliers = rbind(samples_remove_omit, atp_removed_outliers)
  
  rm('conc.temp')
}

samples_final_outliers = atp_removed_outliers %>% 
  unite(Sample_Name, c("Sample_ID", "Rep"), sep = "", remove = FALSE) %>% 
  left_join(samples_final, by = c("Sample_Name")) %>% 
  select(Sample_Name, Sample_ID, ATP_nanomol_per_L.y, ATP_picomol_per_g.y, flag, Methods_Deviation) %>% 
  rename(ATP_nanomol_per_L = ATP_nanomol_per_L.y) %>% 
  rename(ATP_picomol_per_g = ATP_picomol_per_g.y) %>% 
  group_by(Sample_ID) %>% 
  mutate(Methods_Deviation.flag = ifelse(flag == "OMIT", "ATP_CV_030", "N/A")) %>% 
  fill(Methods_Deviation.flag, .direction = 'downup') %>% 
  mutate(Methods_Deviation_new = ifelse(is.na(Methods_Deviation.flag), Methods_Deviation, ifelse(Methods_Deviation != "N/A", paste(Methods_Deviation, Methods_Deviation.flag, sep = ";"), Methods_Deviation.flag))) %>% 
  ungroup() %>% 
  drop_na(Sample_ID) %>% 
  select(c(Sample_Name, ATP_nanomol_per_L, ATP_picomol_per_g, Methods_Deviation_new, flag)) %>% 
           rename(Methods_Deviation = Methods_Deviation_new)

write.csv(samples_final_outliers, paste0(processed.data,"/CM_ATP_ReadyForBoye_04-12-2024.csv"), row.names = F) 

atp_summary = atp_removed_outliers %>% 
  separate(Sample_Name, c("Sample_ID", "Rep"), sep = -1) %>% 
  filter(ATP_nanomol_per_L != -9999) %>% 
  select(c(Sample_ID, ATP_nanomol_per_L, ATP_picomol_per_g)) %>% 
  mutate(ATP_nanomol_per_L = as.numeric(ATP_nanomol_per_L)) %>% 
  mutate(ATP_picomol_per_g = as.numeric(ATP_picomol_per_g)) %>% 
  group_by(Sample_ID) %>% 
  summarise(across(where(is.numeric), list(mean = mean, sd = sd, cv = cv))) %>% 
  ungroup() %>% 
  mutate(round(across(where(is.numeric)),3)) %>% 
  rename(Mean_ATP_nanomol_per_L = ATP_nanomol_per_L_mean) %>% 
  rename(Mean_ATP_picomol_per_g = ATP_picomol_per_g_mean) %>% 
  rename(SD_ATP_nanomol_per_L = ATP_nanomol_per_L_sd) %>% 
  rename(SD_ATP_picomol_per_g = ATP_picomol_per_g_sd) %>% 
  rename(Sample_Name = Sample_ID)

ggplot(atp_summary, aes(x = ATP_picomol_per_g_cv)) +
  geom_histogram()+
  geom_vline(xintercept = 30)

atp_summary_file = atp_summary %>% 
  mutate(Material = "Sediment") %>% 
  select(c(Sample_Name, Material, Mean_ATP_nanomol_per_L, SD_ATP_nanomol_per_L, Mean_ATP_picomol_per_g, SD_ATP_picomol_per_g)) %>% 
  left_join(missing, by = "Sample_Name")

write.csv(atp_summary_file, paste0(processed.data,"/EC_ATP_Summary_ReadyForBoye_04-10-2024.csv"), row.names = F) 
