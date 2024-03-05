library(tidyverse)
library(readxl)
library(ggpmisc)

rm(list=ls());graphics.off()

pnnl.user = 'laan208'

raw.data = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/ATP/01_RawData/")

processed.data = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/ATP/03_ProcessedData/")


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
      low_slope = lm(standard_curve_mg_l~rlu, 
                     subset = standard_curve_mg_l <= 5)$coefficients["rlu"], 
      low_intercept = lm(standard_curve_mg_l~rlu, 
                     subset = standard_curve_mg_l <= 5)$coefficients["(Intercept)"], 
      low_r2 = summary(lm(standard_curve_mg_l~rlu, 
                  subset = standard_curve_mg_l <= 5))$r.squared, 
      med_slope = lm(standard_curve_mg_l~rlu, 
                     subset = between(standard_curve_mg_l, 0, 100))$coefficients["rlu"], 
      med_intercept = lm(standard_curve_mg_l~rlu, 
                         subset = between(standard_curve_mg_l, 0, 100))$coefficients["(Intercept)"], 
      med_r2 = summary(lm(standard_curve_mg_l~rlu, 
                  subset = between(standard_curve_mg_l, 0, 100)))$r.squared, 
      high_slope = lm(standard_curve_mg_l~rlu, 
                     subset = standard_curve_mg_l >= 100)$coefficients["rlu"], 
      high_intercept = lm(standard_curve_mg_l~rlu, 
                         subset = standard_curve_mg_l >= 100)$coefficients["(Intercept)"], 
      high_r2 = summary(lm(standard_curve_mg_l~rlu, 
                          subset = standard_curve_mg_l >= 100))$r.squared
      )
      
  # y = mx + c
  # abs = m*ppm + c
  # ppm = abs-c/m
  
  # Calculate nM value for low, medium, and high range of standards
  data_formatted = atp %>% 
    left_join(calibration_coef) %>% 
  mutate(low_ppm_calculate = (rlu*low_slope) + low_intercept ) %>% 
    mutate(med_ppm_calculate = (rlu*med_slope) + med_intercept ) %>% 
    mutate(high_ppm_calculate = (rlu*high_slope) + high_intercept ) 
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

mean = mean(reference$pmol_per_g)

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
  filter(grepl("EC", sample_id)) %>% 
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


write.csv(samples_final, paste0(processed.data,"/EC_ATP_ReadyForBoye_01-26-2024.csv"), row.names = F) 

### ATP Summary File ####

cv <- function(x) {
  cv_value <- (sd(x) / mean(x)) * 100
  return(cv_value)
}

atp_summary = samples_final %>% 
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
  geom_histogram()

atp_summary_file = atp_summary %>% 
  mutate(Material = "Sediment") %>% 
  select(c(Sample_Name, Material, Mean_ATP_nanomol_per_L, SD_ATP_nanomol_per_L, Mean_ATP_picomol_per_g, SD_ATP_picomol_per_g)) 

write.csv(atp_summary_file, paste0(processed.data,"/EC_ATP_Summary_ReadyForBoye_03-05-2024.csv"), row.names = F) 
