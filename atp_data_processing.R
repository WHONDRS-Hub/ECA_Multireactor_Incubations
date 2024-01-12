library(tidyverse)
library(dplyr)
library(readxl)
library(ggpmisc)

rm(list=ls());graphics.off()

pnnl.user = 'laan208'

raw.data = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/ATP/01_RawData/")

processed.data = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/ATP/03_ProcessedData/")

moisture <- read_csv(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/MOI/03_ProcessedData/ECA_Moisture_Outliers_Removed.csv"))

moisture_clean <- moisture %>%
  filter(is.na(flag)) %>% 
  separate(Sample_Name, c("EC", "Site", "Rep"), sep = "_", remove = FALSE) %>% 
  unite(Sample_Name, c("EC", "Site")) %>% 
  group_by(Sample_Name) %>% 
  mutate(average_grav = mean(percent_water_content_dry)) %>% 
  mutate(cv = (sd(percent_water_content_dry)/mean(percent_water_content_dry))*100)

moisture_average <- moisture_clean %>% 
  distinct(Sample_Name, .keep_all = TRUE) %>% 
  select(c(Sample_Name, average_grav))

#setwd(input.path)

filePaths_atp<- list.files(path = raw.data, pattern = "Data_Raw_ATP", full.names = TRUE, recursive = TRUE) 

atp <- do.call(bind_rows, lapply(filePaths_atp, function(raw.data) {
  df <- read_xlsx(raw.data, skip = 1)%>% mutate_all(as.character) %>% janitor::clean_names()
  df = df %>% mutate(source = basename(raw.data))
  df}))

atp$source <- substr(atp$source, start = 1, stop = 8)

atp <- atp %>% 
  mutate(rlu = as.numeric(rlu)) %>% 
  mutate(standard_curve_mg_l = as.numeric(standard_curve_mg_l)) %>% 
  mutate(n_m_calculation = as.numeric(n_m_calculation)) %>% 
  group_by(source) %>% 
  mutate(low_beg_range = rlu[standard_curve_mg_l == 0.00][1]) %>% 
  mutate(low_end_range = rlu[standard_curve_mg_l == 10][1]) %>% 
  mutate(med_beg_range = rlu[standard_curve_mg_l == 10][1]) %>% 
  mutate(med_end_range = rlu[standard_curve_mg_l == 100][1]) %>% 
  mutate(high_beg_range = rlu[standard_curve_mg_l == 100][1]) %>% 
  mutate(high_end_range = rlu[standard_curve_mg_l == 500][1]) %>% 
  mutate(dilution_factor = as.numeric(dilution_factor))


  
## choose which standard curve to use for Fe ####

calibrate_atp_data = function(data_formatted){
  
  standards = 
    atp %>% 
    filter(grepl("standard", sample_id)) %>% 
    dplyr::select(randomized_id, rlu, standard_curve_mg_l, source) %>% 
    drop_na(standard_curve_mg_l) %>% 
    drop_na(rlu)
   
   
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
  #   filter(standard_curve_mg_l <= 10) %>% 
  #   ggplot(aes(y = standard_curve_mg_l, x = rlu, color = as.character(source)))+
  #   geom_point()+
  #   #stat_poly_line(method = 'lm', 
  #                  #formula = y ~ poly(x,2))+
  #   stat_poly_line()+
  #   stat_poly_eq(#method = 'lm', 
  #                #formula = y ~ poly(x,2),
  #                use_label(c("eq")), label.y = 0.9)+
  #   stat_poly_eq(#method = 'lm', 
  #                #formula = y ~ poly(x,2),
  #                label.y = 0.85, rr.digits = 3)+
  #   facet_wrap(~source) + 
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  # 
  # dev.off()
  # 
  # png(file = paste0(processed.data,"med_standards_linear.png"), width = 15, height = 15, units = "in", res = 300)   
  # 
  # standards %>% 
  #   filter(between(standard_curve_mg_l, 10, 100)) %>% 
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
  # 
  # png(file = paste0(processed.data,"low_med_standards_linear.png"), width = 15, height = 15, units = "in", res = 300)   
  # 
  # standards %>% 
  #   filter(standard_curve_mg_l <= 100) %>% 
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
  
  
  calibration_coef = 
    standards %>% 
    dplyr::group_by(source) %>% 
    dplyr::summarize(
      all_poly_x = lm(standard_curve_mg_l~ poly(rlu,2))$coefficients["poly(rlu, 2)1"],
      all_poly_x2 = lm(standard_curve_mg_l~ poly(rlu,2))$coefficients["poly(rlu, 2)2"],
                          all_poly_intercept = lm(standard_curve_mg_l ~ poly(rlu,2))$coefficients["(Intercept)"], 
                          poly_r2 = summary(lm(standard_curve_mg_l ~ poly(rlu,2)))$r.squared, 
      #stat_poly_line(method = 'lm', 
                    # formula = y ~ poly(x,2))
      all_slope = lm(standard_curve_mg_l~rlu)$coefficients["rlu"], 
                     all_intercept = lm(standard_curve_mg_l~rlu)$coefficients["(Intercept)"], 
                     all_r2 = summary(lm(standard_curve_mg_l~rlu))$r.squared,
      low_slope = lm(standard_curve_mg_l~rlu, 
                     subset = standard_curve_mg_l <= 10)$coefficients["rlu"], 
      low_intercept = lm(standard_curve_mg_l~rlu, 
                     subset = standard_curve_mg_l <= 10)$coefficients["(Intercept)"], 
      low_r2 = summary(lm(standard_curve_mg_l~rlu, 
                  subset = standard_curve_mg_l <= 10))$r.squared, 
      med_slope = lm(standard_curve_mg_l~rlu, 
                     subset = between(standard_curve_mg_l, 10, 100))$coefficients["rlu"], 
      med_intercept = lm(standard_curve_mg_l~rlu, 
                         subset = between(standard_curve_mg_l, 10, 100))$coefficients["(Intercept)"], 
      med_r2 = summary(lm(standard_curve_mg_l~rlu, 
                  subset = between(standard_curve_mg_l, 10, 100)))$r.squared, 
      high_slope = lm(standard_curve_mg_l~rlu, 
                     subset = standard_curve_mg_l >= 100)$coefficients["rlu"], 
      high_intercept = lm(standard_curve_mg_l~rlu, 
                         subset = standard_curve_mg_l >= 100)$coefficients["(Intercept)"], 
      high_r2 = summary(lm(standard_curve_mg_l~rlu, 
                          subset = standard_curve_mg_l >= 100))$r.squared,
      low_med_slope = lm(standard_curve_mg_l~rlu, 
                     subset = standard_curve_mg_l <= 100)$coefficients["rlu"], 
      low_med_intercept = lm(standard_curve_mg_l~rlu, 
                         subset = standard_curve_mg_l <= 100)$coefficients["(Intercept)"], 
      low_med_r2 = summary(lm(standard_curve_mg_l~rlu, 
                          subset = standard_curve_mg_l <= 100))$r.squared,
      )
      
      
  # y = mx + c
  # abs = m*ppm + c
  # ppm = abs-c/m
  
  data_formatted = atp %>% 
    left_join(calibration_coef) %>% 
    mutate(all_ppm_calculated = (rlu * all_slope) + all_intercept) %>% 
    mutate(poly_ppm_calculate = all_poly_intercept + all_poly_x*rlu + all_poly_x2*(rlu^2)) %>% 
  mutate(low_ppm_calculate = (rlu*low_slope) + low_intercept ) %>% 
    mutate(med_ppm_calculate = (rlu*med_slope) + med_intercept ) %>% 
    mutate(high_ppm_calculate = (rlu*high_slope) + high_intercept ) %>% 
    mutate(low_med_ppm_calculate = (rlu*low_med_slope) + low_med_intercept) 
  }

qa = 
  calibrate_atp_data(data_formatted) %>% 
  filter(grepl("standard", sample_id)) %>% 
  filter(is.na(standard_curve_mg_l)) %>% 
  mutate(actual_nM = if_else(rlu <= low_end_range,                              low_ppm_calculate * dilution_factor,
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

reference = 
  calibrate_atp_data(data_formatted) %>% 
  filter(grepl("PE", sample_id)) %>% 
  mutate(dried_sediment_mass_g = as.numeric(dried_sediment_15_m_l_mass_g) - as.numeric(x15_m_l_mass_g)) %>% 
  mutate(water_mass_g = as.numeric(x15_m_l_sediment_edta_wet_mass_g) - as.numeric(x15_m_l_mass_g) - as.numeric(edta_mass_g)) %>% 
  mutate(pmol_per_g_na = (med_ppm_calculate * (0.005/dried_sediment_mass_g))*1000) %>%
  mutate(pmol_per_g_corr = (med_ppm_calculate * ((0.005 + (water_mass_g/1000))/dried_sediment_mass_g))*1000) %>% 
  ungroup() %>% 
  mutate(cv_all_na = (sd(pmol_per_g_na)/mean(pmol_per_g_na))*100) %>% 
  mutate(cv_all_corr = (sd(pmol_per_g_corr)/mean(pmol_per_g_corr))*100) %>%
  group_by(source) %>% 
  mutate(cv_day_na = (sd(pmol_per_g_na)/mean(pmol_per_g_na))*100) %>% 
  mutate(cv_day_corr = (sd(pmol_per_g_corr)/mean(pmol_per_g_corr))*100) %>% 
  ungroup() %>% 
  select(c(sample_id, dilution_factor, rlu, med_ppm_calculate, pmol_per_g_na, pmol_per_g_corr, cv_all_na, cv_all_corr, cv_day_na, cv_day_corr, source))

mean = mean(reference$pmol_per_g_corr)

png(file = paste0(processed.data,"ref_standards_med_g.png"), width = 15, height = 15, units = "in", res = 300) 

ggplot(reference, aes(y = pmol_per_g_corr)) + 
  geom_boxplot()+
  facet_grid(~source) + 
  geom_hline(aes(yintercept = 412))

dev.off()

ggplot(reference, aes(x = cv_day)) + 
  geom_histogram(binwidth = 2) 


## Samples
samples = 
  calibrate_atp_data(data_formatted) %>% 
  ungroup() %>% 
  select(-c(randomized_id, standard_curve_mg_l)) %>% 
  filter(grepl("EC", sample_id)) %>% 
  separate(sample_id, into = c("Sample_Name", "Treat"), remove = FALSE, sep = "_ATP") %>% 
  left_join(moisture_average, by = "Sample_Name") %>% 
  #mutate(dried_sediment_mass_g = as.numeric(dried_sediment_15_m_l_mass_g) - as.numeric(x15_m_l_mass_g)) %>% 
 # mutate(nmol_per_g = med_ppm_calculate * (0.005/dried_sediment_mass_g)) %>% 
  mutate(actual_nM = if_else(rlu <= low_end_range,                              low_ppm_calculate * dilution_factor,
                     if_else(rlu > low_end_range & rlu <= med_end_range, med_ppm_calculate * dilution_factor, high_ppm_calculate * dilution_factor))) %>% 
  mutate(low_med_ppm_calculate = low_med_ppm_calculate * dilution_factor) %>% 
  mutate(sc_used = if_else(rlu <= low_end_range,                              "low",
                           if_else(rlu > low_end_range & rlu <= med_end_range, "med", "high"))) %>% 
  mutate(sc_used = if_else(actual_nM <= 10 & rlu >= med_beg_range, "low_med", sc_used)) %>% 
  mutate(actual_nM = if_else(actual_nM <= 10 & rlu >= med_beg_range, low_med_ppm_calculate * dilution_factor, actual_nM)) %>% 
  separate(sample_id, c("Sample", "Rep"), sep = -1) %>% 
  group_by(Sample) %>% 
  mutate(cv_samp = (sd(actual_nM)/mean(actual_nM))*100) %>% 
  mutate(dried_sed_mass_g = as.numeric(dried_sediment_15_m_l_mass_g) - as.numeric(x15_m_l_mass_g)) %>% 
  mutate(water_mass_g = (as.numeric(x15_m_l_sediment_edta_wet_mass_g) - as.numeric(edta_mass_g) - as.numeric(x15_m_l_mass_g)) - dried_sed_mass_g) %>% 
  mutate(water_mass_moi_g = dried_sed_mass_g * (average_grav/100)) %>% 
  mutate(pmol_g_na = (actual_nM * (0.005/dried_sed_mass_g))*1000) %>% 
  mutate(pmol_g_corr = (actual_nM * ((0.005 + (water_mass_g/1000))/dried_sed_mass_g))*1000) %>% 
  mutate(pmol_g_moi = (actual_nM * ((0.005 + (water_mass_moi_g/1000))/dried_sed_mass_g))*1000) %>% 
  select(c(Sample, Sample_Name, Rep, dilution_factor, rlu, actual_nM, sc_used, source, low_ppm_calculate, med_ppm_calculate, low_med_ppm_calculate, high_ppm_calculate, average_grav, dried_sed_mass_g, water_mass_g, water_mass_moi_g, pmol_g_na, pmol_g_corr, pmol_g_moi))  %>% 
  mutate(error_low_med = ((low_med_ppm_calculate - actual_nM)/actual_nM)*100) %>% 
   relocate(error_low_med, .after = actual_nM) %>% 
  group_by(Sample) %>% 
  mutate(cv_samp_moi = (sd(pmol_g_moi)/mean(pmol_g_moi))*100) %>% 
  mutate(cv_samp_corr = (sd(pmol_g_corr)/mean(pmol_g_corr))*100) %>% 
  mutate(cv_samp_na = (sd(pmol_g_na)/mean(pmol_g_na))*100) %>% 
  ungroup() %>% 
  group_by(Sample_Name) %>% 
  mutate(cv_kit_moi = (sd(pmol_g_moi)/mean(pmol_g_moi))*100) %>% 
  mutate(cv_kit_corr = (sd(pmol_g_corr)/mean(pmol_g_corr))*100) %>% 
  mutate(cv_kit_na = (sd(pmol_g_na)/mean(pmol_g_na))*100)


png(file = paste0(processed.data,"samples_nM.png"), width = 15, height = 15, units = "in", res = 300) 

ggplot(samples, aes(y = actual_nM, x = Sample)) + 
  geom_boxplot()

dev.off()

png(file = paste0(processed.data,"samples_pmol_g_corr.png"), width = 15, height = 15, units = "in", res = 300) 

ggplot(samples, aes(y = pmol_g_corr, x = Sample)) + 
  geom_boxplot()

dev.off()

ggplot(samples, aes(x = cv_samp)) + 
  geom_histogram()

ggplot(samples, aes(x = rlu, y = actual_nM)) + 
  geom_point(aes(color = source))
  