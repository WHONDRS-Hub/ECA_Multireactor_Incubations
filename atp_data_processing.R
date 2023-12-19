library(tidyverse)
library(dplyr)
library(readxl)
library(ggpmisc)

pnnl.user = 'laan208'

raw.data = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/ATP/01_RawData/")

processed.data = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/ATP/03_ProcessedData/")

#setwd(input.path)

filePaths_atp<- list.files(path = raw.data, pattern = "Data_Raw_ATP", full.names = TRUE, recursive = TRUE) 

atp <- do.call(bind_rows, lapply(filePaths_atp, function(raw.data) {
  df <- read_xlsx(raw.data, skip = 1)%>% mutate_all(as.character) %>% janitor::clean_names()
  df = df %>% mutate(source = basename(raw.data))
  df}))

atp$source <- substr(atp$source, start = 1, stop = 8)
  
## choose which standard curve to use for Fe ####

calibrate_atp_data = function(data_formatted){
  
  standards = 
    atp %>% 
    filter(grepl("standard", sample_id)) %>% 
    dplyr::select(randomized_id, rlu, standard_curve_mg_l, source) %>% 
    mutate(rlu = as.numeric(rlu)) %>% 
    mutate(standard_curve_mg_l = as.numeric(standard_curve_mg_l)) %>% 
    drop_na(standard_curve_mg_l) %>% 
    drop_na(rlu)
  

  png(file = paste0(processed.data,"all_standards_linear.png"), width = 15, height = 15, units = "in", res = 300) 
  
  standards %>% 
    ggplot(aes(y = standard_curve_mg_l, x = rlu, color = as.character(source)))+
    geom_point()+
    stat_poly_line()+
    stat_poly_eq(use_label(c("eq")), label.y = 0.9)+
    stat_poly_eq(label.y = 0.85)+
    #geom_smooth(method = "lm", se = F)+
    facet_wrap(~source) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  
  dev.off()
  
  png(file = paste0(processed.data,"all_standards_nonlinear.png"), width = 15, height = 15, units = "in", res = 300) 
  
  standards %>% 
    ggplot(aes(y = standard_curve_mg_l, x = rlu, color = as.character(source)))+
    geom_point()+
    stat_poly_line(method = 'lm', 
                formula = y ~ poly(x,2))+
    #stat_poly_line()+
    stat_poly_eq(method = 'lm', 
                 formula = y ~ poly(x,2),
                 use_label(c("eq")), label.y = 0.9)+
    stat_poly_eq(method = 'lm', 
                 formula = y ~ poly(x,2),
                 label.y = 0.85)+
    facet_wrap(~source) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  
  dev.off()

  png(file = paste0(processed.data,"low_standards_linear.png"), width = 15, height = 15, units = "in", res = 300)   

  standards %>% 
    filter(standard_curve_mg_l <= 10) %>% 
    ggplot(aes(y = standard_curve_mg_l, x = rlu, color = as.character(source)))+
    geom_point()+
    #stat_poly_line(method = 'lm', 
                   #formula = y ~ poly(x,2))+
    stat_poly_line()+
    stat_poly_eq(#method = 'lm', 
                 #formula = y ~ poly(x,2),
                 use_label(c("eq")), label.y = 0.9)+
    stat_poly_eq(#method = 'lm', 
                 #formula = y ~ poly(x,2),
                 label.y = 0.85, rr.digits = 3)+
    facet_wrap(~source) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  
  dev.off()
  
  png(file = paste0(processed.data,"med_standards_linear.png"), width = 15, height = 15, units = "in", res = 300)   
  
  standards %>% 
    filter(between(standard_curve_mg_l, 10, 100)) %>% 
    ggplot(aes(y = standard_curve_mg_l, x = rlu, color = as.character(source)))+
    geom_point()+
    #stat_poly_line(method = 'lm', 
    #formula = y ~ poly(x,2))+
    stat_poly_line()+
    stat_poly_eq(#method = 'lm', 
      #formula = y ~ poly(x,2),
      use_label(c("eq")), label.y = 0.9)+
    stat_poly_eq(#method = 'lm', 
      #formula = y ~ poly(x,2),
      label.y = 0.85, rr.digits = 3)+
    facet_wrap(~source) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  
  dev.off()
  
  png(file = paste0(processed.data,"high_standards_linear.png"), width = 15, height = 15, units = "in", res = 300)   
  
  standards %>% 
    filter(standard_curve_mg_l >= 100) %>% 
    ggplot(aes(y = standard_curve_mg_l, x = rlu, color = as.character(source)))+
    geom_point()+
    #stat_poly_line(method = 'lm', 
    #formula = y ~ poly(x,2))+
    stat_poly_line()+
    stat_poly_eq(#method = 'lm', 
      #formula = y ~ poly(x,2),
      use_label(c("eq")), label.y = 0.9)+
    stat_poly_eq(#method = 'lm', 
      #formula = y ~ poly(x,2),
      label.y = 0.85, rr.digits = 3)+
    facet_wrap(~source) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  
  dev.off()
  
  
  calibration_coef = 
    standards %>% 
    dplyr::group_by(tray_number) %>% 
    dplyr::summarize(slope = lm(absorbance_562 ~ standard_ppm)$coefficients["standard_ppm"], 
                     intercept = lm(absorbance_562 ~ standard_ppm)$coefficients["(Intercept)"], 
                     r2 = summary(lm(absorbance_562 ~ standard_ppm))$r.squared)
  
  # y = mx + c
  # abs = m*ppm + c
  # ppm = abs-c/m
  
  data_formatted = data_formatted %>% 
    left_join(calibration_coef) %>% 
    mutate(ppm_calculated = ((absorbance_562 - intercept) / slope))
  
}
