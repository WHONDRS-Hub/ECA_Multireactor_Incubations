library(tidyverse);library(dplyr);library(readxl)

rm(list=ls());graphics.off()

pnnl.user = 'laan208'

#### Read in Data ####
input.path = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/FE/")

setwd(input.path)

raw.data = ("01_Rawdata//20230519_Data_Raw_SFE_ECA_EC/")

formatted.data = ("02_FormattedData/20230110_Data_Formatted_SFE_ECA_EC_1-270/")

processed.data = ("03_ProcessedData/20230110_Data_Formatted_SFE_ECA_EC_1-270/")


  
  # import map
  ferrozine_map = read_excel("01_Rawdata/20230519_Data_Raw_SFE_ECA_EC/20230519_Mapping_SFE_ECA_EC.xlsx", sheet = "map") %>% mutate_all(as.character) 
  # import data files (plate reader)
  filePaths_ferrozine <- list.files(path = raw.data, pattern = "Tray", full.names = TRUE, recursive = TRUE)
  ferrozine_data <- do.call(bind_rows, lapply(filePaths_ferrozine, function(raw.data) {
    df <- read_xlsx(raw.data, skip = 24) %>% mutate_all(as.character) %>% janitor::clean_names()
    df = df %>% mutate(source = basename(raw.data))
    df}))
  

## clean the map ####
map_processed = 
  ferrozine_map %>% 
  mutate(tray_number = parse_number(tray_number)) %>% 
  filter(!is.na(sample_name) & !is.na(tray_number)) %>% 
  rename(sample_label = sample_name)


## remove samples that had to be diluted or plates that needed to be rerun ####
data_formatted = 
  ferrozine_data %>% 
  mutate_all(na_if,"") %>% 
  #rename("x" = "x1") %>% 
  #rename("x1" = "x1_2") %>% 
  dplyr::select(-x14) %>% 
  #fill(x_1) %>% 
  #filter(x_2 == "562") %>% 
  #dplyr::select(-x_2) %>% 
  pivot_longer(-c(source, x), values_to = "absorbance_562") %>% 
  mutate(name = str_remove(name, "x"),
         well_position = paste0(x, name),
         tray_number = str_extract(source, "Tray[0-9]+"),
         tray_number = parse_number(tray_number),
         absorbance_562 = as.numeric(absorbance_562)) %>% 
  dplyr::select(tray_number, well_position, absorbance_562) %>% 
  right_join(map_processed, by = c("tray_number", "well_position")) %>% 
    filter(!Method_Deviations %in% "OMIT"
 

#### Formatted absorbance data ####
write.csv(data_formatted, paste0(formatted.data,"20230110_Data_Formatted_SFE_ECA_EC_1-270.csv"), row.names = F)
  

####Processed Data #####

## choose which standard curve to use for Fe ####

calibrate_ferrozine_data = function(data_formatted){
  standards = 
    data_formatted %>% 
    filter(grepl("standard", sample_label)) %>% 
    dplyr::select(tray_number, absorbance_562, standard_ppm) %>% 
    mutate(standard_ppm = as.numeric(standard_ppm))
  
  standards %>% 
    ggplot(aes(x = standard_ppm, y = absorbance_562, color = as.character(tray_number)))+
    geom_point()+
    geom_smooth(method = "lm", se = F)+
    facet_wrap(~tray_number)
  
  calibration_coef = 
    standards %>% 
    dplyr::group_by(tray_number) %>% 
    dplyr::summarize(slope = lm(absorbance_562 ~ standard_ppm)$coefficients["standard_ppm"], 
                     intercept = lm(absorbance_562 ~ standard_ppm)$coefficients["(Intercept)"])
  
  # y = mx + c
  # abs = m*ppm + c
  # ppm = abs-c/m
  
  
 data_formatted = data_formatted %>% 
    left_join(calibration_coef) %>% 
    mutate(ppm_calculated = ((absorbance_562 - intercept) / slope))
  
}
## Check Standards CV ####
  standards = calibrate_ferrozine_data(data_formatted) %>% 
    filter(sample_label == "FAS-standard") %>% 
    group_by(standard_ppm) %>% 
    mutate(range = max(ppm_calculated) - min(ppm_calculated)) %>% 
    mutate(CV = (sd(ppm_calculated)/mean(ppm_calculated))*100)
  
## Check Blanks ####
  blanks = 
    calibrate_ferrozine_data(data_formatted) %>% 
    filter(grepl("blank", sample_label)) %>% 
    mutate(CV = ((sd(ppm_calculated)/mean(ppm_calculated))*100)) %>% 
    mutate(range = max(ppm_calculated) - min(ppm_calculated))

## Check Samples ####
samples = 
  calibrate_ferrozine_data(data_formatted) %>% 
  filter(grepl("EC", sample_label)) %>% 
  dplyr::select(sample_label, analysis, ppm_calculated) %>% 
  mutate(ppm_calculated = case_when(
  analysis  == "dilute" ~ ppm_calculated * 2, 
  analysis == "Fe2" ~ ppm_calculated)) %>% 
  rename("mg_Fe_per_L" = "ppm_calculated") %>% 
  mutate(mg_Fe_per_L = if_else(mg_Fe_per_L<0,0,mg_Fe_per_L))

## Flag samples with range > 0.04 and CV > 10% ####
data_flag_conc <- samples %>%
  separate(col = sample_label, into = c("Project", "kit", "analysis"), sep = "_", remove = FALSE) %>% 
  separate(col = analysis, into = c("Analysis", "Replicate"), sep = "-", remove = FALSE) %>% 
  separate(Replicate, into = c("Replicate", "Technical"), sep = "(?<=\\d)(?=[a-z]?)", remove = FALSE) %>% 
  group_by(kit, Replicate) %>% 
  mutate(CV = ((sd(mg_Fe_per_L)/mean(mg_Fe_per_L))*100)) %>% 
  mutate(range = max(mg_Fe_per_L) - min(mg_Fe_per_L)) %>% 
  #distinct(kit_rep, .keep_all = TRUE) %>% 
  mutate(flag = case_when(
    CV <= 10 | range < 0.04 ~ "fine",
    CV >= 10 & range >= 0.04~ "flag"
  )) 
  
## After rerunning with 5 extra samples, remove samples to 3 if range > 0.04 and CV > 10 ####
  
  samples_dist <- data_flag_conc %>% 
    rename(cv.before.removal = CV) %>% 
    select(-c(flag, analysis)) %>% 
    unite(kit_rep, c("kit", "Replicate"), sep = "_", remove = FALSE)
  
  samples_removed_final = as.data.frame(matrix(NA, ncol = 13, nrow = 1))
  
  colnames(samples_removed_final) = c("sample_label", "Project", "kit_rep", "kit", "Analysis", "mg_Fe_per_L", "Replicate", "Technical", "cv.before.removal", "range", "cv.after.removal", "range.after.removal", "flag")
  
  unique.samples = unique(samples_dist$kit_rep)
  
  for (i in 1:length(unique.samples)) {
    
    data_subset = subset(samples_dist, samples_dist$kit_rep == unique.samples[i])
    
    conc.temp = as.numeric(data_subset$mg_Fe_per_L)
    
    conc.temp.sd <- sd(conc.temp)
    conc.temp.mean <- mean(conc.temp)
    CV = (conc.temp.sd/conc.temp.mean)*100
    range = max(conc.temp) - min(conc.temp)
    
    #looping to get 3 out of 8 best samples
    for (sample.reduction in 1:5)  {
      
      if (conc.temp.mean == 0) {
        
        CV = 0
        range = 0
        
      }
      
      else if (length(conc.temp) > 3 & CV >= 10 & range > 0.04) {
        
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
        conc.temp.range <- max(conc.temp) - min(conc.temp)
        range = conc.temp.range
        
      }
      
    }
    
    if(length(conc.temp) >= 3) {
      
      if(CV > 10 & range > 0.04){
        
        samples_removed_combined = as.data.frame(conc.temp)
        
        samples_remove = merge(samples_removed_combined, data_subset, by.x = "conc.temp", by.y = "mg_Fe_per_L", all.x = TRUE)
        
        samples_remove = samples_remove[!duplicated(samples_remove$sample_label), ]
        
        samples_remove_omit = merge(samples_remove, data_subset, by = "sample_label", all = TRUE)
        
        
        samples_remove_omit$cv.after.removal = as.numeric(abs((sd(conc.temp)/mean(conc.temp))*100))
        
        samples_remove_omit$range.after.removal = as.numeric(max(conc.temp) - min(conc.temp))
  
        samples_remove_omit$flag[which(samples_remove_omit$cv.after.removal >= 10 & samples_remove_omit$range.after.removal >= 0.04)] = "Samples too Variable"
        
        samples_remove_omit <- samples_remove_omit %>% 
          dplyr::select(-c(Project.x, kit_rep.x, kit.x, Analysis.x,Replicate.x, Technical.x, cv.before.removal.x, range.x, conc.temp)) %>% 
          rename(Project = Project.y) %>% 
          rename(kit_rep = kit_rep.y) %>% 
          rename(kit = kit.y) %>% 
          rename(Analysis = Analysis.y) %>% 
          rename(Replicate = Replicate.y) %>% 
          rename(Technical = Technical.y) %>% 
          rename(cv.before.removal = cv.before.removal.y) %>% 
          rename(range = range.y) %>% 
          relocate(range.after.removal, .after = cv.after.removal)
        
      }
      
      else {
        
        samples_removed_combined = as.data.frame(conc.temp)
        
        samples_remove = merge(samples_removed_combined, data_subset, by.x = "conc.temp", by.y = "mg_Fe_per_L", all.x = TRUE)
        
        samples_remove = samples_remove[!duplicated(samples_remove$sample_label), ]
        
        samples_remove_omit = merge(samples_remove, data_subset, by = "sample_label", all = TRUE)
        
        
        samples_remove_omit$cv.after.removal = as.numeric(abs((sd(conc.temp)/mean(conc.temp))*100))
        
        samples_remove_omit$range.after.removal = as.numeric(max(conc.temp) - min(conc.temp))
        
        samples_remove_omit$flag[is.na(samples_remove_omit$conc.temp)] = "OMIT"
        
        samples_remove_omit <- samples_remove_omit %>% 
          dplyr::select(-c(Project.x, kit_rep.x, kit.x, Analysis.x,Replicate.x, Technical.x, cv.before.removal.x, range.x, conc.temp)) %>% 
          rename(Project = Project.y) %>% 
          rename(kit_rep = kit_rep.y) %>% 
          rename(kit = kit.y) %>% 
          rename(Analysis = Analysis.y) %>% 
          rename(Replicate = Replicate.y) %>% 
          rename(Technical = Technical.y) %>% 
          rename(cv.before.removal = cv.before.removal.y) %>% 
          rename(range = range.y) %>% 
          relocate(range.after.removal, .after = cv.after.removal)
        
        
      }
      
    }
    
    samples_removed_final = rbind(samples_remove_omit, samples_removed_final)
    
    rm('conc.temp')
  }


  #final data corrected for high CV's - still needs to be processed for correcting for solid/solution ratio
  
  final_data <- samples_removed_final %>% 
    select(c(sample_label, mg_Fe_per_L, flag)) %>% 
    separate(col = sample_label, into = c("Project", "kit", "analysis"), sep = "_", remove = FALSE) %>% 
    separate(col = analysis, into = c("Analysis", "Replicate"), sep = "-", remove = FALSE) %>%
    separate(Replicate, into = c("Replicate", "Technical"), sep = "(?<=\\d)(?=[a-z]?)") %>% 
    unite(col = Sample_ID, c("Project", "kit", "Analysis"), sep = "_") %>% 
    unite(col = Sample_ID, c("Sample_ID", "Replicate"), sep = "-") %>% 
    select(-c(analysis))
  

### pull in moisture data to correct to mg Fe/kg dry sediment

moisture <- read.csv(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/INC/03_ProcessedData/ECA_Drying_Masses_merged_by_laan208_on_2023-05-15.csv"))


#merge moisture and Fe samples
final_data_moi <- moisture %>% 
  mutate(Sample_ID = str_replace(Sample_Name, "INC", "SFE")) %>% 
  distinct(Sample_ID, .keep_all = TRUE) %>% 
dplyr::select(-c(Sample_Name, Date))
  
merged <- merge(final_data_moi, final_data, by = "Sample_ID")


###assuming 40 mL water in 50 mL vial

merged$mg_Fe_per_kg_sediment <- merged$mg_Fe_per_L*(0.04)*(1/merged$dry_wt_sed_g)*1000

processed.data <- merged %>% 
  dplyr::select(sample_label,mg_Fe_per_L, mg_Fe_per_kg_sediment) 


write.csv(processed.data, "C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/FE/03_ProcessedData/20230110_Data_Processed_SFE_ECA_EC_1-270/20230110_Data_Processed_SFE_ECA_EC_1-270.csv")
