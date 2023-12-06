#### Sensitivity Analysis For ECA removals ####

library(dplyr)
library(tidyverse)

rm(list=ls());graphics.off()

#### Read in Data

#Individual samples 
all_data <- read.csv("C:/Github/ECA_Multireactor_Incubations/Data/Cleaned Data/All_ECA_Data.csv",header = TRUE) %>% 
  dplyr::select(-c(X))

#Summary Data 

sum_data <- read.csv("C:/Github/ECA_Multireactor_Incubations/Data/Cleaned Data/Summary_ECA_Data.csv",header = TRUE) %>% 
  dplyr::select(-c(X))

#Effect Size Data
effect_data <- read.csv("C:/Github/ECA_Multireactor_Incubations/Data/Cleaned Data/Effect_ECA_Data.csv",header = TRUE) %>% 
  dplyr::select(-c(X))

#Start SA

cv <- all_data %>% 
  mutate(Treat = case_when(grepl("W",Sample_Name)~"Wet",
                           grepl("D", Sample_Name) ~"Dry")) %>% 
  filter(Respiration_Rate_mg_DO_per_L_per_H > -9999) %>% 
  group_by(Sample_ID, Treat) %>% 
  mutate(cv_before_removal_L = abs(sd(Respiration_Rate_mg_DO_per_L_per_H)/mean(Respiration_Rate_mg_DO_per_L_per_H))*100) %>% 
  mutate(cv_before_removal_kg = abs(sd(Respiration_Rate_mg_DO_per_kg_per_H)/mean(Respiration_Rate_mg_DO_per_kg_per_H))*100) %>% 
  mutate(cv_diff = cv_before_removal_L - cv_before_removal_kg)


slope.outliers <- cv %>% 
  unite(kit_treat, Sample_ID, Treat, sep = "_", remove = FALSE) %>% 
  mutate(Respiration_Rate_mg_DO_per_L_per_H = as.numeric(Respiration_Rate_mg_DO_per_L_per_H)) %>%
  group_by(Sample_ID, Treat) %>% 
  mutate(Mean_Slope_All_L = mean(Respiration_Rate_mg_DO_per_L_per_H)) %>% 
  mutate(Mean_Slope_All_kg = mean(Respiration_Rate_mg_DO_per_kg_per_H)) %>% 
  ungroup() %>% 
  dplyr::select(c(Sample_Name, Respiration_Rate_mg_DO_per_L_per_H,Respiration_Rate_mg_DO_per_kg_per_H, Mean_Slope_All_L, Mean_Slope_All_kg, cv_before_removal_L, cv_before_removal_kg, kit_treat)) 

slope.outliers$flag <- NA

slope.final <- as.data.frame(matrix(NA, ncol = 11, nrow =1))

colnames(slope.final) = c("slope.temp","Sample_Name", "kit_treat", "Respiration_Rate_mg_DO_per_kg_per_H", "Mean_Slope_All_L", "Mean_Slope_All_kg", "cv_before_removal_L", "cv_before_removal_kg", "cv_after_removal", "Mean_Slope_Removed","flag")

unique.samples = unique(slope.outliers$kit_treat)

#try 0, 10, 30, 50, 100, export histograms of removals, effect sizes

cv.threshold = 0

for (i in 1:length(unique.samples)) {
  
  ## Subset replicates
  data_subset = subset(slope.outliers, slope.outliers$kit_treat == unique.samples[i])
  
  ## Pull out Rate values
  slope.temp = as.numeric(data_subset$Respiration_Rate_mg_DO_per_L_per_H)
  
  ## Calculate standard deviation, average, and coefficient of variation of rates
  slope.temp.sd <- sd(slope.temp)
  slope.temp.mean <- mean(slope.temp)
  CV = abs((slope.temp.sd/slope.temp.mean)*100)
  
  #looping to get 4 out of 5 best samples
  for (sample.reduction in 1:5)  {
    
    if (slope.temp.mean == 0) {
      
      CV = 0
      
    }
    
    else if (length(slope.temp) > 4 & CV >= cv.threshold) {
      
      dist.temp = as.matrix(abs(dist(slope.temp)))
      dist.comp = numeric()
      
      for(slope.now in 1:ncol(dist.temp)) {
        
        dist.comp = rbind(dist.comp,c(slope.now,sum(dist.temp[,slope.now])))
        
      }
      
      dist.comp[,2] = as.numeric(dist.comp[,2])
      slope.temp = slope.temp[-which.max(dist.comp[,2])]
      
      slope.temp.sd <- sd(slope.temp)
      slope.temp.mean <- mean(slope.temp)
      slope.temp.cv <- abs((slope.temp.sd/slope.temp.mean)*100)
      CV = slope.temp.cv
      slope.temp.range <- max(slope.temp) - min(slope.temp)
      range = slope.temp.range
      
    }
  }
  
  if (length(slope.temp) >= 3) {
    
    slope.combined <- as.data.frame(slope.temp)
    
    slope.removed <- merge(slope.combined, data_subset, by.x = "slope.temp", by.y = "Respiration_Rate_mg_DO_per_L_per_H", all.x = TRUE)
    
    slope.removed <- slope.removed[!duplicated(slope.removed$Sample_Name), ]
    
    slope.removed$cv_after_removal = as.numeric(abs((sd(slope.temp)/mean(slope.temp))*100))
    
    slope.removed$Mean_Slope_Removed = as.numeric(mean(slope.temp))
    
  }
  
  slope.final = rbind(slope.removed, slope.final)
  
  rm('slope.temp')
}

## This data frame has removed samples
slope.final.flag <- slope.final %>% 
  mutate(flag = if_else(cv_before_removal_L < cv_after_removal, "Issue in dropping sample", "N/A")) %>% 
  rename(Respiration_Rate_mg_DO_per_L_per_H = slope.temp) %>% 
  relocate(Respiration_Rate_mg_DO_per_L_per_H, .after = Sample_Name)

write.csv(slope.final.flag, paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/CV_mg_L_",cv.threshold,"percent_Removed_Respiration_Rates",Sys.Date(),".csv"), row.names = F)

corr_samples <- left_join(slope.final.flag, all_data, by = c("Sample_Name", "Respiration_Rate_mg_DO_per_L_per_H", "Respiration_Rate_mg_DO_per_kg_per_H")) %>% 
  dplyr::select(-c(Respiration_Rate_mg_DO_per_kg_per_H, Mean_Slope_All_L, Mean_Slope_All_kg, cv_before_removal_L, cv_before_removal_kg, kit_treat, flag, Sample_ID, Mean_Slope_Removed, cv_after_removal)) %>% 
  drop_na() %>% 
  mutate(Respiration_Rate_mg_DO_per_L_per_H = abs(Respiration_Rate_mg_DO_per_L_per_H)) %>% 
column_to_rownames("Sample_Name")

all_samples_corr <- cor(corr_samples, method = "spearman")

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/", as.character(Sys.Date()),"_All_Samples_Correlation_Matrix_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

corrplot(all_samples_corr,type = "upper", tl.col = "black", tl.cex = 1.2, cl.cex = 1,  title = "All Samples Correlation")

dev.off()
