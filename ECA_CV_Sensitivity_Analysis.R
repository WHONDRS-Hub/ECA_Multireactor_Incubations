#### Sensitivity Analysis For ECA removals ####

library(dplyr)
library(tidyverse)
library(corrplot)

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
rem.threshold = 3

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
    
    else if (length(slope.temp) > rem.threshold & CV >= cv.threshold) {
      
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

#cv.threshold = "No_Removals"

write.csv(slope.final.flag, paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/", cv.threshold, "_perc_CV/n=3/CV_mg_L_",cv.threshold,"percent_Removed_Respiration_Rates",Sys.Date(),".csv"), row.names = F)

corr_samples <- left_join(slope.final.flag, all_data, by = c("Sample_Name", "Respiration_Rate_mg_DO_per_L_per_H", "Respiration_Rate_mg_DO_per_kg_per_H")) %>% 
  dplyr::select(-c(Respiration_Rate_mg_DO_per_kg_per_H, Mean_Slope_All_L, Mean_Slope_All_kg, cv_before_removal_L, cv_before_removal_kg, kit_treat, flag, Sample_ID, Mean_Slope_Removed, cv_after_removal, Dry_Sediment_Mass_g, Initial_Water_mass_g, Final_Water_mass_g, Fe_mg_per_kg, Percent_Tot_Sand)) %>% 
  drop_na() %>% 
  mutate(Respiration_Rate_mg_DO_per_L_per_H = abs(Respiration_Rate_mg_DO_per_L_per_H)) %>% 
  rename(`Respiration (mg/L)` = Respiration_Rate_mg_DO_per_L_per_H) %>% 
  rename(`Fe (mg/L)` = Fe_mg_per_L) %>% 
    rename(`% Fine Sand` = Percent_Fine_Sand) %>% 
    rename(`% Med. Sand` = Percent_Med_Sand) %>% 
    rename(`% Coarse Sand` = Percent_Coarse_Sand) %>% 
    rename(`% Clay` = Percent_Clay) %>% 
    rename(`% Silt` = Percent_Silt) %>% 
    rename(SSA = average_ssa) %>% 
    rename(`In.Grav.Moi.` = Initial_Gravimetric_Water) %>% 
    rename(`Fin.Grav.Moi.` = Final_Gravimetric_Water) %>% 
    rename(`LostGrav.Moi.` = Lost_Gravimetric_Water) %>% 
column_to_rownames("Sample_Name")

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/", cv.threshold, "_perc_CV/n=", rem.threshold,"/_All_Samples_Histogram_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(gather(corr_samples, cols, value), aes(x = value)) + 
  geom_histogram() + 
  facet_wrap(.~cols, scales = 'free_x')

dev.off()

#all_samples_corr <- cor(corr_samples, method = "spearman")

## Correlation Matrix Function (Spearman) ####
panel.cor.spear <- function(x, y, digits=2, prefix="", cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r = (cor(x, y, method = c("spearman")))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) {cex.cor <- 0.8/strwidth(txt)} else {cex = cex.cor}
  text(0.5, 0.5, txt, cex = cex.cor * (1 + r)/1)
  
  # if(missing(cex.cor)) {cex <- 1.2/strwidth(txt)} else {cex = cex.cor}
  # text(0.5, 0.5, txt, cex = cex * sin(sqrt(abs(r))))
  
  test <- cor.test(x,y)
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " "))
  #text(0.5, 0.5, txt, cex = cex * r)
  text(.5, .8, Signif, cex=cex, col=2)
  
}

## Correlation Matrix Function (Pearson) ####
panel.cor.pear <- function(x, y, digits=2, prefix="", cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r = (cor(x, y, method = c("pearson")))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) {cex.cor <- 0.8/strwidth(txt)} else {cex = cex.cor}
  text(0.5, 0.5, txt, cex = cex.cor * (1 + r)/1)
  
  # if(missing(cex.cor)) {cex <- 1.2/strwidth(txt)} else {cex = cex.cor}
  # text(0.5, 0.5, txt, cex = cex * sin(sqrt(abs(r))))
  
  test <- cor.test(x,y)
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " "))
  #text(0.5, 0.5, txt, cex = cex * r)
  text(.5, .8, Signif, cex=cex, col=2)
  
}

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/", cv.threshold, "_perc_CV/n=",rem.threshold,"/_All_Samples_Correlation_Matrix_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

pairs(corr_samples, lower.panel = panel.smooth,upper.panel = panel.cor, gap = 0, cex.labels = 0.5, cex = .75)

#corrplot(all_samples_corr,type = "upper", tl.col = "black", tl.cex = 1.2, cl.cex = 1,  title = "All Samples Correlation")

dev.off()

log_corr_samples <- corr_samples %>% 
  mutate(`Log Respiration (mg/L)` = log10(`Respiration (mg/L)` + 0.5*min(`Respiration (mg/L)`[`Respiration (mg/L)` != min(`Respiration (mg/L)`)]))) %>% 
  mutate(`Log Fe (mg/L)` = log10(`Fe (mg/L)` + 0.5*min(`Fe (mg/L)`[`Fe (mg/L)` != min(`Fe (mg/L)`)]))) %>% 
  mutate(`Log % Fine Sand` = log10(`% Fine Sand` + 0.5*min(`% Fine Sand`))) %>% 
  mutate(`Log % Med. Sand` = log10(`% Med. Sand` + 0.5*min(`% Med. Sand`))) %>% 
  mutate(`Log % Coarse Sand` = log10(`% Coarse Sand` + 0.5*min(`% Coarse Sand`))) %>% 
  mutate(`Log % Clay` = log10(`% Clay` + 0.5*min(`% Clay`[`% Clay` != min(`% Clay`)])))



       
                
                
               

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/", cv.threshold, "/", as.character(Sys.Date()),"_Log_All_Samples_Correlation_Matrix_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(gather(log_corr_samples, cols, value), aes(x = value)) + 
  geom_histogram() + 
  facet_grid(.~cols)


log_samples_corr <- cor(log_corr_samples, method = "pearson")




png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/", cv.threshold, "/", as.character(Sys.Date()),"_Log_All_Samples_Correlation_Matrix_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

corrplot(log_samples_corr,type = "upper", tl.col = "black", tl.cex = 1.2, cl.cex = 1,  title = "Log All Samples Correlation")

dev.off()

## Wet/Dry Samples ####

all_samples_dry <- left_join(slope.final.flag, all_data, by = c("Sample_Name", "Respiration_Rate_mg_DO_per_L_per_H", "Respiration_Rate_mg_DO_per_kg_per_H")) %>% 
  filter(grepl("D", Sample_Name)) %>% 
  dplyr::select(-c(Respiration_Rate_mg_DO_per_kg_per_H, Mean_Slope_All_L, Mean_Slope_All_kg, cv_before_removal_L, cv_before_removal_kg, kit_treat, flag, Sample_ID, Mean_Slope_Removed, cv_after_removal, Dry_Sediment_Mass_g, Initial_Water_mass_g, Final_Water_mass_g)) %>% 
  drop_na() %>% 
  mutate(Respiration_Rate_mg_DO_per_L_per_H = abs(Respiration_Rate_mg_DO_per_L_per_H))%>% 
  column_to_rownames("Sample_Name")

all_samples_dry_corr <- cor(all_samples_dry,method = "spearman")

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/", as.character(Sys.Date()),"_Dry_Samples_Correlation_Matrix_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

corrplot(all_samples_dry_corr,type = "upper", tl.col = "black", tl.cex = 1.2, cl.cex = 1,  title = "Dry Samples Correlation")

dev.off()

all_samples_wet <- left_join(slope.final.flag, all_data, by = c("Sample_Name", "Respiration_Rate_mg_DO_per_L_per_H", "Respiration_Rate_mg_DO_per_kg_per_H")) %>% 
  filter(grepl("W", Sample_Name)) %>% 
  dplyr::select(-c(Respiration_Rate_mg_DO_per_kg_per_H, Mean_Slope_All_L, Mean_Slope_All_kg, cv_before_removal_L, cv_before_removal_kg, kit_treat, flag, Sample_ID, Mean_Slope_Removed, cv_after_removal, Dry_Sediment_Mass_g, Initial_Water_mass_g, Final_Water_mass_g)) %>% 
  drop_na() %>% 
  mutate(Respiration_Rate_mg_DO_per_L_per_H = abs(Respiration_Rate_mg_DO_per_L_per_H))%>% 
  column_to_rownames("Sample_Name")

all_samples_wet_corr <- cor(all_samples_wet,method = "spearman")

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/", as.character(Sys.Date()),"_Wet_Samples_Correlation_Matrix_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

corrplot(all_samples_wet_corr,type = "upper", tl.col = "black", tl.cex = 1.2, cl.cex = 1,  title = "Wet Samples Correlation")

dev.off()

#####

means <- left_join(slope.final.flag, all_data, by = c("Sample_Name", "Respiration_Rate_mg_DO_per_L_per_H", "Respiration_Rate_mg_DO_per_kg_per_H")) %>% 
  dplyr::select(-c(Respiration_Rate_mg_DO_per_kg_per_H, Mean_Slope_All_L, Mean_Slope_All_kg, cv_before_removal_L, cv_before_removal_kg, kit_treat, flag, Sample_ID, Mean_Slope_Removed, cv_after_removal, Dry_Sediment_Mass_g, Initial_Water_mass_g, Final_Water_mass_g)) %>% 
  drop_na() %>% 
  mutate(Respiration_Rate_mg_DO_per_L_per_H = abs(Respiration_Rate_mg_DO_per_L_per_H)) %>% 
  separate(Sample_Name, c("EC", "kit", "Treat"), remove = FALSE) %>% 
  mutate(Treat = case_when(grepl("W",Sample_Name)~"Wet",
                           grepl("D", Sample_Name) ~"Dry")) %>% 
  group_by(kit, Treat) %>% 
  summarise(across(where(is.numeric), list(mean = mean), na.rm = TRUE))%>% 
  unite(kit_treat, c("kit", "Treat")) %>% 
  column_to_rownames("kit_treat")

mean_samples_corr <- cor(means, method = "spearman")

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/", cv.threshold, "/", as.character(Sys.Date()),"_Mean_Samples_Correlation_Matrix_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

corrplot(mean_samples_corr,type = "upper", tl.col = "black", tl.cex = 1.2, cl.cex = 1,  title = "Mean Samples Correlation")

dev.off()

log_means <- log(means + 1)

log_means_corr <- cor(log_means, method = "pearson")

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/", cv.threshold, "/", as.character(Sys.Date()),"_Log_Mean_Samples_Correlation_Matrix_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

corrplot(log_means_corr,type = "upper", tl.col = "black", tl.cex = 1.2, cl.cex = 1,  title = "Log Mean Samples Correlation")

dev.off()

## Wet/Dry Means ####

wet_means <- left_join(slope.final.flag, all_data, by = c("Sample_Name", "Respiration_Rate_mg_DO_per_L_per_H", "Respiration_Rate_mg_DO_per_kg_per_H")) %>% 
  dplyr::select(-c(Respiration_Rate_mg_DO_per_kg_per_H, Mean_Slope_All_L, Mean_Slope_All_kg, cv_before_removal_L, cv_before_removal_kg, kit_treat, flag, Sample_ID, Mean_Slope_Removed, cv_after_removal, Dry_Sediment_Mass_g, Initial_Water_mass_g, Final_Water_mass_g)) %>% 
  drop_na() %>% 
  mutate(Respiration_Rate_mg_DO_per_L_per_H = abs(Respiration_Rate_mg_DO_per_L_per_H)) %>% 
  separate(Sample_Name, c("EC", "kit", "Treat"), remove = FALSE) %>% 
  mutate(Treat = case_when(grepl("W",Sample_Name)~"Wet",
                           grepl("D", Sample_Name) ~"Dry")) %>% 
  group_by(kit, Treat) %>% 
  summarise(across(where(is.numeric), list(mean = mean), na.rm = TRUE))%>% 
  ungroup() %>% 
  unite(kit_treat, c("kit", "Treat")) %>% 
  filter(grepl("Wet", kit_treat)) %>% 
  column_to_rownames("kit_treat")

dry_means <- left_join(slope.final.flag, all_data, by = c("Sample_Name", "Respiration_Rate_mg_DO_per_L_per_H", "Respiration_Rate_mg_DO_per_kg_per_H")) %>% 
  dplyr::select(-c(Respiration_Rate_mg_DO_per_kg_per_H, Mean_Slope_All_L, Mean_Slope_All_kg, cv_before_removal_L, cv_before_removal_kg, kit_treat, flag, Sample_ID, Mean_Slope_Removed, cv_after_removal, Dry_Sediment_Mass_g, Initial_Water_mass_g, Final_Water_mass_g)) %>% 
  drop_na() %>% 
  mutate(Respiration_Rate_mg_DO_per_L_per_H = abs(Respiration_Rate_mg_DO_per_L_per_H)) %>% 
  separate(Sample_Name, c("EC", "kit", "Treat"), remove = FALSE) %>% 
  mutate(Treat = case_when(grepl("W",Sample_Name)~"Wet",
                           grepl("D", Sample_Name) ~"Dry")) %>% 
  group_by(kit, Treat) %>% 
  summarise(across(where(is.numeric), list(mean = mean), na.rm = TRUE))%>% 
  ungroup() %>% 
  unite(kit_treat, c("kit", "Treat")) %>% 
  filter(grepl("Dry", kit_treat)) %>% 
  column_to_rownames("kit_treat")

#####


effect <- left_join(slope.final.flag, all_data, by = c("Sample_Name", "Respiration_Rate_mg_DO_per_L_per_H", "Respiration_Rate_mg_DO_per_kg_per_H")) %>% 
  dplyr::select(-c(Respiration_Rate_mg_DO_per_kg_per_H, Mean_Slope_All_L, Mean_Slope_All_kg, cv_before_removal_L, cv_before_removal_kg, kit_treat, flag, Sample_ID, Mean_Slope_Removed, cv_after_removal, Dry_Sediment_Mass_g, Initial_Water_mass_g, Final_Water_mass_g)) %>% 
  drop_na() %>% 
  mutate(Respiration_Rate_mg_DO_per_L_per_H = abs(Respiration_Rate_mg_DO_per_L_per_H)) %>% 
  separate(Sample_Name, c("EC", "kit", "Treat"), remove = FALSE) %>% 
  mutate(Treat = case_when(grepl("W",Sample_Name)~"Wet",
                           grepl("D", Sample_Name) ~"Dry")) %>% 
  group_by(kit, Treat) %>% 
  filter(n() >= 4) %>% 
  ungroup() %>% 
  group_by(kit) %>% 
  filter(n() >= 8) %>% 
  ungroup() %>% 
  group_by(kit, Treat) %>% 
  summarise(across(where(is.numeric), list( mean = mean), na.rm = TRUE)) %>% 
  ungroup() %>% 
  group_by(kit) %>% 
  relocate(c(Percent_Fine_Sand_mean:average_ssa_mean), .after = Lost_Gravimetric_Water_mean) %>% 
  summarise(
    across(Respiration_Rate_mg_DO_per_L_per_H_mean:Lost_Gravimetric_Water_mean, diff), 
    across(Percent_Fine_Sand_mean:average_ssa_mean, mean))%>% 
  column_to_rownames("kit") %>% 
  rename(Effect_Size = Respiration_Rate_mg_DO_per_L_per_H_mean) %>%
  rename(Percent_Clay = Percent_Clay_mean) %>% 
  rename(Percent_Silt = Percent_Silt_mean) %>% 
  rename(SSA = average_ssa_mean) %>% 
  rename_with(~gsub("_Sand_mean", "_Sand", .x)) %>% 
  rename_with(~gsub("mean", "diff", .x)) 
  
  
effect_corr <- cor(effect, method = "spearman")

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/", cv.threshold, "/", as.character(Sys.Date()),"_Effect_Samples_Correlation_Matrix_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

corrplot(effect_corr,type = "upper", tl.col = "black", tl.cex = 1.2, cl.cex = 1,  title = "Effect Correlation")

dev.off()
  
log_effect <- log(effect + 1)

log_effect_corr <- cor(log_effect, method = "pearson")

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/", cv.threshold, "/", as.character(Sys.Date()),"_Log_Effect_Samples_Correlation_Matrix_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

corrplot(log_effect_corr,type = "upper", tl.col = "black", tl.cex = 1.2, cl.cex = 1,  title = "Log Effect Correlation")

dev.off()

