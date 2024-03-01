#### Sensitivity Analysis For ECA removals ####

library(dplyr)
library(tidyverse)
library(corrplot)
library(glmnet)
library(readxl)

rm(list=ls());graphics.off()

#### Read in Data

#Individual samples 
all_data <- read.csv("C:/Github/ECA_Multireactor_Incubations/Data/Cleaned Data/All_ECA_Data.csv",header = TRUE) %>% 
  dplyr::select(-c(X))

#ATP data

atp <- read.csv("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/ATP/03_ProcessedData/EC_ATP_ReadyForBoye_01-26-2024.csv") %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "ATP", "INC"))

all_data <- left_join(all_data, atp, by = "Sample_Name")


effect_scale <- effect_data %>% 
  column_to_rownames(var = "Sample_Name") %>% 
  scale()

## Start SA

#Remove RR with -9999, calculate CV for replicates for mg/L and mg/kg values

cv <- all_data %>% 
  mutate(Treat = case_when(grepl("W",Sample_Name)~"Wet",
                           grepl("D", Sample_Name) ~"Dry")) %>% 
  filter(Respiration_Rate_mg_DO_per_L_per_H > -9999) %>% 
  group_by(Sample_ID, Treat) %>% 
  mutate(cv_before_removal_L = abs(sd(Respiration_Rate_mg_DO_per_L_per_H)/mean(Respiration_Rate_mg_DO_per_L_per_H))*100) %>% 
  mutate(cv_before_removal_kg = abs(sd(Respiration_Rate_mg_DO_per_kg_per_H)/mean(Respiration_Rate_mg_DO_per_kg_per_H))*100) %>% 
  mutate(cv_diff = cv_before_removal_L - cv_before_removal_kg)

## Start data frame for removals, calculate means of replicate samples

slope.outliers <- cv %>% 
  unite(kit_treat, Sample_ID, Treat, sep = "_", remove = FALSE) %>% 
  mutate(Respiration_Rate_mg_DO_per_L_per_H = as.numeric(Respiration_Rate_mg_DO_per_L_per_H)) %>%
  group_by(Sample_ID, Treat) %>% 
  mutate(Mean_Slope_All_L = mean(Respiration_Rate_mg_DO_per_L_per_H)) %>% 
  mutate(Mean_Slope_All_kg = mean(Respiration_Rate_mg_DO_per_kg_per_H)) %>% 
  ungroup() %>% 
  dplyr::select(c(Sample_Name, Respiration_Rate_mg_DO_per_L_per_H,Respiration_Rate_mg_DO_per_kg_per_H, Mean_Slope_All_L, Mean_Slope_All_kg, cv_before_removal_L, cv_before_removal_kg, kit_treat)) %>% 
  drop_na(Respiration_Rate_mg_DO_per_kg_per_H)

slope.outliers$flag <- NA

slope.final <- as.data.frame(matrix(NA, ncol = 11, nrow =1))

colnames(slope.final) = c("slope.temp","Sample_Name", "kit_treat", "Respiration_Rate_mg_DO_per_L_per_H", "Mean_Slope_All_L", "Mean_Slope_All_kg", "cv_before_removal_L", "cv_before_removal_kg", "cv_after_removal", "Mean_Slope_Removed","flag")

unique.samples = unique(slope.outliers$kit_treat)

#try 0, 10, 30, 50, 100, export histograms of removals, effect sizes
cv.threshold = 50

#try keeping n = 3,4 samples
rem.threshold = 3

for (i in 1:length(unique.samples)) {
  
  ## Subset replicates
  data_subset = subset(slope.outliers, slope.outliers$kit_treat == unique.samples[i])
  
  ## Pull out Rate values
  #slope.temp = as.numeric(data_subset$Respiration_Rate_mg_DO_per_L_per_H)
  slope.temp = as.numeric(data_subset$Respiration_Rate_mg_DO_per_kg_per_H)
  
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
    
    slope.removed <- merge(slope.combined, data_subset, by.x = "slope.temp", by.y = "Respiration_Rate_mg_DO_per_kg_per_H", all.x = TRUE)
    
    slope.removed <- slope.removed[!duplicated(slope.removed$Sample_Name), ]
    
    slope.removed$cv_after_removal = as.numeric(abs((sd(slope.temp)/mean(slope.temp))*100))
    
    slope.removed$Mean_Slope_Removed = as.numeric(mean(slope.temp))
    
  }
  
  slope.final = rbind(slope.removed, slope.final)
  
  rm('slope.temp')
}

## This data frame has removed samples
slope.final.flag <- slope.final %>% 
  mutate(flag = if_else(cv_before_removal_kg < cv_after_removal, "Issue in dropping sample", "N/A")) %>% 
  rename(Respiration_Rate_mg_DO_per_kg_per_H = slope.temp) %>% 
  relocate(Respiration_Rate_mg_DO_per_kg_per_H, .after = Sample_Name)

cv.threshold = 50
rem.threshold = 3


#Write data frame with removed respiration rates

write.csv(slope.final.flag, paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/DO_per_kg/", cv.threshold, "_perc_CV/n=", rem.threshold,"/CV_mg_kg_",cv.threshold,"percent_Removed_Respiration_Rates",Sys.Date(),".csv"), row.names = F)

## All Samples ####

##Merge removed samples with other data
corr_samples <- left_join(all_data, slope.final.flag, by = c("Sample_Name", "Respiration_Rate_mg_DO_per_L_per_H", "Respiration_Rate_mg_DO_per_kg_per_H")) %>% 
  mutate(flag = if_else(is.na(Mean_Slope_All_L),  "Remove Outlier", "N/A")) %>% 
  filter(flag != "Remove Outlier") %>% 
  dplyr::select(-c(Respiration_Rate_mg_DO_per_kg_per_H, Mean_Slope_All_L, Mean_Slope_All_kg, cv_before_removal_L, cv_before_removal_kg, kit_treat, flag, Sample_ID, Mean_Slope_Removed, cv_after_removal, Dry_Sediment_Mass_g, Initial_Water_mass_g, Final_Water_mass_g, Fe_mg_per_kg, Percent_Tot_Sand, Material, Methods_Deviation, ATP_nanomol_per_L)) %>% 
  drop_na() %>% 
  filter(ATP_picomol_per_g != -9999) %>% 
  filter(Respiration_Rate_mg_DO_per_L_per_H != -9999) %>% 
  mutate(Respiration_Rate_mg_DO_per_L_per_H = abs(Respiration_Rate_mg_DO_per_L_per_H)) %>% 
  rename(`Respiration (mg/L)` = Respiration_Rate_mg_DO_per_L_per_H) %>% 
 # rename(`Fe (mg/kg)` = Fe_mg_per_kg) %>% 
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

#####
## Make histogram of all data
png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/DO_per_kg/", cv.threshold, "_perc_CV/n=", rem.threshold,"/All_Samples_Histogram_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(gather(corr_samples, cols, value), aes(x = value)) + 
  geom_histogram() + 
  facet_wrap(.~cols, scales = 'free_x')

dev.off()

#####
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

#####
png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/DO_per_kg/", cv.threshold, "_perc_CV/n=",rem.threshold,"/All_Samples_Correlation_Matrix_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

pairs(corr_samples, lower.panel = panel.smooth,upper.panel = panel.cor.spear, gap = 0, cex.labels = 0.5, cex = .75)

#corrplot(all_samples_corr,type = "upper", tl.col = "black", tl.cex = 1.2, cl.cex = 1,  title = "All Samples Correlation")

dev.off()

## Log All Samples ####

log_corr_samples <- corr_samples %>% 
  mutate(`Log Respiration (mg/L)` = log10(`Respiration (mg/L)` + 0.5*min(`Respiration (mg/L)`[`Respiration (mg/L)` != min(`Respiration (mg/L)`)]))) %>% 
  mutate(`Log Fe (mg/L)` = log10(Fe_mg_per_L + 0.5*min(Fe_mg_per_L[Fe_mg_per_L != min(Fe_mg_per_L)]))) %>% 
  mutate(`Log % Fine Sand` = log10(`% Fine Sand` + 0.5*min(`% Fine Sand`))) %>% 
  mutate(`Log % Med. Sand` = log10(`% Med. Sand` + 0.5*min(`% Med. Sand`))) %>% 
  mutate(`Log % Coarse Sand` = log10(`% Coarse Sand` + 0.5*min(`% Coarse Sand`))) %>% 
  mutate(`Log % Clay` = log10(`% Clay` + 0.5*min(`% Clay`[`% Clay` != min(`% Clay`)]))) %>% 
  mutate(`Log % Silt` = log10(`% Silt` + 0.5*min(`% Silt`[`% Silt` != min(`% Silt`)]))) %>% 
  mutate(`Log SSA` = log10(`SSA` + 0.5*min(`SSA`))) %>% 
  mutate(`Log SpC` = log10(`SpC` + 0.5*min(`SpC`))) %>% 
  mutate(`Log Temp` = log10(`Temp` + 0.5*min(`Temp`))) %>%
  mutate(`Log pH` = log10(`pH` + 0.5*min(`pH`))) %>% 
  mutate(`Log In.Grav.Moi.` = log10(`In.Grav.Moi.` + 0.5*min(`In.Grav.Moi.`))) %>% 
  mutate(`Log Fin.Grav.Moi.` = log10(`Fin.Grav.Moi.` + 0.5*min(`Fin.Grav.Moi.`))) %>% 
  mutate(`Log Fin.Grav.Moi.` = log10(`Fin.Grav.Moi.` + 0.5*min(`Fin.Grav.Moi.`))) %>% 
  mutate(`Log LostGrav.Moi.` = log10(`LostGrav.Moi.` + 0.5*min(`LostGrav.Moi.`))) %>% 
  mutate(`Log ATP (pmol/g)` = log10(ATP_picomol_per_g + 0.5*min(ATP_picomol_per_g))) %>% 
   dplyr::select(c(`Log Respiration (mg/L)`, `Log Fe (mg/L)`, `Log ATP (pmol/g)`, `Log % Fine Sand`, `Log % Med. Sand`, `Log % Coarse Sand`, `Log % Clay`, `Log % Silt`, `Log SSA`, `Log SpC`, `Log Temp`, `Log pH`, `Log In.Grav.Moi.`, `Log Fin.Grav.Moi.`, `Log LostGrav.Moi.`))

# histogram of log transformed samples (log + 1/2 of minimum to get at 0 values)

#question - do this across the board even if no 0's?

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/DO_per_kg/", cv.threshold, "_perc_CV/n=", rem.threshold,"/Log_All_Samples_Histogram_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(gather(log_corr_samples, cols, value), aes(x = value)) + 
  geom_histogram() + 
  facet_wrap(.~cols, scales = 'free_x')

dev.off()

#log_samples_corr <- cor(log_corr_samples, method = "pearson")

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/DO_per_kg/", cv.threshold, "_perc_CV/n=",rem.threshold,"/Log_All_Samples_Correlation_Matrix_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

pairs(corr_samples, lower.panel = panel.smooth,upper.panel = panel.cor.pear, gap = 0, cex.labels = 0.5, cex = .75)

#corrplot(log_samples_corr,type = "upper", tl.col = "black", tl.cex = 1.2, cl.cex = 1,  title = "Log All Samples Correlation")

dev.off()

## Wet/Dry All Samples ####

all_samples_dry <- corr_samples %>% 
  rownames_to_column("Sample_Name") %>% 
  filter(grepl("D", Sample_Name)) %>% 
  column_to_rownames("Sample_Name")

## Make histogram of all data
png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/DO_per_kg/", cv.threshold, "_perc_CV/n=", rem.threshold,"/All_Dry_Samples_Histogram_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(gather(all_samples_dry, cols, value), aes(x = value)) + 
  geom_histogram() + 
  facet_wrap(.~cols, scales = 'free_x')

dev.off()

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/DO_per_kg/", cv.threshold, "_perc_CV/n=",rem.threshold,"/All_Dry_Samples_Correlation_Matrix_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

pairs(all_samples_dry, lower.panel = panel.smooth,upper.panel = panel.cor.spear, gap = 0, cex.labels = 0.5, cex = .75)

#corrplot(all_samples_corr,type = "upper", tl.col = "black", tl.cex = 1.2, cl.cex = 1,  title = "All Samples Correlation")

dev.off()

log_all_dry <- log_corr_samples %>% 
  rownames_to_column("Sample_Name") %>% 
  filter(grepl("D", Sample_Name)) %>% 
  column_to_rownames("Sample_Name")
  
png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/DO_per_kg/", cv.threshold, "_perc_CV/n=", rem.threshold,"/Log_All_Dry_Samples_Histogram_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(gather(log_all_dry, cols, value), aes(x = value)) + 
  geom_histogram() + 
  facet_wrap(.~cols, scales = 'free_x')

dev.off()

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/DO_per_kg/", cv.threshold, "_perc_CV/n=",rem.threshold,"/Log_All_Dry_Samples_Correlation_Matrix_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

pairs(log_all_dry, lower.panel = panel.smooth,upper.panel = panel.cor.pear, gap = 0, cex.labels = 0.5, cex = .75)

dev.off()

all_samples_wet <-  corr_samples %>% 
  rownames_to_column("Sample_Name") %>% 
  filter(grepl("W", Sample_Name)) %>% 
  column_to_rownames("Sample_Name")

## Make histogram of all data
png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/DO_per_kg/", cv.threshold, "_perc_CV/n=", rem.threshold,"/All_Wet_Samples_Histogram_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(gather(all_samples_wet, cols, value), aes(x = value)) + 
  geom_histogram() + 
  facet_wrap(.~cols, scales = 'free_x')

dev.off()

#all_samples_dry_corr <- cor(all_samples_dry,method = "spearman")

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/DO_per_kg/", cv.threshold, "_perc_CV/n=",rem.threshold,"/All_Wet_Samples_Correlation_Matrix_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

pairs(all_samples_wet, lower.panel = panel.smooth,upper.panel = panel.cor.spear, gap = 0, cex.labels = 0.5, cex = .75)

#corrplot(all_samples_corr,type = "upper", tl.col = "black", tl.cex = 1.2, cl.cex = 1,  title = "All Samples Correlation")

dev.off()

log_all_wet <- log_corr_samples %>% 
  rownames_to_column("Sample_Name") %>% 
  filter(grepl("W", Sample_Name)) %>% 
  dplyr::select(-c(`Log LostGrav.Moi.`)) %>% 
  column_to_rownames("Sample_Name")

## Make histogram of all data
png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/DO_per_kg/", cv.threshold, "_perc_CV/n=", rem.threshold,"/Log_All_Wet_Samples_Histogram_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(gather(log_all_wet, cols, value), aes(x = value)) + 
  geom_histogram() + 
  facet_wrap(.~cols, scales = 'free_x')

dev.off()

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/DO_per_kg/", cv.threshold, "_perc_CV/n=",rem.threshold,"/Log_All_Wet_Samples_Correlation_Matrix_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

pairs(log_all_wet, lower.panel = panel.smooth,upper.panel = panel.cor.pear, gap = 0, cex.labels = 0.5, cex = .75)

dev.off()

## Means ####

mean_resp <- corr_samples %>%
  rownames_to_column("Sample_Name") %>% 
  select(c(Sample_Name, `Respiration (mg/L)`)) %>% 
  separate(Sample_Name, c("EC", "kit", "Treat"), remove = FALSE) %>% 
  mutate(Treat = case_when(grepl("W",Sample_Name)~"Wet",
                           grepl("D", Sample_Name) ~"Dry")) %>% 
  group_by(kit, Treat) %>% 
  summarise(across(where(is.numeric), list(mean = mean), na.rm = TRUE)) %>% 
  unite(kit_treat, c("kit", "Treat")) 

means <- all_data %>% 
  select(c(Sample_Name, Fe_mg_per_L, ATP_picomol_per_g, SpC, Temp, pH, Initial_Gravimetric_Water, Final_Gravimetric_Water, Lost_Gravimetric_Water, Percent_Fine_Sand, Percent_Med_Sand, Percent_Coarse_Sand, Percent_Clay, Percent_Silt, average_ssa)) %>% 
  filter(ATP_picomol_per_g != -9999) %>% 
  separate(Sample_Name, c("EC", "kit", "Treat"), remove = FALSE) %>% 
  mutate(Treat = case_when(grepl("W",Sample_Name)~"Wet",
                           grepl("D", Sample_Name) ~"Dry")) %>% 
  group_by(kit, Treat) %>% 
  summarise(across(where(is.numeric), list(mean = mean), na.rm = TRUE)) %>% 
  unite(kit_treat, c("kit", "Treat")) %>% 
  left_join(mean_resp, by = "kit_treat") %>% 
  drop_na()



  

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/DO_per_kg/", cv.threshold, "_perc_CV/n=", rem.threshold,"/Mean_Samples_Histogram_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(gather(means, cols, value), aes(x = value)) + 
  geom_histogram() + 
  facet_wrap(.~cols, scales = 'free_x')

dev.off()

#mean_samples_corr <- cor(means, method = "spearman")

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/DO_per_kg/", cv.threshold, "_perc_CV/n=",rem.threshold,"/Mean_Samples_Correlation_Matrix_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

pairs(means, lower.panel = panel.smooth,upper.panel = panel.cor.spear, gap = 0, cex.labels = 0.5, cex = .75)

#corrplot(all_samples_corr,type = "upper", tl.col = "black", tl.cex = 1.2, cl.cex = 1,  title = "All Samples Correlation")

dev.off()

## Log Mean Samples ####
log_means <- means %>% 
  mutate(`Log Mean Respiration (mg/L)` = log10(`Mean Respiration (mg/L)` + 0.5*min(`Mean Respiration (mg/kg)`[`Mean Respiration (mg/kg)` != min(`Mean Respiration (mg/kg)`)]))) %>% 
  mutate(`Log Mean Fe (mg/kg)` = log10(`Mean Fe (mg/kg)` + 0.5*min(`Mean Fe (mg/kg)`[`Mean Fe (mg/kg)` != min(`Mean Fe (mg/kg)`)]))) %>% 
  mutate(`Log % Fine Sand` = log10(`% Fine Sand` + 0.5*min(`% Fine Sand`))) %>% 
  mutate(`Log % Med. Sand` = log10(`% Med. Sand` + 0.5*min(`% Med. Sand`))) %>% 
  mutate(`Log % Coarse Sand` = log10(`% Coarse Sand` + 0.5*min(`% Coarse Sand`))) %>% 
  mutate(`Log % Clay` = log10(`% Clay` + 0.5*min(`% Clay`[`% Clay` != min(`% Clay`)]))) %>% 
  mutate(`Log % Silt` = log10(`% Silt` + 0.5*min(`% Silt`[`% Silt` != min(`% Silt`)]))) %>% 
  mutate(`Log SSA` = log10(`SSA` + 0.5*min(`SSA`))) %>% 
  mutate(`Log Mean SpC` = log10(`Mean SpC` + 0.5*min(`Mean SpC`))) %>% 
  mutate(`Log Mean Temp` = log10(`Mean Temp` + 0.5*min(`Mean Temp`))) %>%
  mutate(`Log Mean pH` = log10(`Mean pH` + 0.5*min(`Mean pH`))) %>% 
  mutate(`Log Mean In.Grav.Moi.` = log10(`Mean In.Grav.Moi.` + 0.5*min(`Mean In.Grav.Moi.`))) %>% 
  mutate(`Log Mean Fin.Grav.Moi.` = log10(`Mean Fin.Grav.Moi.` + 0.5*min(`Mean Fin.Grav.Moi.`))) %>% 
  mutate(`Log Mean Fin.Grav.Moi.` = log10(`Mean Fin.Grav.Moi.` + 0.5*min(`Mean Fin.Grav.Moi.`))) %>% 
  mutate(`Log Mean LostGrav.Moi.` = log10(`Mean LostGrav.Moi.` + 0.5*min(`Mean LostGrav.Moi.`))) %>% 
  dplyr::select(c(`Log Mean Respiration (mg/kg)`, `Log Mean Fe (mg/kg)`, `Log % Fine Sand`, `Log % Med. Sand`, `Log % Coarse Sand`, `Log % Clay`, `Log % Silt`, `Log SSA`, `Log Mean SpC`, `Log Mean Temp`, `Log Mean pH`, `Log Mean In.Grav.Moi.`, `Log Mean Fin.Grav.Moi.`, `Log Mean LostGrav.Moi.`))

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/DO_per_kg/", cv.threshold, "_perc_CV/n=", rem.threshold,"/Log_Mean_Samples_Histogram_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(gather(log_means, cols, value), aes(x = value)) + 
  geom_histogram() + 
  facet_wrap(.~cols, scales = 'free_x')

dev.off()

#mean_samples_corr <- cor(means, method = "spearman")

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/DO_per_kg/", cv.threshold, "_perc_CV/n=",rem.threshold,"/Log_Mean_Samples_Correlation_Matrix_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

pairs(log_means, lower.panel = panel.smooth,upper.panel = panel.cor.pear, gap = 0, cex.labels = 0.5, cex = .75)

#corrplot(all_samples_corr,type = "upper", tl.col = "black", tl.cex = 1.2, cl.cex = 1,  title = "All Samples Correlation")

dev.off()

#log_means_corr <- cor(log_means, method = "pearson")

## Wet/Dry Means ####

wet_means <- means %>% 
  rownames_to_column("kit_treat") %>% 
  filter(grepl("W", kit_treat)) %>% 
  dplyr::select(-c(`Mean Fin.Grav.Moi.`, `Mean LostGrav.Moi.`)) %>% 
  column_to_rownames("kit_treat")         
  
png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/DO_per_kg/", cv.threshold, "_perc_CV/n=", rem.threshold,"/Wet_Mean_Samples_Histogram_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(gather(wet_means, cols, value), aes(x = value)) + 
  geom_histogram() + 
  facet_wrap(.~cols, scales = 'free_x')

dev.off()

#mean_samples_corr <- cor(means, method = "spearman")

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/DO_per_kg/", cv.threshold, "_perc_CV/n=",rem.threshold,"/Wet_Mean_Samples_Correlation_Matrix_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

pairs(wet_means, lower.panel = panel.smooth,upper.panel = panel.cor.spear, gap = 0, cex.labels = 0.5, cex = .75)

#corrplot(all_samples_corr,type = "upper", tl.col = "black", tl.cex = 1.2, cl.cex = 1,  title = "All Samples Correlation")

dev.off()


dry_means <- means %>% 
  rownames_to_column("kit_treat") %>% 
  filter(grepl("D", kit_treat)) %>% 
  column_to_rownames("kit_treat")

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/DO_per_kg/", cv.threshold, "_perc_CV/n=", rem.threshold,"/Dry_Mean_Samples_Histogram_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(gather(dry_means, cols, value), aes(x = value)) + 
  geom_histogram() + 
  facet_wrap(.~cols, scales = 'free_x')

dev.off()

#mean_samples_corr <- cor(means, method = "spearman")

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/DO_per_kg/", cv.threshold, "_perc_CV/n=",rem.threshold,"/Dry_Mean_Samples_Correlation_Matrix_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

pairs(dry_means, lower.panel = panel.smooth,upper.panel = panel.cor.spear, gap = 0, cex.labels = 0.5, cex = .75)

#corrplot(all_samples_corr,type = "upper", tl.col = "black", tl.cex = 1.2, cl.cex = 1,  title = "All Samples Correlation")

dev.off()

## Log Mean Samples

log_wet_means <- log_means %>% 
  rownames_to_column("kit_treat") %>% 
  filter(grepl("W", kit_treat)) %>% 
  dplyr::select(-c(`Log Mean LostGrav.Moi.`, `Log Mean In.Grav.Moi.`, `Log Mean Fin.Grav.Moi.`)) %>% 
  column_to_rownames("kit_treat")

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/DO_per_kg/", cv.threshold, "_perc_CV/n=", rem.threshold,"/Log_Wet_Mean_Samples_Histogram_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(gather(log_wet_means, cols, value), aes(x = value)) + 
  geom_histogram() + 
  facet_wrap(.~cols, scales = 'free_x')

dev.off()

#mean_samples_corr <- cor(means, method = "spearman")

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/DO_per_kg/", cv.threshold, "_perc_CV/n=",rem.threshold,"/Log_Wet_Mean_Samples_Correlation_Matrix_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

pairs(log_wet_means, lower.panel = panel.smooth,upper.panel = panel.cor.pear, gap = 0, cex.labels = 0.5, cex = .75)

#corrplot(all_samples_corr,type = "upper", tl.col = "black", tl.cex = 1.2, cl.cex = 1,  title = "All Samples Correlation")

dev.off()


log_dry_means <- log_means %>% 
  rownames_to_column("kit_treat") %>% 
  filter(grepl("D", kit_treat)) %>% 
  column_to_rownames("kit_treat")
  
png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/DO_per_kg/", cv.threshold, "_perc_CV/n=", rem.threshold,"/Log_Dry_Mean_Samples_Histogram_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(gather(log_dry_means, cols, value), aes(x = value)) + 
  geom_histogram() + 
  facet_wrap(.~cols, scales = 'free_x')

dev.off()

#mean_samples_corr <- cor(means, method = "spearman")

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/DO_per_kg/", cv.threshold, "_perc_CV/n=",rem.threshold,"/Log_Dry_Mean_Samples_Correlation_Matrix_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

pairs(log_dry_means, lower.panel = panel.smooth,upper.panel = panel.cor.pear, gap = 0, cex.labels = 0.5, cex = .75)

#corrplot(all_samples_corr,type = "upper", tl.col = "black", tl.cex = 1.2, cl.cex = 1,  title = "All Samples Correlation")

dev.off()

#####
effect <- means %>% 
 # rownames_to_column("kit_treat") %>% 
  separate(kit_treat, c("kit", "treat"), remove = FALSE) %>% 
  group_by(kit) %>% 
  filter(n() >= 2) %>% 
  ungroup() %>% 
  group_by(kit) %>% 
 # relocate(c(`% Fine Sand`:SSA), .after = `ATP_picomol_per_g_mean`) %>% 
  relocate(c(Percent_Fine_Sand_mean:average_ssa_mean), .after = `Respiration (mg/L)_mean`) %>% 
  rename(`Mean Respiration (mg/L)` = `Respiration (mg/L)_mean`) %>% 
  relocate(`Mean Respiration (mg/L)`, .after= treat) %>% 
  summarise(
    across(`Mean Respiration (mg/L)`:`Lost_Gravimetric_Water_mean`, diff), 
    #across(`% Fine Sand`:SSA, mean))%>% 
    across(Percent_Fine_Sand_mean:average_ssa_mean, mean))%>% 
  column_to_rownames("kit") %>% 
  rename(`Effect Size` = `Mean Respiration (mg/L)`) %>% 
  rename(Fe_mg_per_L_diff = Fe_mg_per_L_mean) %>% 
  rename(`SpC Diff.` = SpC_mean) %>% 
  rename(`Temp Diff.` = Temp_mean) %>% 
  rename(`pH Diff.` = pH_mean) %>% 
  rename(`In.Grav.Moi. Diff.` =Initial_Gravimetric_Water_mean) %>% 
  rename(`Fin.Grav.Moi. Diff.` = Final_Gravimetric_Water_mean) %>%
  rename(`LostGrav.Moi. Diff.` = Lost_Gravimetric_Water_mean)

#effect_corr <- cor(effect, method = "spearman")

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/DO_per_kg/", cv.threshold, "_perc_CV/n=", rem.threshold,"/Effect_Histogram_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(gather(effect, cols, value), aes(x = value)) + 
  geom_histogram() + 
  facet_wrap(.~cols, scales = 'free_x')

dev.off()

#mean_samples_corr <- cor(means, method = "spearman")

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/DO_per_kg/", cv.threshold, "_perc_CV/n=",rem.threshold,"/Effect_Correlation_Matrix_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

pairs(effect, lower.panel = panel.smooth,upper.panel = panel.cor.spear, gap = 0, cex.labels = 0.5, cex = .75)

#corrplot(effect_corr,type = "upper", tl.col = "black", tl.cex = 1.2, cl.cex = 1,  title = "Effect Correlation")

dev.off()

## Trying to deal with negative values by adding minimum + 1 to all columns 

neg_cols <- c("Effect Size", "Fe_mg_per_L_diff", "SpC Diff.", "Temp Diff.", "pH Diff.", "ATP_picomol_per_g_mean")

silt <- "Percent_Silt_mean"

all_others <- c("Effect Size", "Fe_mg_per_L_diff", "SpC Diff.", "Temp Diff.", "pH Diff.", "ATP_picomol_per_g_mean","Percent_Fine_Sand_mean", "Percent_Med_Sand_mean", "Percent_Coarse_Sand_mean", "Percent_Clay_mean", "average_ssa_mean", "In.Grav.Moi. Diff.", "Fin.Grav.Moi. Diff.")
  
  #c("Effect Size", "Fe_mg_per_L_diff", "SpC Diff.", "Temp Diff.", "pH Diff.", "ATP_picomol_per_g_mean","% Fine Sand", "% Med. Sand", "% Coarse Sand", "% Clay", "SSA", "In.Grav.Moi. Diff.", "Fin.Grav.Moi. Diff.")

log_effect <- effect %>% 
  mutate(across(all_of(neg_cols), ~. + (abs((min(.))))+ 1)) %>% 
  mutate(`% Silt` = log10(`Percent_Silt_mean` + 0.5*min(`Percent_Silt_mean`[`Percent_Silt_mean` != min(`Percent_Silt_mean`)]))) %>% 
  mutate(across(all_of(all_others), ~log10(.))) %>% 
  select(-c(`LostGrav.Moi. Diff.`)) %>% 
  rename_with(~paste0("log_", .), everything())

cube_root <- function(x) sign(x) * (abs(x))^(1/3)
  
cube_root_effect = effect %>% 
  mutate_all(cube_root) %>% 
  rename_all(~ paste0(., "_cube")) 

ggplot(gather(cube_root_effect, cols, value), aes(x = value)) + 
  geom_histogram() + 
  facet_wrap(.~cols, scales = 'free_x')

ggplot(gather(effect, cols, value), aes(x = value)) + 
  geom_histogram() + 
  facet_wrap(.~cols, scales = 'free_x')

cube_root_comp = cube_root_effect %>% 
  rownames_to_column("Sample_Name") %>% 
  left_join(effect_comp, by = "Sample_Name")

log_effect_only = effect %>% 
  mutate(log_effect = log10(`Effect Size`))

ggplot(log_effect_only, aes(x = `Effect Size`, y = log_effect)) +
  geom_point()

ggplot(cube_root_comp, aes(x = `Effect Size`, y = `Effect Size_cube`)) + 
  geom_point()

ggplot(cube_root_comp, aes(x = `Effect Size_cube`))+
  geom_histogram()

cube_root_hist = cube_root_comp %>% 
  column_to_rownames("Sample_Name")




png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/DO_per_kg/", cv.threshold, "_perc_CV/n=", rem.threshold,"/Log_Effect_Histogram_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(gather(log_effect, cols, value), aes(x = value)) + 
  geom_histogram() + 
  facet_wrap(.~cols, scales = 'free_x')

dev.off()

#mean_samples_corr <- cor(means, method = "spearman")

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/DO_per_kg/", cv.threshold, "_perc_CV/n=",rem.threshold,"/Log_Effect_Correlation_Matrix_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

pairs(log_effect, lower.panel = panel.smooth,upper.panel = panel.cor.pear, gap = 0, cex.labels = 0.5, cex = .75)

#corrplot(effect_corr,type = "upper", tl.col = "black", tl.cex = 1.2, cl.cex = 1,  title = "Effect Correlation")

dev.off()




wet_dry_effect <- means %>% 
  rownames_to_column("kit_treat") %>% 
  separate(kit_treat, c("kit", "treat"), remove = FALSE) %>% 
  group_by(kit) %>% 
  filter(n() >= 2) %>% 
  ungroup() %>% 
  group_by(kit) %>% 
  relocate(c(`% Fine Sand`:SSA), .after = `Mean LostGrav.Moi.`) %>% 
  mutate(`Effect Size` = (`Mean Respiration (mg/kg)`[treat == "Wet"])/(`Mean Respiration (mg/kg)`[treat == "Dry"])) %>%
  mutate(`Fe (mg/kg)` = `Mean Fe (mg/kg)`[treat == "Wet"]/`Mean Fe (mg/kg)`[treat == "Dry"]) %>%
  mutate(`SpC` = `Mean SpC`[treat == "Wet"]/`Mean SpC`[treat == "Dry"]) %>%
  mutate(`Temp` = `Mean Temp`[treat == "Wet"]/`Mean Temp`[treat == "Dry"]) %>%
  mutate(`pH` = `Mean pH`[treat == "Wet"]/`Mean pH`[treat == "Dry"]) %>%
  mutate(`In.Grav.Moi.` = `Mean In.Grav.Moi.`[treat == "Wet"]/`Mean In.Grav.Moi.`[treat == "Dry"]) %>%
  mutate(`Fin.Grav.Moi.` = `Mean Fin.Grav.Moi.`[treat == "Wet"]/`Mean Fin.Grav.Moi.`[treat == "Dry"]) %>%
  mutate(`LostGrav.Moi.` = `Mean LostGrav.Moi.`[treat == "Wet"]/`Mean LostGrav.Moi.`[treat == "Dry"]) %>%
  mutate(across(`% Fine Sand`:SSA, mean)) %>% 
  filter(treat != "Wet") %>% 
  select(c(kit, `Effect Size`, `Fe (mg/kg)`, SpC, Temp, pH, In.Grav.Moi., Fin.Grav.Moi., LostGrav.Moi., `% Fine Sand`, `% Med. Sand`, `% Coarse Sand`, `% Clay`, `% Silt`, SSA)) %>% 
  filter(kit != "018") %>% 
  filter(kit != "024") %>% 
  filter(kit != "039")%>% 
  filter(kit != "056") %>% 
  column_to_rownames("kit") #%>% 
# rename(`Effect Size` = `Mean Respiration (mg/kg)`) %>% 
# rename(`Fe (mg/kg) Diff.` = `Mean Fe (mg/kg)`) %>% 
# rename(`SpC Diff.` = `Mean SpC`) %>% 
# rename(`Temp Diff.` = `Mean Temp`) %>% 
# rename(`pH Diff.` = `Mean pH`) %>% 
# rename(`In.Grav.Moi. Diff.` = `Mean In.Grav.Moi.`) %>% 
# rename(`Fin.Grav.Moi. Diff.` = `Mean Fin.Grav.Moi.`) %>%
# rename(`LostGrav.Moi. Diff.` = `Mean LostGrav.Moi.`)
# 

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/DO_per_kg/", cv.threshold, "_perc_CV/n=", rem.threshold,"/Wet_Div_Dry_Effect_Histogram_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(gather(wet_dry_effect, cols, value), aes(x = value)) + 
  geom_histogram() + 
  facet_wrap(.~cols, scales = 'free_x')

dev.off()

#mean_samples_corr <- cor(means, method = "spearman")

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/DO_per_kg/", cv.threshold, "_perc_CV/n=",rem.threshold,"/Wet_Div_Dry_Effect_Correlation_Matrix_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

pairs(wet_dry_effect, lower.panel = panel.smooth,upper.panel = panel.cor.spear, gap = 0, cex.labels = 0.5, cex = .75)

#corrplot(effect_corr,type = "upper", tl.col = "black", tl.cex = 1.2, cl.cex = 1,  title = "Effect Correlation")

dev.off()



### LASSO REGRESSION
#scale because large differences in magnitude between pred. variables

lasso <- all_data %>% 
  drop_na(Fe_mg_per_L) %>% 
  drop_na(average_ssa) %>% 
  drop_na(Initial_Gravimetric_Water) %>% 
  mutate(Respiration_Rate_mg_DO_per_L_per_H = abs(Respiration_Rate_mg_DO_per_L_per_H))

resp <- lasso$Respiration_Rate_mg_DO_per_L_per_H

resp <- scale(resp)

pred <- data.matrix(lasso[, c('Fe_mg_per_L', 'Percent_Fine_Sand', 'Percent_Coarse_Sand', 'Percent_Tot_Sand', 'Percent_Silt', 'Percent_Clay', 'average_ssa', 'SpC', 'Temp', 'pH', 'Initial_Gravimetric_Water', 'Final_Gravimetric_Water', 'Lost_Gravimetric_Water')])

pred <- scale(pred)

cv_model <- cv.glmnet(pred, resp, alpha = 1)

best_lambda <- cv_model$lambda.min
best_lambda

plot(cv_model)

best_model <- glmnet(pred, resp, alpha = 1, lambda = best_lambda)
coef(best_model)

resp_predict <- predict(best_model, s = best_lambda, newx = pred)

sst <- sum((resp - mean(resp))^2)
sse <- sum((resp_predict - resp)^2)

rsq = 1 - sse/sst

rsq

## LASSO ####

## Log LASSO

# lasso <- all_data %>% 
#   drop_na(Fe_mg_per_L) %>% 
#   drop_na(average_ssa) %>% 
#   drop_na(Initial_Gravimetric_Water) %>% 
#   mutate(Respiration_Rate_mg_DO_per_L_per_H = abs(Respiration_Rate_mg_DO_per_L_per_H))

log_lasso <- log_corr_samples %>% 
  rename(log_fe = `Log Fe (mg/L)`) %>% 
  rename(log_atp = `Log ATP (pmol/g)`) %>% 
  rename(log_fine_sand = `Log % Fine Sand`) %>% 
  rename(log_med_sand = `Log % Med. Sand`) %>% 
  rename(log_coarse_sand = `Log % Coarse Sand`) %>% 
  rename(log_clay = `Log % Clay`) %>% 
  rename(log_silt = `Log % Silt`) %>% 
  rename(log_ssa = `Log SSA`) %>% 
  rename(log_spc = `Log SpC`) %>% 
  rename(log_temp = `Log Temp`) %>% 
  rename(log_ph = `Log pH`) %>% 
  rename(log_in_grav_moi = `Log In.Grav.Moi.`) %>% 
  rename(log_fin_grav_moi = `Log Fin.Grav.Moi.`)

log_resp <- log_lasso$`Log Respiration (mg/L)`

#log_resp <- scale(log_resp)

log_pred <- data.matrix(log_lasso[, c("log_fe", "log_atp",  "log_fine_sand", "log_med_sand", "log_coarse_sand", "log_clay","log_silt", "log_ssa", "log_spc", "log_temp", "log_ph","log_in_grav_moi", "log_fin_grav_moi")])

#log_pred <- scale(log_pred)

log_cv_model <- cv.glmnet(log_pred, log_resp, alpha = 1)

log_best_lambda <- log_cv_model$lambda.min
log_best_lambda

plot(log_cv_model)

log_best_model <- glmnet(log_pred, log_resp, alpha = 1, lambda = log_best_lambda)
coef(log_best_model)

log_resp_predict <- predict(log_best_model, s = log_best_lambda, newx = log_pred)

log_sst <- sum((log_resp - mean(log_resp))^2)
log_sse <- sum((log_resp_predict - log_resp)^2)

log_rsq = 1 - log_sse/log_sst

log_rsq

## Log FS, SSA, Fe, Resp, Moi

log_pred_co <- data.matrix(log_lasso[, c("log_fe", "log_atp",  "log_fine_sand", "log_ssa", "log_spc", "log_temp", "log_ph", "log_fin_grav_moi")])

#log_pred <- scale(log_pred)

log_cv_model_co <- cv.glmnet(log_pred_co, log_resp, alpha = 1)

log_best_lambda_co <- log_cv_model_co$lambda.min
log_best_lambda_co

plot(log_cv_model_co)

log_best_model_co <- glmnet(log_pred_co, log_resp, alpha = 1, lambda = log_best_lambda_co)
coef(log_best_model_co)

log_resp_predict_co <- predict(log_best_model_co, s = log_best_lambda_co, newx = log_pred_co)

log_sst_co <- sum((log_resp - mean(log_resp))^2)
log_sse_co <- sum((log_resp_predict_co - log_resp)^2)

log_rsq_co = 1 - log_sse_co/log_sst_co

log_rsq_co

## Log LASSO Dry Samples

log_lasso_dry <- log_lasso %>% 
  rownames_to_column("Sample_Name") %>% 
  filter(grepl("D", Sample_Name)) %>%
  dplyr::select(c(`Log Respiration (mg/L)`, log_fe, log_atp, log_fine_sand, log_ssa, log_spc, log_temp, log_ph, log_fin_grav_moi))

log_resp_dry <- log_lasso_dry$`Log Respiration (mg/L)`

log_pred_dry <- data.matrix(log_lasso_dry[, c('log_fe', 'log_atp', 'log_fine_sand', 'log_ssa', 'log_spc', 'log_temp', 'log_ph', 'log_fin_grav_moi')])

#log_pred <- scale(log_pred)

log_cv_model_dry <- cv.glmnet(log_pred_dry, log_resp_dry, alpha = 1)

log_best_lambda_dry <- log_cv_model_dry$lambda.min
log_best_lambda_dry

plot(log_cv_model_dry)

log_best_model_dry <- glmnet(log_pred_dry, log_resp_dry, alpha = 1, lambda = log_best_lambda_dry)
coef(log_best_model_dry)

log_resp_predict_dry <- predict(log_best_model_dry, s = log_best_lambda_dry, newx = log_pred_dry)

log_sst_dry <- sum((log_resp_dry - mean(log_resp_dry))^2)
log_sse_dry <- sum((log_resp_predict_dry - log_resp_dry)^2)

log_rsq_dry = 1 - log_sse_dry/log_sst_dry

log_rsq_dry

## Log Wet LASSO
log_lasso_wet <- log_lasso %>% 
  rownames_to_column("Sample_Name") %>% 
  filter(grepl("W", Sample_Name)) %>%
  dplyr::select(c(`Log Respiration (mg/L)`, log_fe, log_atp, log_fine_sand, log_ssa, log_spc, log_temp, log_ph, log_fin_grav_moi))

log_resp_wet <- log_lasso_wet$`Log Respiration (mg/L)`

log_pred_wet <- data.matrix(log_lasso_wet[, c('log_fe', 'log_atp', 'log_fine_sand', 'log_ssa', 'log_spc', 'log_temp', 'log_ph', 'log_fin_grav_moi')])

#log_pred <- scale(log_pred)

log_cv_model_wet <- cv.glmnet(log_pred_wet, log_resp_wet, alpha = 1)

log_best_lambda_wet <- log_cv_model_wet$lambda.min
log_best_lambda_wet

plot(log_cv_model_wet)

log_best_model_wet <- glmnet(log_pred_wet, log_resp_wet, alpha = 1, lambda = log_best_lambda_wet)
coef(log_best_model_wet)

log_resp_predict_wet <- predict(log_best_model_wet, s = log_best_lambda_wet, newx = log_pred_wet)

log_sst_wet <- sum((log_resp_wet - mean(log_resp_wet))^2)
log_sse_wet <- sum((log_resp_predict_wet - log_resp_wet)^2)

log_rsq_wet = 1 - log_sse_wet/log_sst_wet

log_rsq_wet


## Effect LASSO 
effect_lasso <- log_effect

eff <- effect_lasso$`log_Effect Size`

#eff <- scale(eff)

pred <- data.matrix(effect_lasso[, c("log_Fe_mg_per_L_diff", "log_ATP_picomol_per_g_mean", "log_SpC Diff.", "log_Temp Diff.", "log_pH Diff.", "log_Fin.Grav.Moi. Diff.", "log_Percent_Fine_Sand_mean", "log_Percent_Med_Sand_mean", "log_Percent_Coarse_Sand_mean", "log_Percent_Clay_mean", "log_% Silt", "log_average_ssa_mean")]) 
  
#pred <- scale(pred)

cv_model <- cv.glmnet(pred, eff, alpha = 1)

best_lambda <- cv_model$lambda.min
best_lambda

plot(cv_model)

best_model <- glmnet(pred, eff, alpha = 1, lambda = best_lambda)
coef(best_model)

eff_predict <- predict(best_model, s = best_lambda, newx = pred)

sst <- sum((eff - mean(eff))^2)
sse <- sum((eff_predict - eff)^2)

rsq = 1 - sse/sst

rsq

eff <- effect_lasso$`log_Effect Size`

#eff <- scale(eff)

pred <- data.matrix(effect_lasso[, c("log_Fe_mg_per_L_diff", "log_ATP_picomol_per_g_mean", "log_SpC Diff.", "log_Temp Diff.", "log_pH Diff.", "log_Fin.Grav.Moi. Diff.", "log_Percent_Fine_Sand_mean", "log_average_ssa_mean")]) 

#pred <- scale(pred)

cv_model <- cv.glmnet(pred, eff, alpha = 1)

best_lambda <- cv_model$lambda.min
best_lambda

plot(cv_model)

best_model <- glmnet(pred, eff, alpha = 1, lambda = best_lambda)
coef(best_model)

eff_predict <- predict(best_model, s = best_lambda, newx = pred)

sst <- sum((eff - mean(eff))^2)
sse <- sum((eff_predict - eff)^2)

rsq = 1 - sse/sst

rsq


corr_samples_treat_real <- corr_samples %>%
  rownames_to_column("Sample_Name") %>% 
  mutate(Treat = case_when(grepl("W",Sample_Name)~"Wet",
                           grepl("D", Sample_Name) ~"Dry")) %>% 
  filter(`Respiration (mg/L)` < 90)

ggplot(corr_samples_treat_real, aes(x = `Respiration (mg/L)`, y = ATP_picomol_per_g, color = Treat)) + 
  geom_point()

ggplot(corr_samples_treat_real, aes(x = `Respiration (mg/L)`, y = Fe_mg_per_L, color = Treat)) + 
  geom_point()

ggplot(corr_samples_treat, aes(x = `Respiration (mg/L)`, y = Fe_mg_per_L, color = Treat)) + 
  geom_point()

dry_samples = corr_samples_treat %>% 
  filter()

ggplot(effect, aes(x = `Effect Size`, y = ATP_picomol_per_g_mean)) + 
  geom_point() + 
  ylab("ATP Diff.")

ggplot(effect, aes(x = `Effect Size`, y = Percent_Fine_Sand_mean)) + 
  geom_point() + 
  ylab("% Fine Sand")

ggplot(effect, aes(x = `Effect Size`, y = Fe_mg_per_L_diff)) + 
  geom_point() + 
  ylab("Fe Diff.")

ggplot(effect, aes(x = `Effect Size`, y = `Fin.Grav.Moi. Diff.`)) + 
  geom_point() + 
  ylab("Final Grav Moi. Diff.")

