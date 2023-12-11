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
  dplyr::select(c(Sample_Name, Respiration_Rate_mg_DO_per_L_per_H,Respiration_Rate_mg_DO_per_kg_per_H, Mean_Slope_All_L, Mean_Slope_All_kg, cv_before_removal_L, cv_before_removal_kg, kit_treat)) 

slope.outliers$flag <- NA

slope.final <- as.data.frame(matrix(NA, ncol = 11, nrow =1))

colnames(slope.final) = c("slope.temp","Sample_Name", "kit_treat", "Respiration_Rate_mg_DO_per_kg_per_H", "Mean_Slope_All_L", "Mean_Slope_All_kg", "cv_before_removal_L", "cv_before_removal_kg", "cv_after_removal", "Mean_Slope_Removed","flag")

unique.samples = unique(slope.outliers$kit_treat)

#try 0, 10, 30, 50, 100, export histograms of removals, effect sizes
cv.threshold = 250

#try keeping n = 3,4 samples
rem.threshold = 4

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

cv.threshold = "No_Removals"
rem.threshold = 3

#Write data frame with removed respiration rates

write.csv(slope.final.flag, paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/", cv.threshold, "_perc_CV/n=", rem.threshold,"/CV_mg_L_",cv.threshold,"percent_Removed_Respiration_Rates",Sys.Date(),".csv"), row.names = F)

## All Samples ####

##Merge removed samples with other data
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

## Make histogram of all data
png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/", cv.threshold, "_perc_CV/n=", rem.threshold,"/All_Samples_Histogram_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

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
png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/", cv.threshold, "_perc_CV/n=",rem.threshold,"/All_Samples_Correlation_Matrix_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

pairs(corr_samples, lower.panel = panel.smooth,upper.panel = panel.cor.spear, gap = 0, cex.labels = 0.5, cex = .75)

#corrplot(all_samples_corr,type = "upper", tl.col = "black", tl.cex = 1.2, cl.cex = 1,  title = "All Samples Correlation")

dev.off()

## Log All Samples ####

log_corr_samples <- corr_samples %>% 
  mutate(`Log Respiration (mg/L)` = log10(`Respiration (mg/L)` + 0.5*min(`Respiration (mg/L)`[`Respiration (mg/L)` != min(`Respiration (mg/L)`)]))) %>% 
  mutate(`Log Fe (mg/L)` = log10(`Fe (mg/L)` + 0.5*min(`Fe (mg/L)`[`Fe (mg/L)` != min(`Fe (mg/L)`)]))) %>% 
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
   dplyr::select(c(`Log Respiration (mg/L)`, `Log Fe (mg/L)`, `Log % Fine Sand`, `Log % Med. Sand`, `Log % Coarse Sand`, `Log % Clay`, `Log % Silt`, `Log SSA`, `Log SpC`, `Log Temp`, `Log pH`, `Log In.Grav.Moi.`, `Log Fin.Grav.Moi.`, `Log LostGrav.Moi.`))

# histogram of log transformed samples (log + 1/2 of minimum to get at 0 values)

#question - do this across the board even if no 0's?

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/", cv.threshold, "_perc_CV/n=", rem.threshold,"/Log_All_Samples_Histogram_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(gather(log_corr_samples, cols, value), aes(x = value)) + 
  geom_histogram() + 
  facet_wrap(.~cols, scales = 'free_x')

dev.off()

#log_samples_corr <- cor(log_corr_samples, method = "pearson")

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/", cv.threshold, "_perc_CV/n=",rem.threshold,"/Log_All_Samples_Correlation_Matrix_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

pairs(corr_samples, lower.panel = panel.smooth,upper.panel = panel.cor.pear, gap = 0, cex.labels = 0.5, cex = .75)

#corrplot(log_samples_corr,type = "upper", tl.col = "black", tl.cex = 1.2, cl.cex = 1,  title = "Log All Samples Correlation")

dev.off()

## Wet/Dry All Samples ####

all_samples_dry <- corr_samples %>% 
  rownames_to_column("Sample_Name") %>% 
  filter(grepl("D", Sample_Name)) %>% 
  column_to_rownames("Sample_Name")

## Make histogram of all data
png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/", cv.threshold, "_perc_CV/n=", rem.threshold,"/All_Dry_Samples_Histogram_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(gather(all_samples_dry, cols, value), aes(x = value)) + 
  geom_histogram() + 
  facet_wrap(.~cols, scales = 'free_x')

dev.off()

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/", cv.threshold, "_perc_CV/n=",rem.threshold,"/All_Dry_Samples_Correlation_Matrix_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

pairs(all_samples_dry, lower.panel = panel.smooth,upper.panel = panel.cor.spear, gap = 0, cex.labels = 0.5, cex = .75)

#corrplot(all_samples_corr,type = "upper", tl.col = "black", tl.cex = 1.2, cl.cex = 1,  title = "All Samples Correlation")

dev.off()

log_all_dry <- log_corr_samples %>% 
  rownames_to_column("Sample_Name") %>% 
  filter(grepl("D", Sample_Name)) %>% 
  column_to_rownames("Sample_Name")
  
png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/", cv.threshold, "_perc_CV/n=", rem.threshold,"/Log_All_Dry_Samples_Histogram_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(gather(log_all_dry, cols, value), aes(x = value)) + 
  geom_histogram() + 
  facet_wrap(.~cols, scales = 'free_x')

dev.off()

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/", cv.threshold, "_perc_CV/n=",rem.threshold,"/Log_All_Dry_Samples_Correlation_Matrix_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

pairs(log_all_dry, lower.panel = panel.smooth,upper.panel = panel.cor.pear, gap = 0, cex.labels = 0.5, cex = .75)

dev.off()

all_samples_wet <-  corr_samples %>% 
  rownames_to_column("Sample_Name") %>% 
  filter(grepl("W", Sample_Name)) %>% 
  column_to_rownames("Sample_Name")

## Make histogram of all data
png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/", cv.threshold, "_perc_CV/n=", rem.threshold,"/All_Wet_Samples_Histogram_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(gather(all_samples_wet, cols, value), aes(x = value)) + 
  geom_histogram() + 
  facet_wrap(.~cols, scales = 'free_x')

dev.off()

#all_samples_dry_corr <- cor(all_samples_dry,method = "spearman")

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/", cv.threshold, "_perc_CV/n=",rem.threshold,"/All_Wet_Samples_Correlation_Matrix_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

pairs(all_samples_wet, lower.panel = panel.smooth,upper.panel = panel.cor.spear, gap = 0, cex.labels = 0.5, cex = .75)

#corrplot(all_samples_corr,type = "upper", tl.col = "black", tl.cex = 1.2, cl.cex = 1,  title = "All Samples Correlation")

dev.off()

log_all_wet <- log_corr_samples %>% 
  rownames_to_column("Sample_Name") %>% 
  filter(grepl("W", Sample_Name)) %>% 
  dplyr::select(-c(`Log LostGrav.Moi.`)) %>% 
  column_to_rownames("Sample_Name")

## Make histogram of all data
png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/", cv.threshold, "_perc_CV/n=", rem.threshold,"/Log_All_Wet_Samples_Histogram_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(gather(log_all_wet, cols, value), aes(x = value)) + 
  geom_histogram() + 
  facet_wrap(.~cols, scales = 'free_x')

dev.off()

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/", cv.threshold, "_perc_CV/n=",rem.threshold,"/Log_All_Wet_Samples_Correlation_Matrix_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

pairs(log_all_wet, lower.panel = panel.smooth,upper.panel = panel.cor.pear, gap = 0, cex.labels = 0.5, cex = .75)

dev.off()

## Means ####

means <- corr_samples %>% 
  rownames_to_column("Sample_Name") %>% 
  separate(Sample_Name, c("EC", "kit", "Treat"), remove = FALSE) %>% 
  mutate(Treat = case_when(grepl("W",Sample_Name)~"Wet",
                           grepl("D", Sample_Name) ~"Dry")) %>% 
  group_by(kit, Treat) %>% 
  summarise(across(where(is.numeric), list(mean = mean), na.rm = TRUE))%>% 
  unite(kit_treat, c("kit", "Treat")) %>% 
  rename(`Mean Respiration (mg/L)` = `Respiration (mg/L)_mean`) %>% 
  rename(`Mean Fe (mg/L)` = `Fe (mg/L)_mean`) %>% 
  rename(`% Fine Sand` = `% Fine Sand_mean`) %>% 
  rename(`% Med. Sand` = `% Med. Sand_mean`) %>% 
  rename(`% Coarse Sand` = `% Coarse Sand_mean`) %>% 
  rename(`% Clay` = `% Clay_mean`) %>% 
  rename(`% Silt` = `% Silt_mean`) %>% 
  rename(SSA = SSA_mean) %>% 
  rename(`Mean pH` = pH_mean) %>% 
  rename(`Mean SpC` = SpC_mean) %>% 
  rename(`Mean Temp` = Temp_mean) %>% 
  rename(`Mean In.Grav.Moi.` = In.Grav.Moi._mean) %>% 
  rename(`Mean Fin.Grav.Moi.` = Fin.Grav.Moi._mean) %>% 
  rename(`Mean LostGrav.Moi.` = LostGrav.Moi._mean) %>% 
  column_to_rownames("kit_treat")

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/", cv.threshold, "_perc_CV/n=", rem.threshold,"/Mean_Samples_Histogram_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(gather(means, cols, value), aes(x = value)) + 
  geom_histogram() + 
  facet_wrap(.~cols, scales = 'free_x')

dev.off()

#mean_samples_corr <- cor(means, method = "spearman")

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/", cv.threshold, "_perc_CV/n=",rem.threshold,"/Mean_Samples_Correlation_Matrix_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

pairs(means, lower.panel = panel.smooth,upper.panel = panel.cor.spear, gap = 0, cex.labels = 0.5, cex = .75)

#corrplot(all_samples_corr,type = "upper", tl.col = "black", tl.cex = 1.2, cl.cex = 1,  title = "All Samples Correlation")

dev.off()

## Log Mean Samples ####
log_means <- means %>% 
  mutate(`Log Mean Respiration (mg/L)` = log10(`Mean Respiration (mg/L)` + 0.5*min(`Mean Respiration (mg/L)`[`Mean Respiration (mg/L)` != min(`Mean Respiration (mg/L)`)]))) %>% 
  mutate(`Log Mean Fe (mg/L)` = log10(`Mean Fe (mg/L)` + 0.5*min(`Mean Fe (mg/L)`[`Mean Fe (mg/L)` != min(`Mean Fe (mg/L)`)]))) %>% 
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
  dplyr::select(c(`Log Mean Respiration (mg/L)`, `Log Mean Fe (mg/L)`, `Log % Fine Sand`, `Log % Med. Sand`, `Log % Coarse Sand`, `Log % Clay`, `Log % Silt`, `Log SSA`, `Log Mean SpC`, `Log Mean Temp`, `Log Mean pH`, `Log Mean In.Grav.Moi.`, `Log Mean Fin.Grav.Moi.`, `Log Mean LostGrav.Moi.`))

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/", cv.threshold, "_perc_CV/n=", rem.threshold,"/Log_Mean_Samples_Histogram_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(gather(log_means, cols, value), aes(x = value)) + 
  geom_histogram() + 
  facet_wrap(.~cols, scales = 'free_x')

dev.off()

#mean_samples_corr <- cor(means, method = "spearman")

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/", cv.threshold, "_perc_CV/n=",rem.threshold,"/Log_Mean_Samples_Correlation_Matrix_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

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
  
png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/", cv.threshold, "_perc_CV/n=", rem.threshold,"/Wet_Mean_Samples_Histogram_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(gather(wet_means, cols, value), aes(x = value)) + 
  geom_histogram() + 
  facet_wrap(.~cols, scales = 'free_x')

dev.off()

#mean_samples_corr <- cor(means, method = "spearman")

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/", cv.threshold, "_perc_CV/n=",rem.threshold,"/Wet_Mean_Samples_Correlation_Matrix_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

pairs(wet_means, lower.panel = panel.smooth,upper.panel = panel.cor.spear, gap = 0, cex.labels = 0.5, cex = .75)

#corrplot(all_samples_corr,type = "upper", tl.col = "black", tl.cex = 1.2, cl.cex = 1,  title = "All Samples Correlation")

dev.off()


dry_means <- means %>% 
  rownames_to_column("kit_treat") %>% 
  filter(grepl("D", kit_treat)) %>% 
  column_to_rownames("kit_treat")

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/", cv.threshold, "_perc_CV/n=", rem.threshold,"/Dry_Mean_Samples_Histogram_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(gather(dry_means, cols, value), aes(x = value)) + 
  geom_histogram() + 
  facet_wrap(.~cols, scales = 'free_x')

dev.off()

#mean_samples_corr <- cor(means, method = "spearman")

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/", cv.threshold, "_perc_CV/n=",rem.threshold,"/Dry_Mean_Samples_Correlation_Matrix_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

pairs(dry_means, lower.panel = panel.smooth,upper.panel = panel.cor.spear, gap = 0, cex.labels = 0.5, cex = .75)

#corrplot(all_samples_corr,type = "upper", tl.col = "black", tl.cex = 1.2, cl.cex = 1,  title = "All Samples Correlation")

dev.off()

## Log Mean Samples

log_wet_means <- log_means %>% 
  rownames_to_column("kit_treat") %>% 
  filter(grepl("W", kit_treat)) %>% 
  dplyr::select(-c(`Log Mean LostGrav.Moi.`, `Log Mean In.Grav.Moi.`, `Log Mean Fin.Grav.Moi.`)) %>% 
  column_to_rownames("kit_treat")

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/", cv.threshold, "_perc_CV/n=", rem.threshold,"/Log_Wet_Mean_Samples_Histogram_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(gather(log_wet_means, cols, value), aes(x = value)) + 
  geom_histogram() + 
  facet_wrap(.~cols, scales = 'free_x')

dev.off()

#mean_samples_corr <- cor(means, method = "spearman")

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/", cv.threshold, "_perc_CV/n=",rem.threshold,"/Log_Wet_Mean_Samples_Correlation_Matrix_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

pairs(log_wet_means, lower.panel = panel.smooth,upper.panel = panel.cor.pear, gap = 0, cex.labels = 0.5, cex = .75)

#corrplot(all_samples_corr,type = "upper", tl.col = "black", tl.cex = 1.2, cl.cex = 1,  title = "All Samples Correlation")

dev.off()


log_dry_means <- log_means %>% 
  rownames_to_column("kit_treat") %>% 
  filter(grepl("D", kit_treat)) %>% 
  column_to_rownames("kit_treat")
  
png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/", cv.threshold, "_perc_CV/n=", rem.threshold,"/Log_Dry_Mean_Samples_Histogram_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(gather(log_dry_means, cols, value), aes(x = value)) + 
  geom_histogram() + 
  facet_wrap(.~cols, scales = 'free_x')

dev.off()

#mean_samples_corr <- cor(means, method = "spearman")

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/", cv.threshold, "_perc_CV/n=",rem.threshold,"/Log_Dry_Mean_Samples_Correlation_Matrix_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

pairs(log_dry_means, lower.panel = panel.smooth,upper.panel = panel.cor.pear, gap = 0, cex.labels = 0.5, cex = .75)

#corrplot(all_samples_corr,type = "upper", tl.col = "black", tl.cex = 1.2, cl.cex = 1,  title = "All Samples Correlation")

dev.off()

#####


effect <- means %>% 
  rownames_to_column("kit_treat") %>% 
  separate(kit_treat, c("kit", "treat"), remove = FALSE) %>% 
  group_by(kit) %>% 
  filter(n() >= 2) %>% 
  ungroup() %>% 
  group_by(kit) %>% 
  relocate(c(`% Fine Sand`:SSA), .after = `Mean LostGrav.Moi.`) %>% 
  summarise(
    across(`Mean Respiration (mg/L)`:`Mean LostGrav.Moi.`, diff), 
    across(`% Fine Sand`:SSA, mean))%>% 
  column_to_rownames("kit") %>% 
  rename(`Effect Size` = `Mean Respiration (mg/L)`) %>% 
  rename(`Fe (mg/L) Diff.` = `Mean Fe (mg/L)`) %>% 
  rename(`SpC Diff.` = `Mean SpC`) %>% 
  rename(`Temp Diff.` = `Mean Temp`) %>% 
  rename(`pH Diff.` = `Mean pH`) %>% 
  rename(`In.Grav.Moi. Diff.` = `Mean In.Grav.Moi.`) %>% 
  rename(`Fin.Grav.Moi. Diff.` = `Mean Fin.Grav.Moi.`) %>%
  rename(`LostGrav.Moi. Diff.` = `Mean LostGrav.Moi.`)

#effect_corr <- cor(effect, method = "spearman")

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/", cv.threshold, "_perc_CV/n=", rem.threshold,"/Effect_Histogram_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(gather(effect, cols, value), aes(x = value)) + 
  geom_histogram() + 
  facet_wrap(.~cols, scales = 'free_x')

dev.off()

#mean_samples_corr <- cor(means, method = "spearman")

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/", cv.threshold, "_perc_CV/n=",rem.threshold,"/Effect_Correlation_Matrix_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

pairs(effect, lower.panel = panel.smooth,upper.panel = panel.cor.spear, gap = 0, cex.labels = 0.5, cex = .75)

#corrplot(effect_corr,type = "upper", tl.col = "black", tl.cex = 1.2, cl.cex = 1,  title = "Effect Correlation")

dev.off()

log_effect <- log_means %>% 
  rownames_to_column("kit_treat") %>% 
  separate(kit_treat, c("kit", "treat")) %>% 
  group_by(kit) %>%   filter(n() >= 2) %>% 
  ungroup() %>% 
  group_by(kit) %>% 
  relocate(c(`Log % Fine Sand`:`Log SSA`), .after = `Log Mean LostGrav.Moi.`) %>% 
  summarise(
    across(`Log Mean Respiration (mg/L)`:`Log Mean LostGrav.Moi.`, diff), 
    across(`Log % Fine Sand`:`Log SSA`, mean))%>% 
  column_to_rownames("kit") %>% 
  rename(`Log Effect Size` = `Log Mean Respiration (mg/L)`) %>% 
  rename(`Log Fe (mg/L) Diff.` = `Log Mean Fe (mg/L)`) %>% 
  rename(`Log SpC Diff.` = `Log Mean SpC`) %>% 
  rename(`Log Temp Diff.` = `Log Mean Temp`) %>% 
  rename(`Log pH Diff.` = `Log Mean pH`) %>% 
  rename(`Log In.Grav.Moi. Diff.` = `Log Mean In.Grav.Moi.`) %>% 
  rename(`Log Fin.Grav.Moi. Diff.` = `Log Mean Fin.Grav.Moi.`) %>%
  dplyr::select(-c(`Log Mean LostGrav.Moi.`))

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/", cv.threshold, "_perc_CV/n=", rem.threshold,"/Log_Effect_Histogram_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(gather(log_effect, cols, value), aes(x = value)) + 
  geom_histogram() + 
  facet_wrap(.~cols, scales = 'free_x')

dev.off()

#mean_samples_corr <- cor(means, method = "spearman")

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Data/Effect Size Sensitivity Analysis/", cv.threshold, "_perc_CV/n=",rem.threshold,"/Log_Effect_Correlation_Matrix_CV_", cv.threshold, "percent_Removed.png"), width = 8, height = 8, units = "in", res = 300)

pairs(log_effect, lower.panel = panel.smooth,upper.panel = panel.cor.pear, gap = 0, cex.labels = 0.5, cex = .75)

#corrplot(effect_corr,type = "upper", tl.col = "black", tl.cex = 1.2, cl.cex = 1,  title = "Effect Correlation")

dev.off()

