#### Load Library ####

library(lubridate);library(writexl);library(raster);library(tidyverse);library(devtools)
library(patchwork)

#### Load data #####
rm(list=ls());graphics.off()

# Set working directory to data file

pnnl.user = 'laan208'

input.path <- paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/rates/Plots")

setwd(input.path)

import_data = function(input.path){
  
  # pull a list of file names in the target folder with the target pattern
  # then read all files and combine
  
  filePaths <- list.files(path = input.path, recursive = T, pattern = "Range_2023", full.names = TRUE)
  
  # dat <- 
  do.call(rbind, lapply(filePaths, function(input.path){
    # then add a new column `source` to denote the file name
    df <- read.csv(input.path)
    df[["source_file"]] <- rep(input.path, nrow(df)) # add a column for the source file
    
    df %>%
      mutate(source_file = str_remove_all(source_file, paste0(input.path, "/")))
  }
  ))
}
data = import_data(input.path)

data_long = 
  data %>% 
  mutate(source_file = str_remove_all(source_file, ".*Rates_"),source_file = str_remove_all(source_file, "_2023.*")) %>% 
  drop_na(Sample_Name)


data_clean <- data_long %>% 
  separate(Sample_Name, c("ECA", "kit", "rep"), sep = "_", remove = FALSE) %>% 
  mutate(Treat = case_when(grepl("W",rep)~"Wet",
                           grepl("D", rep) ~"Dry")) %>% 
  relocate(Treat, .after = rep) %>% 
  dplyr::select(c(Sample_Name,ECA,kit,rep,slope_of_the_regression,rate_mg_per_L_per_min,rate_mg_per_L_per_h,Treat, source_file)) %>% 
  mutate(slope_of_the_regression = if_else(slope_of_the_regression>0,0,slope_of_the_regression)) %>%
  mutate(rate_mg = if_else(slope_of_the_regression>=0,0,rate_mg_per_L_per_min)) %>%  
  group_by(kit, Treat, source_file) %>% 
  mutate(Mean_Slope = mean(rate_mg))

data_clean$source_file_ordered = factor(data_clean$source_file, levels=c('0_Range','0.5_Range','1_Range','1.4_Range', '1.5_Range', '2_Range', '5_Range', '10_Range'))


wet_no <- ggplot(subset(data_clean, Treat %in% "Wet"), aes(x = rate_mg)) +
  geom_histogram()+
  facet_grid(~source_file_ordered)+
  ggtitle("Wet Rates, No Removals")+
  theme(strip.text = element_text(
    size = 4))

dry_no <- ggplot(subset(data_clean, Treat %in% "Dry"), aes(x = rate_mg)) +
  geom_histogram()+
  facet_grid(~source_file_ordered)+
  ggtitle("Dry Rates, No Removals")+
  theme(strip.text = element_text(
    size = 4))

data_clean_effect <- data_clean %>% 
  distinct(kit, Treat, source_file, .keep_all = TRUE) %>% 
  filter(!grepl(011, kit)) %>% 
  filter(!grepl(012, kit)) %>% 
  group_by(kit, source_file) %>% 
  mutate(effect = (Mean_Slope[Treat == "Wet"] - Mean_Slope[Treat == "Dry"])) %>% 
  mutate(log_effect = log10(abs(effect+1))) %>% 
  distinct(kit, .keep_all = TRUE)
 
effect_no <- ggplot(data_clean_effect, aes(x = effect)) + 
  geom_histogram() +
  facet_grid(~source_file_ordered)+
  ggtitle("Effect Size, No Removals")+
  theme(strip.text = element_text(
    size = 4))

#choosing best 4 out of 5 samples to keep using dist matrix


slope.new <- data_long %>% 
  separate(Sample_Name, c("ECA", "kit", "rep"), sep = "_", remove = FALSE) %>% 
  mutate(Treat = case_when(grepl("W",rep)~"Wet",
                           grepl("D", rep) ~"Dry")) %>% 
  relocate(Treat, .after = rep) %>% 
  dplyr::select(c(Sample_Name,ECA,kit,rep,slope_of_the_regression,rate_mg_per_L_per_min,rate_mg_per_L_per_h,Treat, source_file)) %>% 
  group_by(kit, Treat, source_file) %>%
  mutate(Mean_Slope = mean(slope_of_the_regression)) %>% 
  mutate(cv_before_removal = (sd(slope_of_the_regression)/mean(slope_of_the_regression))*100) %>% 
  ungroup() %>% 
  unite(kit_treat_sat, c("kit", "Treat", "source_file"),  sep = "_", remove = FALSE)
  
slope.new$flag <- NA

slope.final <- as.data.frame(matrix(NA, ncol = 15, nrow =1))

colnames(slope.final) = c("slope.temp","Sample_Name", "ECA", "kit_treat_sat", "kit", "rep", "rate_mg_per_L_per_min","rate_mg_per_L_per_h", "Treat", "source_file", "Mean_Slope","cv_before_removal", "cv_after_removal", "Mean_Slope_Removed","flag")


##if more than 4 samples and CV > 10%, then remove 1 sample

unique.samples = unique(slope.new$kit_treat_sat)

for (i in 1:length(unique.samples)) {
  
  data_subset = subset(slope.new, slope.new$kit_treat_sat == unique.samples[i])
  
  slope.temp = as.numeric(data_subset$slope_of_the_regression)
  
  slope.temp.sd <- sd(slope.temp)
  slope.temp.mean <- mean(slope.temp)
  CV = abs((slope.temp.sd/slope.temp.mean)*100)
  
  #looping to get 4 out of 5 best samples
  for (sample.reduction in 1:5)  {
    
    if (slope.temp.mean == 0) {
      
      CV = 0
      
    }
    
    else if (length(slope.temp) > 4 & CV >= 10) {
      
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
    
    # if(CV > 10 ) {
    #   
    #   slope.new$Slope_Removed_Mean[which(slope.new$kit_treat == unique.samples[i])] = "Samples too Variable"
    #   
    # }
    
    # else {
    
    
    slope.combined <- as.data.frame(slope.temp)
    
    slope.removed <- merge(slope.combined, data_subset, by.x = "slope.temp", by.y = "slope_of_the_regression", all.x = TRUE)
    
    slope.removed <- slope.removed[!duplicated(slope.removed$Sample_Name), ]
    
    slope.removed$cv_after_removal = as.numeric(abs((sd(slope.temp)/mean(slope.temp))*100))
    
    slope.removed$Mean_Slope_Removed = as.numeric(mean(slope.temp))
    
    #slope.new$cv_after_removal[which(slope.new$kit_treat == unique.samples[i])] = as.numeric(abs((sd(slope.temp)/mean(slope.temp))*100))
    
    
    # slope.new$Slope_Removed_Mean[which(slope.new$kit_treat == unique.samples[i])] = as.numeric(mean(slope.temp))
    
    #  }
    
  }
  
  slope.final = rbind(slope.removed, slope.final)
  
  rm('slope.temp')
}


## This data frame has removed samples 
slope.final$flag <- ifelse(slope.final$cv_before_removal < slope.final$cv_after_removal, "Issue in dropping sample", NA)

slope.final <- rename(slope.final, "slope_of_the_regression" = "slope.temp")

slope.final$rem <- abs(slope.final$slope_of_the_regression) - slope.final$rate_mg_per_L_per_min

data_removal_clean <- slope.final %>% 
  mutate(slope_of_the_regression = if_else(slope_of_the_regression>0,0,slope_of_the_regression)) %>%
  mutate(rate_mg = if_else(slope_of_the_regression>=0,0,rate_mg_per_L_per_min)) %>%  
  group_by(kit, Treat, source_file) %>% 
  mutate(Mean_Rate_Removed = mean(rate_mg))

data_removal_clean$source_file_ordered = factor(data_removal_clean$source_file, levels=c('0_Range','0.5_Range','1_Range','1.4_Range', '1.5_Range', '2_Range', '5_Range', '10_Range'))

wet_rem <- ggplot(subset(data_removal_clean, Treat %in% "Wet"), aes(x = rate_mg)) +
  geom_histogram()+
  facet_grid(~source_file_ordered)+
  ggtitle("Wet Rates, Removals")+
  theme(strip.text = element_text(
    size = 4))

dry_rem <- ggplot(subset(data_removal_clean, Treat %in% "Dry"), aes(x = rate_mg)) +
  geom_histogram()+
  facet_grid(~source_file_ordered)+
  ggtitle("Dry Rates, Removals")+
  theme(strip.text = element_text(
    size = 4))

data_removal_effect <- data_removal_clean %>% 
  distinct(kit, Treat, source_file, .keep_all = TRUE) %>% 
  filter(!grepl(011, kit)) %>% 
  filter(!grepl(012, kit)) %>%
  group_by(kit, source_file) %>% 
  mutate(effect = (Mean_Rate_Removed[Treat == "Wet"] - Mean_Rate_Removed[Treat == "Dry"])) %>% 
  mutate(log_effect = log10(abs(effect+1))) %>% 
  distinct(kit, .keep_all = TRUE) %>% 
  drop_na(Sample_Name)

effect_rem <- ggplot(data_removal_effect, aes(x = effect)) + 
  geom_histogram() +
  facet_grid(~source_file_ordered)+
  ggtitle("Effect Size, Removals")+
  theme(strip.text = element_text(
    size = 4))


png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/", as.character(Sys.Date()),"_Range_Sensitivity_Analysis_Histograms.png"), width = 8, height = 8, units = "in", res = 300)

all <- (wet_no + wet_rem + dry_no + dry_rem + effect_no + effect_rem) +
  plot_layout(widths = c(3,3))

print(all)

dev.off()
