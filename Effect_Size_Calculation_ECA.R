###### Load Library ######

library(lubridate);library(writexl);library(raster);library(tidyverse);library(devtools)

##### Load data ######
rm(list=ls());graphics.off()

# Set working directory to data file
#Example:
pnnl.user = 'laan208'

input.path <- paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/rates")

setwd(input.path)

path <- ("Plots")

#change date to most recent respiration rate csv

#read in all files, remove csv that are not rate data files, bind all results files together, 
import_data = function(path){
  
  # pull a list of file names in the target folder with the target pattern
  # then read all files and combine
  
  files <- list.files(path, pattern = ".csv", recursive = T, full.names = T)
  
  all <- files[grep("2023-04-06", files)]
  
  data <- lapply(all, read.table, sep = ",", header = T)
  
  data <- do.call(rbind, data)
}

data = import_data(path)


#separate names, add Treatment (wet or dry), replace NAs with blanks

all.data <- data %>% 
  separate(Sample_ID, c("ECA", "kit", "rep"), sep = "_", remove = FALSE) %>% 
  mutate(Treat = case_when(
    endsWith(rep, "W1") ~ "Wet",
    endsWith(rep, "W2") ~ "Wet",
    endsWith(rep, "W3") ~ "Wet",
    endsWith(rep, "W4") ~ "Wet", 
    endsWith(rep, "W5") ~ "Wet",
    endsWith(rep, "D1") ~ "Dry",
    endsWith(rep, "D2") ~ "Dry",
    endsWith(rep, "D3") ~ "Dry",
    endsWith(rep, "D4") ~ "Dry",
    endsWith(rep, "D5") ~ "Dry"
  ))


#choosing best 4 out of 5 samples to keep using dist matrix

slope.new <- all.data %>% 
  unite(kit_treat, kit, Treat, remove = FALSE) %>% 
  group_by(kit_treat) %>% 
  mutate(Mean.Slope.All = mean(slope_of_the_regression)) %>% 
  ungroup() %>% 
  dplyr::select(-slope_beginning,-Initial_R_squared, -Final_R_squared, -R_squared_adj,-residuals,-initial_p_value,-final_p_value,-total_incubation_time_min,-number_of_points,-removed_points_high,-removed_points_beg,-removed_points_end,-breusch_p_value,-flag_r2,-flag_pos_slope,-flag_heteroscedastic)


slope.new <- slope.new %>% 
  #add_column(newColname = "cv") %>% 
  group_by(kit_treat) %>% 
  mutate(cv.before.removal = cv(slope_of_the_regression)) %>% 
ungroup()

slope.new$Slope_Removed_Mean <- NA
slope.new$cv_after_removal <- NA
slope.new$Slope_Removed <- NA
slope.new$flag <- NA

#ggplot(slope.new, aes(x = cv.before.removal)) +
  #geom_histogram(binwidth = 0.5)

#dist matrix doesn't work when there are 0's present and we set all of the positive slopes to 0 - breaks the dist matrix because of dividing by 0. hence making 0's non-zero below

# all.data = all.data %>% 
#   mutate(slope_of_the_regression = if_else(slope_of_the_regression>0,0,slope_of_the_regression)) %>% 
#   mutate(rate_mg = if_else(slope_of_the_regression>0,0,slope_of_the_regression))
# 
# slope.new.na <- slope.new %>% 
#   mutate(slope_of_the_regression = if_else(slope_of_the_regression==0,-0.0001,slope_of_the_regression)) %>% 
#   mutate(rate_mg = if_else(slope_of_the_regression==0,0.0001,slope_of_the_regression)) 

unique.samples = unique(slope.new$kit_treat)

for (i in 1:length(unique.samples)) {
  
  slope.temp = as.numeric(slope.new$slope_of_the_regression[which(slope.new$kit_treat == unique.samples[i])])
  
  slope.temp.sd <- sd(slope.temp)
  slope.temp.mean <- mean(slope.temp)
  CV = abs((slope.temp.sd/slope.temp.mean)*100)
  
  #why 1:10? looping to get 4 out of 5 best samples
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
  
  if (length(slope.temp) >= 4) {
    
    # if(CV > 10 ) {
    #   
    #   slope.new$Slope_Removed_Mean[which(slope.new$kit_treat == unique.samples[i])] = "Samples too Variable"
    #   
    # }
    
   # else {
    
   
    slope.new$Slope_Removed[which(slope.new$kit_treat == unique.samples[i])] = as.numeric(slope.temp)
      
    slope.new$cv_after_removal[which(slope.new$kit_treat == unique.samples[i])] = as.numeric(abs((sd(slope.temp)/mean(slope.temp))*100))
    
    
   slope.new$Slope_Removed_Mean[which(slope.new$kit_treat == unique.samples[i])] = as.numeric(mean(slope.temp))
   
      #  }
    
  }
  
  rm('slope.temp')
}

slope.new$flag <- ifelse(slope.new$cv.before.removal < slope.new$cv_after_removal, "Issue in dropping sample", NA)


#slope.new$rem.sample = abs(slope.new$slope_of_the_regression) - abs(slope.new$Slope_Removed)  

slope.new <- slope.new %>% 
  separate(kit_treat, c("kit", "Treat"), sep = "_", remove = FALSE)


#Histogram of all slopes from 40 mL vials with facet by wet vs. dry treatment

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/", as.character(Sys.Date()),"_all_slope_facet_histogram.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(slope.new, aes(x = slope_of_the_regression, fill = Treat))+
  geom_histogram(binwidth = 0.15)+
  facet_grid(~Treat)+
  scale_fill_brewer(palette="Set2")+
  theme_bw()+
  #ggtitle("Histogram of All Slopes")+
  theme(axis.title.x = element_text(size = 24, margin = margin(b = 5)),
        axis.title.y = element_text(size = 24),
        axis.text.x = element_text(size = 16, angle = 90, hjust = 1, vjust = 1),
        axis.text.y = element_text(size =18),
        title = element_text(size = 24),
        legend.text = element_text(size = 18))+
  ylab("Count\n") +
  xlab(expression(paste("Slope (mg " ~L^-1 ~ min^-1~")")))

dev.off()

slope.new.data = slope.new %>%
  mutate(slope_of_the_regression = if_else(slope_of_the_regression>0,0,slope_of_the_regression)) %>%
  mutate(rate_mg = if_else(slope_of_the_regression>0,0,rate_mg_per_L_per_min))

#log10 transformation of positive rate data - can't log transform negative data

slope.new.data$log_rate_mg_per_L_per_min = log10(slope.new.data$rate_mg+1)

slope.new.na$log_slope_mg_per_L_per_min = log10(slope.new.na$slope_of_the_regression)

#log10 transformation of positive rate data faceted by treatment

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/ESS-PI_EGU/", as.character(Sys.Date()),"_Log_All_Slopes_Histogram.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(slope.new.data, aes(x = log_rate_mg_per_L_per_min, fill = Treat))+
  geom_histogram(binwidth = 0.04)+
  facet_grid(~Treat)+
  scale_fill_brewer(palette="Set2")+
  theme_bw()+
  #ggtitle("Histogram of All Slopes")+
  theme(axis.title.x = element_text(size = 24, margin = margin(b = 5)),
        axis.title.y = element_text(size = 24),
        axis.text.x = element_text(size = 16, angle = 90, hjust = 1, vjust = 1),
        axis.text.y = element_text(size =18),
        title = element_text(size = 24),
        legend.text = element_text(size = 18))+
  ylab("Count\n") +
  xlab(expression(paste("Log of Rate (mg " ~L^-1 ~ min^-1~")")))

dev.off()

wet <- slope.new.data %>% 
  filter(Treat == "Wet")
#log10 histogram of wet treatments (rate + 1)

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/ESS-PI_EGU/", as.character(Sys.Date()),"_Log_All_Wet_Slopes_Histogram.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(wet, aes(x = log_rate_mg_per_L_per_min, fill = Treat))+
  geom_histogram(fill = "#FC8D62", binwidth = 0.04)+
  #facet_grid(~Treat)+
  #scale_fill_brewer(fill="66C2A5")+
  theme_bw()+
  #ggtitle("Histogram of All Slopes")+
  theme(axis.title.x = element_text(size = 24, margin = margin(b = 5)),
        axis.title.y = element_text(size = 24),
        axis.text.x = element_text(size = 16, angle = 90, hjust = 1, vjust = 1),
        axis.text.y = element_text(size =18),
        title = element_text(size = 24),
        legend.text = element_text(size = 18))+
  ylab("Count\n") +
  xlab(expression(paste("Log of Rate (mg " ~L^-1 ~ min^-1~")")))

dev.off()

dry <- slope.new.data %>% 
  filter(Treat == "Dry")

#log10 histogram of dry treatments (rate + 1)

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/ESS-PI_EGU/", as.character(Sys.Date()),"_Log_All_Dry_Slopes_Histogram.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(dry, aes(x = log_rate_mg_per_L_per_min, fill = Treat))+
  geom_histogram(fill = "#66C2A5", binwidth = 0.04)+
  #facet_grid(~Treat)+
  #scale_fill_brewer(palette="Set2")+
  theme_bw()+
  #ggtitle("Histogram of All Slopes")+
  theme(axis.title.x = element_text(size = 24, margin = margin(b = 5)),
        axis.title.y = element_text(size = 24),
        axis.text.x = element_text(size = 16, angle = 90, hjust = 1, vjust = 1),
        axis.text.y = element_text(size =18),
        title = element_text(size = 24),
        legend.text = element_text(size = 18))+
  ylab("Count\n") +
  xlab(expression(paste("Log of Rate (mg " ~L^-1 ~ min^-1~")")))

dev.off()

#Histogram of all slopes from 40 mL vials without facet by wet vs. dry treatment

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/ECA/EC 2022 Experiment/Optode multi reactor/Optode_multi_reactor_incubation/effect size/", as.character(Sys.Date()),"_all_slope_hist_no_facet.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(slope.new.na,aes(x = slope_of_the_regression))+
  geom_histogram(fill = "cornflowerblue", color = "black", binwidth = 0.08)+
  scale_fill_brewer(palette="Set2")+
  theme_bw()+
  #ggtitle("Histogram of All Slopes")+
  theme(axis.title.x = element_text(size = 24, margin = margin(b = 5)),
        axis.title.y = element_text(size = 24),
        axis.text.x = element_text(size = 16, angle = 90, hjust = 1, vjust = 1),
        axis.text.y = element_text(size =18),
        title = element_text(size = 24),
        legend.text = element_text(size = 18))+
  ylab("Count\n") +
  xlab(expression(paste("Rate (mg " ~L^-1 ~ min^-1~")")))

dev.off()

#log10 Histogram of all slopes from 40 mL vials without facet 

png(file = paste0("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/ECA/EC 2022 Experiment/Optode multi reactor/Optode_multi_reactor_incubation/effect size/", as.character(Sys.Date()),"_log_all_slope_hist_no_facet.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(slope.new.na, aes(x = log_rate_mg_per_L_per_min, fill = Treat))+
  geom_histogram(fill = "cornflowerblue", color = "black", binwidth = 0.08)+
  scale_fill_brewer(palette="Set2")+
  theme_bw()+
  #ggtitle("Histogram of All Slopes")+
  theme(axis.title.x = element_text(size = 24, margin = margin(b = 5)),
        axis.title.y = element_text(size = 24),
        axis.text.x = element_text(size = 16, angle = 90, hjust = 1, vjust = 1),
        axis.text.y = element_text(size =18),
        title = element_text(size = 24),
        legend.text = element_text(size = 18))+
  ylab("Count\n") +
  xlab(expression(paste("Log of Rate (mg " ~L^-1 ~ min^-1~")")))

dev.off()

#log10 Histogram of all slopes from 40 mL vials with facet 

png(file = paste0("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/ECA/EC 2022 Experiment/Optode multi reactor/Optode_multi_reactor_incubation/effect size/", as.character(Sys.Date()),"_log_all_slope_hist_facet.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(slope.new.na, aes(x = log_rate_mg_per_L_per_min, fill = Treat))+
  geom_histogram(binwidth = 0.08)+
  facet_grid(~Treat)+
  scale_fill_brewer(palette="Set2")+
  theme_bw()+
  #ggtitle("Histogram of All Slopes")+
  theme(axis.title.x = element_text(size = 24, margin = margin(b = 5)),
        axis.title.y = element_text(size = 24),
        axis.text.x = element_text(size = 16, angle = 90, hjust = 1, vjust = 1),
        axis.text.y = element_text(size =18),
        title = element_text(size = 24),
        legend.text = element_text(size = 18))+
  ylab("Count\n") +
  xlab(expression(paste("Log of Rate (mg " ~L^-1 ~ min^-1~")")))

dev.off()

#calculate mean slopes by Kit and Treatment (wet or dry)

slope.means <- slope.new.data %>% 
  distinct(kit_treat, .keep_all = TRUE) %>% 
  filter(Slope_Removed_Mean != "Too Variable") 

slope.means$Slope_Removed_Mean <- as.numeric(slope.means$Slope_Removed_Mean)

#Histogram of mean slopes
ggplot(slope.means, aes(x = Slope.Rem, fill = Treat))+
  geom_histogram(binwidth = 0.05)+
  facet_wrap(~Treat)+
  ggtitle("Histogram of Slopes averaged by kit and treatment")

#Bar chart of mean slopes
ggplot(slope.means, aes(x = reorder(kit,Slope.Rem), y = Slope.Rem, fill = Treat))+
  geom_bar(stat = "identity",position = position_dodge())+
  ggtitle("Mean slope by kit and treatment")+
  xlab("Kit")

#these kits were from the same site? and 27 only has 3 dry replicates
eca <- slope.means %>% 
  filter(kit != "27") 

%>% 
  filter(kit != "14")

#effect size wet - dry
eca <- eca %>% 
  group_by(kit) %>% 
  mutate(effect = (Slope_Removed_Mean[Treat == "Wet"] - Slope_Removed_Mean[Treat == "Dry"])) %>% 
  mutate(log_effect = log10(abs(effect)))



##Effect Size graph without 27 and 14
png(file = paste0("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/ESS-PI_EGU/", as.character(Sys.Date()),"_effect_size.png"), width = 10, height = 10, units = "in", res = 300)

ggplot(eca, aes(x = reorder(kit, effect), y = effect))+
  geom_bar(fill = "cornflowerblue",col = "black", stat = "summary") +
  theme_bw()+
  theme(axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text.x = element_text(size = 18, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size =18))+
  ylab("Wet - Dry Rate\n")+
  xlab("\nSite Number")

dev.off()


##Effect Size graph with 27 and 14

png(file = paste0("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/", as.character(Sys.Date()),"_wet-dry_pts_rem.png"), width = 9, height = 9, units = "in", res = 300)

ggplot(slope.means, aes(x = reorder(kit, effect), y = sub))+
  geom_bar(fill = "cornflowerblue",col = "black", stat = "summary") +
  theme_bw()+
  theme(axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text.x = element_text(size = 18, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size =18))+
  ylab("Wet - Dry Rate\n")+
  xlab("\nSite Number")

dev.off()


#histogram of effect size
png(file = paste0("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/", as.character(Sys.Date()),"_wd hist.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(slope.means, aes(x = effect))+
  geom_histogram(binwidth = 0.2, col = "black", fill = "cornflowerblue")+
  theme_bw()+
  theme(axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size =18))+
  ylab("Count\n")+
  xlab("\nWet - Dry Rate")

dev.off()

png(file = paste0("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/", as.character(Sys.Date()),"_wd_hist_log.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(slope.means, aes(x = log_effect))+
  geom_histogram(binwidth = 0.2, col = "black", fill = "cornflowerblue")+
  theme_bw()+
  theme(axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size =18))+
  ylab("Count\n")+
  xlab("\nWet - Dry Rate")

dev.off()

#write effect size results to csv

write.csv(eca, paste0(input.path,"/Effect_Size_merged_by_",pnnl.user,"_on_",Sys.Date(),".csv"), row.names = F)
