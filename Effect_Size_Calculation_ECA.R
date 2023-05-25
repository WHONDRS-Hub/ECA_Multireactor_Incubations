#### Load Library ####

library(lubridate);library(writexl);library(raster);library(tidyverse);library(devtools)

#### Load data #####
rm(list=ls());graphics.off()

# Set working directory to data file

pnnl.user = 'laan208'

input.path <- paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/rates")

setwd(input.path)

path <- ("Plots")

#change date to most recent respiration rate csv

date = '2023-05-25'

#read in all files, remove csv that are not rate data files, bind all results files together, 
import_data = function(input.path){
  
  # pull a list of file names in the target folder with the target pattern
  # then read all files and combine
  
  files <- list.files(input.path, pattern = ".csv", recursive = T, full.names = T)
  
  all <- files[grep(paste0(date), files)]
  
  data <- lapply(all, read.table, sep = ",", header = T)
  
  data <- do.call(rbind, data)
}

data = import_data(input.path)

#### Clean Data ####

#separate names, add Treatment (wet or dry), replace NAs with blanks

all.data <- data %>% 
  separate(Sample_ID, c("ECA", "kit", "rep"), sep = "_", remove = FALSE) %>% 
  mutate(Treat = case_when(grepl("W",rep)~"Wet",
                           grepl("D", rep) ~"Dry")) %>% 
  relocate(Treat, .after = rep) %>% 
  dplyr::select(c(Sample_ID,ECA,kit,rep,slope_of_the_regression,rate_mg_per_L_per_min,rate_mg_per_L_per_h,Treat))
    
          
#choosing best 4 out of 5 samples to keep using dist matrix

slope.new <- all.data %>% 
  unite(kit_treat, kit, Treat, remove = FALSE) %>% 
  group_by(kit_treat) %>% 
  mutate(Mean_Slope_All = mean(slope_of_the_regression)) %>% 
  ungroup()


slope.new <- slope.new %>% 
   group_by(kit_treat) %>% 
  mutate(cv_before_removal = cv(slope_of_the_regression)) %>% 
ungroup()

slope.new$flag <- NA

slope.final <- as.data.frame(matrix(NA, ncol = 14, nrow =1))

colnames(slope.final) = c("slope.temp","Sample_ID", "ECA", "kit_treat", "kit", "rep", "rate_mg_per_L_per_min","rate_mg_per_L_per_h", "Treat", "Mean_Slope_All","cv_before_removal", "cv_after_removal", "Mean_Slope_Removed","flag")




##if more than 4 samples and CV > 10%, then remove 1 sample

unique.samples = unique(slope.new$kit_treat)

for (i in 1:length(unique.samples)) {
  
  data_subset = subset(slope.new, slope.new$kit_treat == unique.samples[i])
    
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
  
  if (length(slope.temp) >= 4) {
    
    # if(CV > 10 ) {
    #   
    #   slope.new$Slope_Removed_Mean[which(slope.new$kit_treat == unique.samples[i])] = "Samples too Variable"
    #   
    # }
    
   # else {
   
    
    slope.combined <- as.data.frame(slope.temp)
    
    slope.removed <- merge(slope.combined, data_subset, by.x = "slope.temp", by.y = "slope_of_the_regression", all.x = TRUE)
    
    slope.removed <- slope.removed[!duplicated(slope.removed$Sample_ID), ]
    
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

#for some reason, 66-D is being added twice
slope.final <- slope.final %>% 
  distinct()

#Histogram of all slopes from 40 mL vials with facet by wet vs. dry treatment ####
png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/", as.character(Sys.Date()),"_all_slope_facet_histogram.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(slope.final, aes(x = slope_of_the_regression, fill = Treat))+
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

# Finalize respiration data with removals ####
#turn all positive rates to 0

slope.final.clean = slope.final %>%
  mutate(slope_of_the_regression = if_else(slope_of_the_regression>0,0,slope_of_the_regression)) %>%
  mutate(rate_mg = if_else(slope_of_the_regression>=0,0,rate_mg_per_L_per_min)) %>% 
  mutate(Mean_Slope_Removed = if_else(Mean_Slope_Removed>0,0,Mean_Slope_Removed)) %>% 
  dplyr::select(c(Sample_ID,rate_mg)) %>% 
  rename(rate_mg_per_L_per_min = rate_mg) %>% 
  mutate(rate_mg_per_L_per_h = rate_mg_per_L_per_min*60) %>% 
  separate(Sample_ID, into = c("EC", "kit", "INC"), remove = FALSE) %>% 
  separate(Sample_ID, into = c("ID", "rep"), sep = "-", remove = FALSE) %>% 
  mutate(Treat = case_when(grepl("W",rep)~"Wet",
                           grepl("D", rep) ~"Dry")) %>% 
  unite(kit_treat, c("kit", "Treat"),  sep = "_", remove = FALSE)


write.csv(slope.final.clean,paste0(input.path,"/removed_respiration_merged_by_",pnnl.user,"_on_",Sys.Date(),".csv"), row.names = F)

# Log Transformations and Slope Histograms####
#log10 transformation of positive rate data +1 - can't log transform negative data or 0

slope.final.clean$log_rate_mg_per_L_per_min = log10(slope.final.clean$rate_mg_per_L_per_min+1)

slope.final.clean <- na.omit(slope.final.clean)



#log10 transformation of positive rate data faceted by treatment

cbPalette <- c("#D55E00","#0072B2","#999999","#E69F00", "#56B4E9","#009E73","#F0E442","#CC79A7")

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/", as.character(Sys.Date()),"_Log_All_Slopes_Histogram.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(slope.final.clean, aes(x = log_rate_mg_per_L_per_min, fill = Treat))+
  geom_histogram(binwidth = 0.04)+
  facet_grid(~Treat)+
  scale_fill_manual(values=cbPalette, labels = c("Dry", "Wet"))+
  theme_bw()+
  #ggtitle("Histogram of All Slopes")+
  theme(axis.title.x = element_text(size = 24, margin = margin(b = 5)),
        axis.title.y = element_text(size = 24),
        axis.text.x = element_text(size = 16, angle = 90, hjust = 1, vjust = 1),
        axis.text.y = element_text(size =18),
        title = element_text(size = 24),
        legend.text = element_text(size = 18))+
  ylab("Count\n") +
  xlab(expression(paste("Log of Rate (mg O"[2]* " L"^-1* "min"^-1*")")))

dev.off()

#

wet <- slope.final.clean %>% 
  filter(Treat == "Wet")


#log10 histogram of wet treatments with dist removals(rate + 1)

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/", as.character(Sys.Date()),"_Log_All_Wet_Slopes_Histogram.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(wet, aes(x = log_rate_mg_per_L_per_min, fill = Treat))+
  geom_histogram(fill = "#0072B2", binwidth = 0.04)+
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
  xlab(expression(paste("Log of Rate (mg O"[2]* " L"^-1* " min"^-1*")")))+
  scale_y_continuous(limits = c(0,150), breaks = c(0, 25, 50, 75, 100, 125, 150))

dev.off()

dry <- slope.final.clean %>% 
  filter(Treat == "Dry")

#log10 histogram of dry treatments with removals (rate + 1)

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/", as.character(Sys.Date()),"_Log_All_Dry_Slopes_Histogram.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(dry, aes(x = log_rate_mg_per_L_per_min, fill = Treat))+
  geom_histogram(fill = "#D55E00", binwidth = 0.04)+
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
  xlab(expression(paste("Log of Rate (mg O"[2]* " L"^-1* " min"^-1*")")))+
  scale_y_continuous(limits = c(0,150), breaks = c(0, 25, 50, 75, 100, 125, 150))

dev.off()



#log10 Histogram of all slopes from 40 mL vials with facet 

png(file = paste0("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/ECA/EC 2022 Experiment/Optode multi reactor/Optode_multi_reactor_incubation/effect size/", as.character(Sys.Date()),"_log_all_slope_hist_facet.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(slope.final.clean, aes(x = log_rate_mg_per_L_per_min, fill = Treat))+
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

#Effect Size Calc ####
#calculate mean slopes by Kit and Treatment (wet or dry)

slope.means <- slope.final.clean %>% 
  separate(Sample_ID, into = c("kit", "rep"), sep = "-") %>% 
  mutate(Treat = case_when(grepl("W", rep)~"Wet",
                           grepl("D", rep) ~"Dry")) %>% group_by(kit, Treat) %>% 
  mutate(Mean_Slope_Removed = mean(rate_mg_per_L_per_min)) %>% 
  distinct(kit, Treat, .keep_all = TRUE) %>% 
  dplyr::select(c(kit,Treat,Mean_Slope_Removed))

slope.means$Mean_Slope_Removed <- as.numeric(slope.means$Mean_Slope_Removed)



#27 and 39 from same site, 56 and 57 from same site


#effect size wet - dry
eca <- slope.means %>% 
  filter(kit != "EC_27_INC") %>% 
  group_by(kit) %>% 
  
  mutate(effect = (Mean_Slope_Removed[Treat == "Wet"] - Mean_Slope_Removed[Treat == "Dry"])) %>% 
  mutate(log_effect = log10(abs(effect+1))) %>% 
  separate(kit, into = c("EC",  "Site", "INC"), sep = "_") %>% 
  dplyr::select(-c(EC, INC))


png(file = paste0("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/", as.character(Sys.Date()),"_effect_size.png"), width = 10, height = 10, units = "in", res = 300)

ggplot(eca, aes(x = reorder(Site, pos_effect), y = pos_effect))+
  geom_bar(fill = "cornflowerblue",col = "black", stat = "summary") +
  theme_bw()+
  theme(axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text.x = element_text(size = 18, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size =18))+
  ylab(expression(paste("Effect Size (Wet - Dry Rate) (mg O"[2]*" L"^-1*" min"^-1*")")))+
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

effect <- "C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Effect_Size_Data"

write.csv(eca, paste0(effect,"/Effect_Size_merged_by_",pnnl.user,"_on_",Sys.Date(),".csv"), row.names = F)
