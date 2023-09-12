library(readxl); library(tidyverse); library(stringr)

pnnl.user = 'laan208'

setwd("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/rates/Plots/")

df <- read.csv("ECA_Sediment_Incubations_Respiration_Rates_merged_by_laan208_on_2023-09-05.csv")

real <- df %>% 
  separate(Sample_ID, sep = "_", c("EC", "Site", "Treat"), remove = F) %>% 
  separate(Treat, sep = "-", c("INC", "Treat"), remove = F) %>% 
  mutate(Treat = str_sub(Treat, end = -2)) %>% 
  select(c(Site, Treat, Sample_ID, rate_mg_per_L_per_min, first_concentration, theoretical)) %>% 
  rename(th_or_real = theoretical)

unique.incubations = unique(real$Site)

all_summary<- as.data.frame(matrix(NA, ncol = 4, nrow =1))

colnames(all_summary) = c("Site", "Treat", "theoretical", "real")

summary <- as.data.frame(matrix(NA, ncol = 4, nrow = length(unique(real$Site))))

colnames(summary) = c("Site", "Treat", "theoretical", "real")



for (i in 1:length(unique.incubations)) {
  
   data_site_subset = subset(real, real$Site == unique.incubations[i])
  
  unique.treatment = unique(data_site_subset$Treat)
  
  summary <- as.data.frame(matrix(NA, ncol = 4, nrow = length(unique(real$Site))))
                          
colnames(summary) = c("Site", "Treat", "theoretical", "real")
  
  for (j in 1:length(unique.treatment)) {
    
    data_treat_subset = subset(data_site_subset, data_site_subset$Treat == unique.treatment[j])
    
    data_treat_subset$theoretical= length(which(data_treat_subset$th_or_real == "yes"))
    
    data_treat_subset$real = length(data_treat_subset$th_or_real[is.na(data_treat_subset$th_or_real)])
    
    summary$Site[j] = as.character(data_treat_subset$Site[1])
    
    summary$Treat[j] = as.character(data_treat_subset$Treat[1])

    summary$theoretical[j] = as.character(data_treat_subset$theoretical[1])
    
    summary$real[j] = as.character(data_treat_subset$real[1])
    
    }
  
 all_summary = rbind(summary, all_summary)
  
}

all_summary$real <- as.numeric(all_summary$real) 
all_summary$theoretical <- as.numeric(all_summary$theoretical)

all_summary = all_summary %>% 
  drop_na(Site) %>% 
  mutate(Difference = (real - theoretical))

write.csv(all_summary,paste0("real_or_theoretical_summary.csv"), row.names = F)
