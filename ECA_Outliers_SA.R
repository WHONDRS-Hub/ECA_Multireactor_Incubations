#### ECA LASSO ####
library(tidyverse)
library(glmnet)
library(factoextra)
library(janitor)

rm(list=ls());graphics.off()

## Read in and clean all data ####
# remove NEON sites - 052, 053, 057
# remove samples with no respiration rate
# Total: 512 samples

# icon_resp = read.csv("C:/GitHub/ECA_Multireactor_Incubations/Data/v2_CM_SSS_Sediment_Sample_Data_Summary.csv") %>% 
#   filter(!row_number() %in% c(1, 3:13)) %>% 
#   janitor::row_to_names(row_number = 1) %>% 
#   dplyr::select(c(Sample_Name,Mean_Normalized_Respiration_Rate_mg_DO_per_H_per_L_sediment)) %>% 
#   filter(!grepl("SSS", Sample_Name)) %>% 
#   filter(Sample_Name != "") %>% 
#   mutate(Sample_Name = str_replace(Sample_Name, "CM", "EC")) %>% 
#   mutate(Sample_ID = str_remove(Sample_Name, "_Sediment")) %>% 
#   filter(Mean_Normalized_Respiration_Rate_mg_DO_per_H_per_L_sediment != -9999) %>% 
#   mutate(Mean_Normalized_Respiration_Rate_mg_DO_per_H_per_L_sediment = abs(as.numeric(Mean_Normalized_Respiration_Rate_mg_DO_per_H_per_L_sediment))) %>% 
#   select(-c(Sample_Name))

all_data = read.csv("C:/GitHub/ECA_Multireactor_Incubations/Data/Cleaned Data/All_ECA_Data_03-06-2024.csv")  %>% 
  filter(!grepl("052", Sample_Name)) %>% 
  filter(!grepl("053", Sample_Name)) %>% 
  filter(!grepl("057", Sample_Name)) %>% 
  filter(Respiration_Rate_mg_DO_per_L_per_H != -9999) %>% 
  left_join(icon_resp, by = "Sample_ID")

#remove site 023 because no moisture data
all_data_mg_L = all_data %>% 
  select(-c(INC)) %>%
  filter(!grepl("023", Sample_Name)) %>% 
  mutate(Fe_mg_per_L = if_else(grepl("SFE_Below", Fe_mg_per_L), "0.002", Fe_mg_per_L)) %>% 
  mutate(Fe_mg_per_kg = if_else(grepl("SFE_Below", Fe_mg_per_kg), "0.006", Fe_mg_per_kg)) %>% 
  mutate(Fe_mg_per_L = as.numeric(Fe_mg_per_L)) %>% 
  mutate(Respiration_Rate_mg_DO_per_L_per_H = abs(Respiration_Rate_mg_DO_per_L_per_H)) %>% 
  mutate(Fe_mg_per_kg = as.numeric(Fe_mg_per_kg)) %>% 
  mutate(Respiration_Rate_mg_DO_per_kg_per_H = abs(Respiration_Rate_mg_DO_per_kg_per_H)) %>% 
  rename(ICON_resp = Mean_Normalized_Respiration_Rate_mg_DO_per_H_per_L_sediment) %>% 
  mutate(Treat = ifelse(grepl("W", Sample_Name), "wet", "dry")) #%>% 
#mutate(type = ifelse(Respiration_Rate_mg_DO_per_L_per_H > 90, "theoretical", "real"))


## Keep all means
all_data_means = all_data_mg_L %>% 
  select(c(Sample_Name, Sample_ID, Treat, SpC, pH, Temp, Respiration_Rate_mg_DO_per_kg_per_H, Respiration_Rate_mg_DO_per_L_per_H, Fe_mg_per_kg, Fe_mg_per_L, ATP_nanomol_per_L, ATP_picomol_per_g, Percent_Clay, Percent_Silt, Percent_Fine_Sand, Percent_Med_Sand, Percent_Coarse_Sand, Percent_Tot_Sand, mean_ssa, Initial_Gravimetric_Water, Final_Gravimetric_Water, Lost_Gravimetric_Water)) %>%
  group_by(Sample_ID, Treat) %>% 
  summarise(across(where(is.numeric), mean)) %>% 
  mutate(across(where(is.numeric), round, 3))
 
## Keep all medians 
all_data_medians = all_data_mg_L %>% 
select(c(Sample_Name, Sample_ID, Treat, SpC, pH, Temp, Respiration_Rate_mg_DO_per_kg_per_H, Respiration_Rate_mg_DO_per_L_per_H, Fe_mg_per_kg, Fe_mg_per_L, ATP_nanomol_per_L, ATP_picomol_per_g, Percent_Clay, Percent_Silt, Percent_Fine_Sand, Percent_Med_Sand, Percent_Coarse_Sand, Percent_Tot_Sand, mean_ssa, Initial_Gravimetric_Water, Final_Gravimetric_Water, Lost_Gravimetric_Water)) %>%
  group_by(Sample_ID, Treat) %>% 
  summarise(across(where(is.numeric), median)) %>% 
  mutate(across(where(is.numeric), round, 3))


## Remove same outliers as respiration 

resp_out = read.csv("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/INC/03_ProcessedData/ECA_Sediment_Incubations_Respiration_Rates_ReadyForBoye_2024-04-05.csv") %>% 
  left_join(all_data) %>% 
  filter(!grepl("RATE_OUTLIER_000", Methods_Deviation))

resp_out_means = resp_out %>% 
  mutate(Fe_mg_per_L = if_else(grepl("SFE_Below", Fe_mg_per_L), "0.002", Fe_mg_per_L)) %>% 
  mutate(Fe_mg_per_kg = if_else(grepl("SFE_Below", Fe_mg_per_kg), "0.006", Fe_mg_per_kg)) %>% 
  mutate(Fe_mg_per_L = as.numeric(Fe_mg_per_L)) %>% 
  mutate(Respiration_Rate_mg_DO_per_L_per_H = abs(Respiration_Rate_mg_DO_per_L_per_H)) %>% 
  mutate(Fe_mg_per_kg = as.numeric(Fe_mg_per_kg)) %>% 
  mutate(Respiration_Rate_mg_DO_per_kg_per_H = abs(Respiration_Rate_mg_DO_per_kg_per_H)) %>% 
  filter(Sample_ID != "EC_023") %>% 
  mutate(Treat = ifelse(grepl("W", Sample_Name), "wet", "dry")) %>% 
  select(c(Sample_Name, Sample_ID, Treat, SpC, pH, Temp, Respiration_Rate_mg_DO_per_kg_per_H, Respiration_Rate_mg_DO_per_L_per_H, Fe_mg_per_kg, Fe_mg_per_L, ATP_nanomol_per_L, ATP_picomol_per_g, Percent_Clay, Percent_Silt, Percent_Fine_Sand, Percent_Med_Sand, Percent_Coarse_Sand, Percent_Tot_Sand, mean_ssa, Initial_Gravimetric_Water, Final_Gravimetric_Water, Lost_Gravimetric_Water)) %>%
  group_by(Sample_ID, Treat) %>% 
  summarise(across(where(is.numeric), mean)) %>% 
  mutate(across(where(is.numeric), round, 3))

##Remove outliers independently

cv <- function(x) {
  cv_value <- (sd(x) / mean(x)) * 100
  return(cv_value)
}



eca_summary_all = all_data_mg_L %>% 
  select(c(Sample_Name, Sample_ID, Treat, SpC, pH, Temp, Respiration_Rate_mg_DO_per_kg_per_H, Respiration_Rate_mg_DO_per_L_per_H, Fe_mg_per_kg, Fe_mg_per_L, ATP_nanomol_per_L, ATP_picomol_per_g, Percent_Clay, Percent_Silt, Percent_Fine_Sand, Percent_Med_Sand, Percent_Coarse_Sand, Percent_Tot_Sand, mean_ssa, Initial_Gravimetric_Water, Final_Gravimetric_Water, Lost_Gravimetric_Water)) %>%
  group_by(Sample_ID, Treat) %>% 
  #summarise(across(c(SpC:ATP_picomol_per_g), list(mean = mean, sd = sd, cv = cv))) %>%
  summarise(across(where(is.numeric), list(mean = mean, sd = sd, cv = cv))) %>% 
  mutate(round(across(where(is.numeric)), 3)) %>% 
  select(-c(Percent_Clay_sd, Percent_Clay_cv, Percent_Silt_sd, Percent_Silt_cv, Percent_Fine_Sand_sd, Percent_Fine_Sand_cv, Percent_Med_Sand_sd, Percent_Med_Sand_cv, Percent_Coarse_Sand_sd, Percent_Coarse_Sand_cv, Percent_Tot_Sand_sd, Percent_Tot_Sand_cv, mean_ssa_sd, mean_ssa_cv))

ggplot(eca_summary_all, aes(x = Final_Gravimetric_Water_cv)) +
  geom_histogram()

cv_outliers = all_data_mg_L %>% 
  select(c(Sample_Name, Sample_ID, Treat, SpC, Temp, pH, Fe_mg_per_kg, ATP_picomol_per_g, Final_Gravimetric_Water)) %>% 
  unite(kit_treat, c(Sample_ID, Treat), sep = "-")

cv_outliers$flag <- NA

removal_final <- as.data.frame(matrix(NA, ncol = 17, nrow =1))

colnames(removal_final) = c("spc.temp", "ph.temp", "fe.temp", "atp.temp", "grav.temp", "Sample_Name", "kit_treat", "Mean_SpC_Removed", "spc_cv_after_removal", "Mean_pH_Removed", "ph_cv_after_removal", "Mean_Fe_Removed", "fe_cv_after_removal", "Mean_ATP_Removed", "atp_cv_after_removal", "Mean_Final_Grav_Removed", "grav_cv_after_removal")

unique.samples = unique(cv_outliers$kit_treat)

#respiration = 50%
#spc = 20% (MilliQ vs. real)
#ph = 7.5%
#temp = N/A
#Fe = 50% - lower than suggested from histogram
#ATP = 30% - same as ICON, lots being removed
#Gravimetric Water = 15%

spc.threshold = 20
ph.threshold = 7.5
fe.threshold = 50
atp.threshold = 30
grav.threshold = 15

#try keeping n = 3,4 samples
rem.threshold = 3

for (i in 1:length(unique.samples)) {
  
  ## Subset replicates
  data_subset = subset(cv_outliers, cv_outliers$kit_treat == unique.samples[i])
  
  ## Pull out Rate values
  #slope.temp = as.numeric(data_subset$Respiration_Rate_mg_DO_per_L_per_H)
  spc.temp = as.numeric(data_subset$SpC)
  ph.temp = as.numeric(data_subset$pH)
  fe.temp = as.numeric(data_subset$Fe_mg_per_kg)
  atp.temp = as.numeric(data_subset$ATP_picomol_per_g)
  grav.temp = as.numeric(data_subset$Final_Gravimetric_Water)
  
  ## Calculate standard deviation, average, and coefficient of variation of rates
  spc.temp.sd <- sd(spc.temp)
  spc.temp.mean <- mean(spc.temp)
  spc.cv = abs((spc.temp.sd/spc.temp.mean)*100)
  
  ph.temp.sd <- sd(ph.temp)
  ph.temp.mean <- mean(ph.temp)
  ph.cv = abs((ph.temp.sd/ph.temp.mean)*100)
  
  fe.temp.sd <- sd(fe.temp)
  fe.temp.mean <- mean(fe.temp)
  fe.cv = abs((fe.temp.sd/fe.temp.mean)*100)
  
  atp.temp.sd <- sd(atp.temp)
  atp.temp.mean <- mean(atp.temp)
  atp.cv = abs((atp.temp.sd/atp.temp.mean)*100)
  
  grav.temp.sd <- sd(grav.temp)
  grav.temp.mean <- mean(grav.temp)
  grav.cv = abs((grav.temp.sd/grav.temp.mean)*100)
  
  #looping to get 3 out of 5 best samples for SpC
  for (sample.reduction in 1:5)  {
    
    if (spc.temp.mean == 0) {
      
      spc.cv = 0
      
    }
    
    else if (length(spc.temp) > rem.threshold & spc.cv >= spc.threshold) {
      
      dist.temp = as.matrix(abs(dist(spc.temp)))
      dist.comp = numeric()
      
      for(spc.now in 1:ncol(dist.temp)) {
        
        dist.comp = rbind(dist.comp,c(spc.now,sum(dist.temp[,spc.now])))
        
      }
      
      dist.comp[,2] = as.numeric(dist.comp[,2])
      spc.temp = spc.temp[-which.max(dist.comp[,2])]
      
      spc.temp.sd <- sd(spc.temp)
      spc.temp.mean <- mean(spc.temp)
      spc.temp.cv <- abs((spc.temp.sd/spc.temp.mean)*100)
      spc.cv = spc.temp.cv

    }
  }
  
  #looping to get 3 out of 5 best samples for pH
  for (sample.reduction in 1:5)  {
    
    if (ph.temp.mean == 0) {
      
      ph.cv = 0
      
    }
    
    else if (length(ph.temp) > rem.threshold & ph.cv >= ph.threshold) {
      
      dist.temp = as.matrix(abs(dist(ph.temp)))
      dist.comp = numeric()
      
      for(ph.now in 1:ncol(dist.temp)) {
        
        dist.comp = rbind(dist.comp,c(ph.now,sum(dist.temp[,ph.now])))
        
      }
      
      dist.comp[,2] = as.numeric(dist.comp[,2])
      ph.temp = ph.temp[-which.max(dist.comp[,2])]
      
      ph.temp.sd <- sd(ph.temp)
      ph.temp.mean <- mean(ph.temp)
      ph.temp.cv <- abs((ph.temp.sd/ph.temp.mean)*100)
      ph.cv = ph.temp.cv
      
    }
  }
  
  #looping to get 3 out of 5 best samples for Fe
  for (sample.reduction in 1:5)  {
    
    if (fe.temp.mean == 0) {
      
      fe.cv = 0
      
    }
    
    else if (length(fe.temp) > rem.threshold & fe.cv >= fe.threshold) {
      
      dist.temp = as.matrix(abs(dist(fe.temp)))
      dist.comp = numeric()
      
      for(fe.now in 1:ncol(dist.temp)) {
        
        dist.comp = rbind(dist.comp,c(fe.now,sum(dist.temp[,fe.now])))
        
      }
      
      dist.comp[,2] = as.numeric(dist.comp[,2])
      fe.temp = fe.temp[-which.max(dist.comp[,2])]
      
      fe.temp.sd <- sd(fe.temp)
      fe.temp.mean <- mean(fe.temp)
      fe.temp.cv <- abs((fe.temp.sd/fe.temp.mean)*100)
      fe.cv = fe.temp.cv
      
    }
  }
  
  #looping to get 3 out of 5 best samples for ATP
  for (sample.reduction in 1:5)  {
    
    if (atp.temp.mean == 0) {
      
      atp.cv = 0
      
    }
    
    else if (length(atp.temp) > rem.threshold & atp.cv >= atp.threshold) {
      
      dist.temp = as.matrix(abs(dist(atp.temp)))
      dist.comp = numeric()
      
      for(atp.now in 1:ncol(dist.temp)) {
        
        dist.comp = rbind(dist.comp,c(atp.now,sum(dist.temp[,atp.now])))
        
      }
      
      dist.comp[,2] = as.numeric(dist.comp[,2])
      atp.temp = atp.temp[-which.max(dist.comp[,2])]
      
      atp.temp.sd <- sd(atp.temp)
      atp.temp.mean <- mean(atp.temp)
      atp.temp.cv <- abs((atp.temp.sd/atp.temp.mean)*100)
      atp.cv = atp.temp.cv
      
    }
  }
  
  #looping to get 3 out of 5 best samples for Final Gravimetric Water
  for (sample.reduction in 1:5)  {
    
    if (grav.temp.mean == 0) {
      
      grav.cv = 0
      
    }
    
    else if (length(grav.temp) > rem.threshold & grav.cv >= grav.threshold) {
      
      dist.temp = as.matrix(abs(dist(grav.temp)))
      dist.comp = numeric()
      
      for(grav.now in 1:ncol(dist.temp)) {
        
        dist.comp = rbind(dist.comp,c(grav.now,sum(dist.temp[,grav.now])))
        
      }
      
      dist.comp[,2] = as.numeric(dist.comp[,2])
      grav.temp = grav.temp[-which.max(dist.comp[,2])]
      
      grav.temp.sd <- sd(grav.temp)
      grav.temp.mean <- mean(grav.temp)
      grav.temp.cv <- abs((grav.temp.sd/grav.temp.mean)*100)
      grav.cv = grav.temp.cv
      
    }
  }
  
  
  if (length(spc.temp) >= 3) {
    
    spc.combined <- as.data.frame(spc.temp)
    
    spc.removed <- merge(spc.combined, data_subset, by.x = "spc.temp", by.y = "SpC", all.x = TRUE) %>% 
      select(c(spc.temp, Sample_Name, kit_treat))
    
    spc.removed <- spc.removed[!duplicated(spc.removed$Sample_Name), ]
    
    spc.removed$spc_cv_after_removal = as.numeric(abs((sd(spc.temp)/mean(spc.temp))*100))
    
    spc.removed$Mean_SpC_Removed = as.numeric(mean(spc.temp))
    
  }
  
  if (length(ph.temp) >= 3) {
    
    ph.combined <- as.data.frame(ph.temp)
    
    ph.removed <- merge(ph.combined, data_subset, by.x = "ph.temp", by.y = "pH", all.x = TRUE) %>% 
      select(c(ph.temp, Sample_Name, kit_treat))
    
    ph.removed <- ph.removed[!duplicated(ph.removed$Sample_Name), ]
    
    ph.removed$ph_cv_after_removal = as.numeric(abs((sd(ph.temp)/mean(ph.temp))*100))
    
    ph.removed$Mean_pH_Removed = as.numeric(mean(ph.temp))
    
  }
  
  if (length(fe.temp) >= 3) {
    
    fe.combined <- as.data.frame(fe.temp)
    
    fe.removed <- merge(fe.combined, data_subset, by.x = "fe.temp", by.y = "Fe_mg_per_kg", all.x = TRUE) %>% 
      select(c(fe.temp, Sample_Name, kit_treat))
    
    fe.removed <- fe.removed[!duplicated(fe.removed$Sample_Name), ]
    
    fe.removed$fe_cv_after_removal = as.numeric(abs((sd(fe.temp)/mean(fe.temp))*100))
    
    fe.removed$Mean_Fe_Removed = as.numeric(mean(fe.temp))
    
  }
  
  if (length(atp.temp) >= 3) {
    
    atp.combined <- as.data.frame(atp.temp)
    
    atp.removed <- merge(atp.combined, data_subset, by.x = "atp.temp", by.y = "ATP_picomol_per_g", all.x = TRUE)%>% 
      select(c(atp.temp, Sample_Name, kit_treat))
    
    atp.removed <- atp.removed[!duplicated(atp.removed$Sample_Name), ]
    
    atp.removed$atp_cv_after_removal = as.numeric(abs((sd(atp.temp)/mean(atp.temp))*100))
    
    atp.removed$Mean_ATP_Removed = as.numeric(mean(atp.temp))
    
  }
  
  if (length(grav.temp) >= 3) {
    
    grav.combined <- as.data.frame(grav.temp)
    
    grav.removed <- merge(grav.combined, data_subset, by.x = "grav.temp", by.y = "Final_Gravimetric_Water", all.x = TRUE) %>% 
      select(c(grav.temp, Sample_Name, kit_treat))
    
    grav.removed <- grav.removed[!duplicated(grav.removed$Sample_Name), ]
    
    grav.removed$grav_cv_after_removal = as.numeric(abs((sd(grav.temp)/mean(grav.temp))*100))
    
    grav.removed$Mean_Final_Grav_Removed = as.numeric(mean(grav.temp))
    
  }
  
  all_df = list(spc.removed, ph.removed, fe.removed, atp.removed, grav.removed)
  
  merged = Reduce(function(x,y) merge(x, y, by = c("Sample_Name", "kit_treat"), all = TRUE), all_df)
  
  removal_final = rbind(removal_final, merged)
  
}

## This data frame has removed samples
removal_flag <- removal_final %>% 
  drop_na(Sample_Name) %>% 
  mutate(flag = ifelse(is.na(spc.temp), "SPC_OUT", "N/A")) %>% 
  mutate(flag = ifelse(is.na(ph.temp) & flag == "N/A", "pH_OUT", ifelse(is.na(ph.temp) & flag != "N/A", paste0(flag, "; pH_OUT"), flag))) %>% 
  mutate(flag = ifelse(is.na(fe.temp) & flag == "N/A", "Fe_OUT", ifelse(is.na(fe.temp) & flag != "N/A", paste0(flag, "; Fe_OUT"), flag))) %>% 
  mutate(flag = ifelse(is.na(atp.temp) & flag == "N/A", "ATP_OUT", ifelse(is.na(atp.temp) & flag != "N/A", paste0(flag, "; ATP_OUT"), flag))) %>% 
  mutate(flag = ifelse(is.na(grav.temp) & flag == "N/A", "Grav_OUT", ifelse(is.na(grav.temp) & flag != "N/A", paste0(flag, "; Grav_OUT"), flag))) 
