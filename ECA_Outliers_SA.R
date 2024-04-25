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

##Functions for Correlation Matrices ####
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
  filter(Respiration_Rate_mg_DO_per_L_per_H != -9999) 

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
  mutate(Treat = ifelse(grepl("W", Sample_Name), "wet", "dry")) #%>% 
#mutate(type = ifelse(Respiration_Rate_mg_DO_per_L_per_H > 90, "theoretical", "real"))


## Keep all data ####
all_data_long = all_data_mg_L %>% 
  pivot_longer(cols = -c(Sample_Name, Sample_ID, Treat), names_to = "Variable", values_to = "Value")

ggplot(all_data_long, aes(x = Value)) +
  geom_histogram() +
  facet_wrap(~Variable, scales = "free")

all_data_corr = all_data_mg_L %>% 
  column_to_rownames("Sample_Name") %>% 
  select(-c(Sample_ID, Treat, Respiration_Rate_mg_DO_per_L_per_H, Fe_mg_per_L, Initial_Water_Mass_g, Initial_Gravimetric_Water, Lost_Gravimetric_Water, ATP_nanomol_per_L, Final_Water_Mass_g, Dry_Sediment_Mass_g, Incubation_Water_Mass_g)) %>% 
  relocate(Respiration_Rate_mg_DO_per_kg_per_H, .before = SpC)

pairs(all_data_corr, lower.panel = panel.smooth,upper.panel = panel.cor.spear, gap = 0, cex.labels = 0.5, cex = .75)

## No Removals Means
all_data_means = all_data_mg_L %>% 
  select(c(Sample_Name, Sample_ID, Treat, SpC, pH, Temp, Respiration_Rate_mg_DO_per_kg_per_H, Respiration_Rate_mg_DO_per_L_per_H, Fe_mg_per_kg, Fe_mg_per_L, ATP_nanomol_per_L, ATP_picomol_per_g, Percent_Clay, Percent_Silt, Percent_Fine_Sand, Percent_Med_Sand, Percent_Coarse_Sand, Percent_Tot_Sand, mean_ssa, Initial_Gravimetric_Water, Final_Gravimetric_Water, Lost_Gravimetric_Water)) %>%
  group_by(Sample_ID, Treat) %>% 
  summarise(across(where(is.numeric), mean)) %>% 
  mutate(across(where(is.numeric), round, 3)) %>% 
  relocate(c(Initial_Gravimetric_Water:Lost_Gravimetric_Water), .after = ATP_picomol_per_g) %>% 
  rename_with(.cols = c(SpC:Lost_Gravimetric_Water), .fn = ~ paste0("Mean_", .x)) 

all_means_long = all_data_means %>% 
  pivot_longer(cols = -c(Sample_ID, Treat), names_to = "Variable", values_to = "Value")

ggplot(all_means_long, aes(x = Value)) +
  geom_histogram() +
  facet_wrap(~Variable, scales = "free") +
  theme(strip.text = element_text(size = 6.5))

all_means_corr = all_data_means %>% 
  unite(kit_treat, c("Sample_ID", "Treat")) %>% 
  column_to_rownames("kit_treat") %>% 
  select(-c(Mean_Respiration_Rate_mg_DO_per_L_per_H, Mean_Fe_mg_per_L,  Mean_Initial_Gravimetric_Water, Mean_Lost_Gravimetric_Water, Mean_ATP_nanomol_per_L)) %>% 
  relocate(Mean_Respiration_Rate_mg_DO_per_kg_per_H, .before = Mean_SpC)

pairs(all_means_corr, lower.panel = panel.smooth,upper.panel = panel.cor.spear, gap = 0, cex.labels = 0.5, cex = .75)

## No Removals Effect

all_effect = all_means_corr %>% 
  rownames_to_column("Sample_ID") %>% 
  separate(Sample_ID, c("Sample_ID", "Treat"), sep = -4) %>% 
  filter(!grepl("EC_011|EC_012", Sample_ID)) %>% 
  group_by(Sample_ID) %>% 
mutate(across(c(Mean_Respiration_Rate_mg_DO_per_kg_per_H:Mean_Final_Gravimetric_Water), ~. [Treat == "_wet"] - .[Treat  == "_dry"])) %>% 
  rename_with(.cols = c(Mean_Respiration_Rate_mg_DO_per_kg_per_H:Mean_Final_Gravimetric_Water), .fn = ~ paste0("diff_", .x)) %>% 
  distinct(Sample_ID, .keep_all = TRUE) %>% 
  select(-c(Treat))

all_effect_long = all_effect %>% 
  pivot_longer(cols = -c(Sample_ID), names_to = "Variable", values_to = "Value")

ggplot(all_effect_long, aes(x = Value)) +
  geom_histogram() +
  facet_wrap(~Variable, scales = "free") +
  theme(strip.text = element_text(size = 6.5))

all_effect_corr = all_effect %>% 
  column_to_rownames("Sample_ID") 

pairs(all_effect_corr, lower.panel = panel.smooth,upper.panel = panel.cor.spear, gap = 0, cex.labels = 0.5, cex = .75)
 
## Keep all medians ####
all_data_medians = all_data_mg_L %>% 
select(c(Sample_Name, Sample_ID, Treat, SpC, pH, Temp, Respiration_Rate_mg_DO_per_kg_per_H, Respiration_Rate_mg_DO_per_L_per_H, Fe_mg_per_kg, Fe_mg_per_L, ATP_nanomol_per_L, ATP_picomol_per_g, Percent_Clay, Percent_Silt, Percent_Fine_Sand, Percent_Med_Sand, Percent_Coarse_Sand, Percent_Tot_Sand, mean_ssa, Initial_Gravimetric_Water, Final_Gravimetric_Water, Lost_Gravimetric_Water)) %>%
  group_by(Sample_ID, Treat) %>% 
  summarise(across(where(is.numeric), median)) %>% 
  mutate(across(where(is.numeric), round, 3))%>% 
  rename_with(.cols = c(SpC:Lost_Gravimetric_Water), .fn = ~ paste0("Median_", .x)) 

all_medians_long = all_data_medians %>% 
  pivot_longer(cols = -c(Sample_ID, Treat), names_to = "Variable", values_to = "Value")

ggplot(all_medians_long, aes(x = Value)) +
  geom_histogram() +
  facet_wrap(~Variable, scales = "free")

medians_corr = all_data_medians %>% 
  unite(kit_treat, c("Sample_ID", "Treat")) %>% 
  column_to_rownames("kit_treat") %>% 
  select(-c(Median_ATP_nanomol_per_L, Median_Fe_mg_per_L, Median_Initial_Gravimetric_Water, Median_Lost_Gravimetric_Water, Median_Respiration_Rate_mg_DO_per_L_per_H)) %>% 
  relocate(Median_Respiration_Rate_mg_DO_per_kg_per_H, .before = Median_SpC)

pairs(medians_corr, lower.panel = panel.smooth,upper.panel = panel.cor.spear, gap = 0, cex.labels = 0.5, cex = .75)


## Medians Effect

median_effect = medians_corr %>% 
  rownames_to_column("Sample_ID") %>% 
  separate(Sample_ID, c("Sample_ID", "Treat"), sep = -4) %>% 
  filter(!grepl("EC_011|EC_012", Sample_ID)) %>%
  relocate(Median_Final_Gravimetric_Water, .after = Median_ATP_picomol_per_g) %>% 
  group_by(Sample_ID) %>% 
  mutate(across(c(Median_Respiration_Rate_mg_DO_per_kg_per_H:Median_Final_Gravimetric_Water), ~. [Treat == "_wet"] - .[Treat  == "_dry"])) %>% 
  rename_with(.cols = c(Median_Respiration_Rate_mg_DO_per_kg_per_H:Median_Final_Gravimetric_Water), .fn = ~ paste0("diff_", .x)) %>% 
  distinct(Sample_ID, .keep_all = TRUE) %>% 
  select(-c(Treat))

median_effect_long = median_effect %>% 
  pivot_longer(cols = -c(Sample_ID), names_to = "Variable", values_to = "Value")

ggplot(median_effect_long, aes(x = Value)) +
  geom_histogram() +
  facet_wrap(~Variable, scales = "free") +
  theme(strip.text = element_text(size = 6.5))

median_effect_corr = median_effect %>% 
  column_to_rownames("Sample_ID") 

pairs(median_effect_corr, lower.panel = panel.smooth,upper.panel = panel.cor.spear, gap = 0, cex.labels = 0.5, cex = .75)

## Remove same outliers as respiration ####

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
  mutate(across(where(is.numeric), round, 3)) %>% 
  rename_with(.cols = c(SpC:Lost_Gravimetric_Water), .fn = ~ paste0("Best_Effect_Mean_", .x)) 

resp_out_long = resp_out_means %>% 
  pivot_longer(cols = -c(Sample_ID, Treat), names_to = "Variable", values_to = "Value")

ggplot(resp_out_long, aes(x = Value)) +
  geom_histogram() +
  facet_wrap(~Variable, scales = "free") +
  theme(strip.text = element_text(size = 6.5))

resp_out_corr = resp_out_means %>% 
  unite(kit_treat, c("Sample_ID", "Treat")) %>% 
  column_to_rownames("kit_treat") %>% 
  select(-c(Best_Effect_Mean_ATP_nanomol_per_L, Best_Effect_Mean_Fe_mg_per_L, Best_Effect_Mean_Initial_Gravimetric_Water, Best_Effect_Mean_Lost_Gravimetric_Water, Best_Effect_Mean_Respiration_Rate_mg_DO_per_L_per_H)) %>% 
  relocate(Best_Effect_Mean_Respiration_Rate_mg_DO_per_kg_per_H, .before = Best_Effect_Mean_SpC)

pairs(resp_out_corr, lower.panel = panel.smooth,upper.panel = panel.cor.spear, gap = 0, cex.labels = 0.5, cex = .75)

## Best Means Effect

resp_out_effect = resp_out_corr %>% 
  rownames_to_column("Sample_ID") %>% 
  separate(Sample_ID, c("Sample_ID", "Treat"), sep = -4) %>% 
  filter(!grepl("EC_011|EC_012", Sample_ID)) %>%
  relocate(Best_Effect_Mean_Final_Gravimetric_Water, .after = Best_Effect_Mean_ATP_picomol_per_g) %>% 
  group_by(Sample_ID) %>% 
  mutate(across(c(Best_Effect_Mean_Respiration_Rate_mg_DO_per_kg_per_H:Best_Effect_Mean_Final_Gravimetric_Water), ~. [Treat == "_wet"] - .[Treat  == "_dry"])) %>% 
  rename_with(.cols = c(Best_Effect_Mean_Respiration_Rate_mg_DO_per_kg_per_H:Best_Effect_Mean_Final_Gravimetric_Water), .fn = ~ paste0("diff_", .x)) %>% 
  distinct(Sample_ID, .keep_all = TRUE) %>% 
  select(-c(Treat))

resp_out_effect_long = resp_out_effect %>% 
  pivot_longer(cols = -c(Sample_ID), names_to = "Variable", values_to = "Value")

ggplot(resp_out_effect_long, aes(x = Value)) +
  geom_histogram() +
  facet_wrap(~Variable, scales = "free") +
  theme(strip.text = element_text(size = 6.5))

resp_out_effect_corr = resp_out_effect %>% 
  column_to_rownames("Sample_ID") 

pairs(resp_out_effect_corr, lower.panel = panel.smooth,upper.panel = panel.cor.spear, gap = 0, cex.labels = 0.5, cex = .75)

ggplot(resp_out_effect_corr, aes(y = diff_Best_Effect_Mean_Respiration_Rate_mg_DO_per_kg_per_H, x = diff_Best_Effect_Mean_Fe_mg_per_kg)) + 
  geom_point()

##Remove outliers independently ####

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

ggplot(eca_summary_all, aes(x = ATP_picomol_per_g_cv)) +
  geom_histogram()+
  geom_vline(xintercept = 30)

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

resp_means = resp_out_corr %>% 
  rownames_to_column("kit_treat") %>% 
  select(c(kit_treat, Best_Effect_Mean_Respiration_Rate_mg_DO_per_kg_per_H, Best_Effect_Mean_mean_ssa, Best_Effect_Mean_Percent_Clay, Best_Effect_Mean_Percent_Silt, Best_Effect_Mean_Percent_Fine_Sand, Best_Effect_Mean_Percent_Med_Sand, Best_Effect_Mean_Percent_Coarse_Sand, Best_Effect_Mean_Percent_Tot_Sand)) %>% 
  mutate(kit_treat = str_replace(kit_treat, "_d", "-d")) %>% 
  mutate(kit_treat = str_replace(kit_treat, "_w", "-w"))

removals_means = removal_final %>% 
  select(c(kit_treat, Mean_SpC_Removed, Mean_pH_Removed, Mean_Fe_Removed, Mean_ATP_Removed, Mean_Final_Grav_Removed)) %>% 
  drop_na() %>% 
  distinct(kit_treat, .keep_all = TRUE) %>% 
  left_join(resp_means, by = "kit_treat")

removals_long = removals_means %>% 
  pivot_longer(cols = -c(kit_treat), names_to = "Variable", values_to = "Value")

ggplot(removals_long, aes(x = Value)) +
  geom_histogram() +
  facet_wrap(~Variable, scales = "free")


removals_corr = removals_means %>% 
  column_to_rownames("kit_treat") %>% 
  relocate(Best_Effect_Mean_Respiration_Rate_mg_DO_per_kg_per_H, .before = Mean_SpC_Removed)

pairs(removals_corr, lower.panel = panel.smooth,upper.panel = panel.cor.spear, gap = 0, cex.labels = 0.5, cex = .75)

## Best Analyte Effect

removals_effect = removals_corr %>% 
  rownames_to_column("Sample_ID") %>% 
  separate(Sample_ID, c("Sample_ID", "Treat"), sep = -4) %>% 
  filter(!grepl("EC_011|EC_012", Sample_ID)) %>% 
  group_by(Sample_ID) %>% 
  mutate(across(c(Best_Effect_Mean_Respiration_Rate_mg_DO_per_kg_per_H:Mean_Final_Grav_Removed), ~. [Treat == "-wet"] - .[Treat  == "-dry"])) %>% 
  rename_with(.cols = c(Best_Effect_Mean_Respiration_Rate_mg_DO_per_kg_per_H:Mean_Final_Grav_Removed), .fn = ~ paste0("diff_", .x)) %>% 
  distinct(Sample_ID, .keep_all = TRUE) %>% 
  select(-c(Treat))

removals_effect_long = removals_effect %>% 
  pivot_longer(cols = -c(Sample_ID), names_to = "Variable", values_to = "Value")

ggplot(removals_effect_long, aes(x = Value)) +
  geom_histogram() +
  facet_wrap(~Variable, scales = "free") +
  theme(strip.text = element_text(size = 6.5))

removals_effect_corr = removals_effect %>% 
  column_to_rownames("Sample_ID") 

pairs(removals_effect_corr, lower.panel = panel.smooth,upper.panel = panel.cor.spear, gap = 0, cex.labels = 0.5, cex = .75)


ggplot(removals_effect_corr, aes(y = diff_Best_Effect_Mean_Respiration_Rate_mg_DO_per_kg_per_H, x = diff_Mean_Fe_Removed)) + 
  geom_point()

## LASSO for SA ####

# all_effect_corr
# median_effect_corr
# resp_out_effect_corr
# removals_effect_corr

## Cube Effect Size

cube_root <- function(x) sign(x) * (abs(x))^(1/3)

# To get theoretical, effect size > 3.5
cube_best_effect = median_effect_corr %>% 
  mutate(across(where(is.numeric), cube_root)) %>% 
  rename_with( ~ paste0("cube_", .), everything())# %>% 
  #select(-c(Treat)) %>% 
  #filter(Mean_WithOutliers_Respiration_Rate_mg_DO_per_L_per_H < 3.5) %>% 
  #column_to_rownames("Sample_ID") #%>% 
#mutate(Th = ifelse(Mean_OutliersRemoved_Respiration_Rate_mg_DO_per_L_per_H > 3.5, "theoretical", "real"))

cube_best_effect_long = cube_best_effect %>% 
  rownames_to_column("Sample_ID") %>% 
  pivot_longer(cols = -c(Sample_ID), names_to = "Variable", values_to = "Value")

ggplot(cube_best_effect_long, aes(x = Value)) +
  geom_histogram() +
  facet_wrap(~Variable, scales = "free") +
  theme(strip.text = element_text(size = 5))

## Scale cube effect size
z_cube_best_effect = cube_best_effect %>% 
  mutate(across(where(is.numeric), function(x) ((x - mean(x)) / sd(x))))

## EFFECT SIZE LASSO ####
## Set response variable
eff <- z_cube_best_effect$cube_diff_Best_Effect_Mean_Respiration_Rate_mg_DO_per_kg_per_H

## Set predictor variables
#Try all data
#z_cube_effect_pred <- data.matrix(z_cube_best_effect[, c('Mean_Fe_mg_per_L', 'Mean_ATP_nanomol_per_L', 'Percent_Fine_Sand', 'Percent_Coarse_Sand', 'Percent_Tot_Sand', 'Percent_Silt', 'Percent_Clay', 'mean_ssa', 'Mean_SpC', 'Mean_Temp', 'Mean_pH', 'Initial_Gravimetric_Water', 'Final_Gravimetric_Water')])

## Chosen variables: FS, SSA, Fe, Moi, ATP, SpC, Temp, pH
z_cube_effect_pred <- data.matrix(z_cube_best_effect[, c("cube_diff_Mean_Fe_Removed", "cube_diff_Mean_ATP_Removed",  "cube_Best_Effect_Mean_Percent_Fine_Sand", "cube_Best_Effect_Mean_mean_ssa", "cube_diff_Mean_SpC_Removed", "cube_diff_Mean_pH_Removed",  "cube_diff_Mean_Final_Grav_Removed")])

## Start LASSO
#alpha = 1 is for LASSO regression
cv_model <- cv.glmnet(z_cube_effect_pred, eff, alpha = 1)

#Find lambda for lowest MSE using k-fold cross-validation
best_lambda <- cv_model$lambda.min
best_lambda

plot(cv_model)

# Rerun LASSO with best lambda and get coefficients
best_model <- glmnet(z_cube_effect_pred, eff, alpha = 1, lambda = best_lambda)
coef(best_model)

## Get R Sq. of best model
eff_predict <- predict(best_model, s = best_lambda, newx = z_cube_effect_pred)

sst <- sum((eff - mean(eff))^2)
sse <- sum((eff_predict - eff)^2)

rsq = 1 - sse/sst
rsq

# Get model residuals
residuals = eff - eff_predict

par(mfrow=c(2,2)) # Set up a 2x2 grid of plots
hist(residuals, main="Histogram of Residuals")
qqnorm(residuals, main="QQ Plot of Residuals")
qqline(residuals)
plot(density(residuals), main="Density Plot of Residuals")
plot(residuals ~ eff_predict)

dev.off()

## Plots of individual vs. final water content

cube_root <- function(x) sign(x) * (abs(x))^(1/3)

cube_data = all_data_mg_L %>% 
  mutate(across(where(is.numeric), cube_root))

final_th = cube_data %>% 
  filter(Respiration_Rate_mg_DO_per_L_per_H > 4.5)

final_real = cube_data %>% 
  filter(Respiration_Rate_mg_DO_per_L_per_H < 4.5)

grav_p = ggplot(cube_data, aes(x = Final_Gravimetric_Water, y = Respiration_Rate_mg_DO_per_L_per_H)) +
  geom_point(aes(color = Treat)) +
  theme_bw() +
  stat_cor(data = final_th, label.y = 6.5, digits = 2, aes(label = paste(..rr.label.., ..p.label.., sep = "~`;`~")))+
  stat_poly_line(data = final_real, se = FALSE)+
  stat_cor(data = final_real, digits = 2, label.y = 4, aes(label = paste(..rr.label.., ..p.label.., sep = "~`;`~")))+
  stat_poly_line(data = final_th, se = FALSE)+
  xlab("Cube Root Final Grav. Water") +
  ylab("Cube Root Respiration")+ 
  theme(text = element_text(size = 9)) 

grav_p

