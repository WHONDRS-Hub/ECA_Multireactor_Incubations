#read in libraries

library(lubridate);library(writexl);library(raster);library(tidyverse);library(devtools)
library(readxl)
library(corrplot)

##### Load data ######
rm(list=ls());graphics.off()

# Set working directory to data file
#Example:
pnnl.user = 'laan208'


#Read in all data
setwd(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/"))

###effect size - change date to most recent
effect_size <- read_csv(paste0("Optode multi reactor/Optode_multi_reactor_incubation/rates/Effect_Size_merged_by_laan208_on_2023-04-11.csv"))


#Respiration rates 
respiration <- read.csv(paste0("Optode multi reactor/Optode_multi_reactor_incubation/rates/removed_respiration_merged_by_laan208_on_2023-04-12.csv"))
#
 


#ECA Iron
iron <- read_csv(paste0("FE/03_ProcessedData/20230110_Data_Processed_SFE_ECA_EC_1-270/20230110_Data_Processed_SFE_ECA_EC_1-270.csv"))


#ICON Grain Size
grain <- read_csv(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ICON_ModEx_SSS/09_Grain_Size/03_ProcessedData/20221230_Grain_Size_SBR_RC4_CM_1-42/20221230_Data_Processed_Grain_Size_SBR_RC4_CM_1-42.csv"))


#All incubation pH, SpC, temp
chemistry <- paste0("Optode multi reactor/Optode_multi_reactor_incubation/")

import_data = function(chemistry){ 
  
  #pull a list of files in target folder with correct pattern
  #read all files and combine
  
  map.file <-  list.files(chemistry, recursive = T, pattern = "\\.xlsx$",full.names = T)
  
  map.file <- map.file[!grepl("red|Red|EC_01_|EC_02_|EC_03_|EC_06_07|EC_10_15|EC_04_08|fast|combined|SPC|IC", map.file)]
  
  mapping <- lapply(map.file, read_xlsx)
  
  for (i in 1:length(mapping)){mapping[[i]] <- cbind(mapping[[i]], map.file[i])}
  
  all.map <- 
    do.call(rbind,mapping)
}

map = import_data(chemistry)

#Gravimetric Moisture

grav_inc <- read.csv(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/INC/03_ProcessedData/ECA_Drying_Masses_merged_by_laan208_on_2023-04-13.csv"))

wet_wt <- read.csv(paste0("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/MOI/03_ProcessedData/EC_Moisture_Content_2022.csv"))

###EFFECT SIZE/RESPIRATION###

#Use this for individual samples correlation matrix

resp <- respiration %>% 
  mutate(slope_of_the_regression = if_else(slope_of_the_regression>0,0,slope_of_the_regression)) %>%
  mutate(rate_mg = if_else(slope_of_the_regression>=0,0,rate_mg_per_L_per_min)) %>% 
  mutate(Slope_Removed_Mean = if_else(Slope_Removed_Mean>0,0,Slope_Removed_Mean))


#Use this for wet/dry correlation matrix
effect_all <- effect_size %>%
  mutate(Average_Rate = abs(Slope_Removed_Mean)) %>% 
  dplyr::select(c(kit_treat,kit,Treat,Average_Rate,pos_effect))
 
#Use this for effect size correlation matrix
effect_diff <- effect_all %>% 
  distinct(kit, .keep_all = TRUE) %>% 
  dplyr::select(-c(Average_Rate, kit_treat,Treat))

###IRON DATA

#look at individual samples

#calculate mean Fe for kit/treatment from analytical replicates

mean_fe <- iron %>% 
  separate(Sample_Name, into = c("Sample_ID", "rep"), sep = -1, convert = TRUE) %>% 
  group_by(Sample_ID) %>% 
  mutate(Mean_Fe_mg_per_L = mean(Fe_mg_per_L)) %>% 
  mutate(Fe_mg_kg= mean(Fe_mg_per_kg_sediment)) %>% 
 # summarise_at(vars(Fe_mg_per_kg_sediment), list(mean = mean)) %>% 
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
  )) %>% 
 # rename("Fe_mg_kg" = "mean") %>% 
  mutate(Log_Fe_mg_kg = log10(Fe_mg_kg)) %>% 
  unite(kit_treat, kit, Treat, remove = FALSE) %>% 
  mutate(rep = str_replace(rep, "SFE", "INC")) %>% 
  mutate(Sample_ID = str_replace(Sample_ID, "SFE", "INC")) %>% 
  distinct(Sample_ID, .keep_all = TRUE)

mean_fe_check <- mean_fe %>% 
  group_by(kit_treat) %>% 
  mutate(CV = (sd(Fe_mg_kg)/mean(Fe_mg_kg))*100)

mean_fe_check$Fe_mg_kg <- format(mean_fe_check$Fe_mg_kg, scientific = FALSE)

##All Fe data
# ggplot(mean_fe, aes(x = Fe_mg_kg))+
#   geom_histogram(binwidth = 0.1, fill = "cornflowerblue", col = "black")+
#   theme_bw()+
#   theme(axis.title.x = element_text(size = 18),
#         axis.title.y = element_text(size = 18),
#         axis.text.x = element_text(size = 16),
#         axis.text.y = element_text(size = 16))+
#   xlab("\n Fe (II) (mg/L)")+
#   ylab("Count\n")

ggplot(mean_fe, aes(y = Fe_mg_kg, x = kit_treat))+
    geom_boxplot()+
    theme_bw()+
    theme(axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16))+
    xlab("\n Fe (II) (mg/L)")+
    ylab("Count\n")
  


#All Fe data, faceted by wet vs dry
# ggplot(mean_fe, aes(x = Mean_Fe_mg_kg, fill = Treat))+
#   geom_histogram(binwidth = 0.1)+
#   facet_wrap(~Treat)+
#   theme_bw()+
#   theme(axis.title.x = element_text(size = 18), 
#         axis.title.y = element_text(size = 18), 
#         axis.text.x = element_text(size = 16),
#         axis.text.y = element_text(size = 16))+
#   xlab("\n Fe (II) (mg/L)")+
#   ylab("Count\n")



#png(file = paste0("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/ECA/EC 2022 Experiment/Optode multi reactor/Optode_multi_reactor_incubation/effect size/", as.character(Sys.Date()),"_log_fe_hist.png"), width = 8, height = 8, units = "in", res = 300)

# ggplot(mean_fe, aes(x = log_iron_mean))+
#   geom_histogram(binwidth = 0.1, fill = "cornflowerblue", col = "black")+
#   theme_bw()+
#   theme(axis.title.x = element_text(size = 18), 
#         axis.title.y = element_text(size = 18), 
#         axis.text.x = element_text(size = 16),
#         axis.text.y = element_text(size = 16))+
#   xlab("\n Log of Fe (II) (mg/L)")+
#   ylab("Count\n")
# 
# #dev.off()

#Log Fe, facetted by wet vs. dry
# ggplot(mean_fe, aes(x = log_iron_mean, fill = Treat))+
#   geom_histogram(binwidth = 0.1)+
#   facet_wrap(~Treat)+
#   theme_bw()+
#   theme(axis.title.x = element_text(size = 18), 
#         axis.title.y = element_text(size = 18), 
#         axis.text.x = element_text(size = 16),
#         axis.text.y = element_text(size = 16))+
#   xlab("\n Log of Fe (II) (mg/L)")+
#   ylab("Count\n")

#Individual Fe Samples
fe_all <- mean_fe %>% 
  mutate(rep = str_replace(rep, "SFE", "INC")) %>% 
  dplyr::select(-c(rep, ECA))

#Fe Means for use in wet/dry correlation matrix
mean_fe_treat <- mean_fe %>% 
  group_by(kit_treat) %>%
  mutate(Mean_Fe_mg_kg = mean(Fe_mg_kg)) %>% 
  mutate(Log_Mean_Fe_mg_kg = log10(Mean_Fe_mg_kg)) %>% 
  dplyr::select(c(kit_treat, kit, Treat, Mean_Fe_mg_kg)) %>% 
  distinct(kit_treat, .keep_all = TRUE)

#Differences for Effect Size Correlation Matrix
mean_fe_diff <- mean_fe_treat %>% 
  group_by(kit) %>% 
  mutate(Fe_Difference = (Mean_Fe_mg_kg[Treat == "Wet"] - Mean_Fe_mg_kg[Treat == "Dry"])) %>% 
  dplyr::select(-c(Mean_Fe_mg_kg, kit_treat ,Treat)) %>% 
  distinct(kit, .keep_all = TRUE)
 


###Grain Size Data

#this is the only data type that is not vial specific

grn <- grain %>% 
  separate(Sample_ID, c("CM", "kit", "an"), sep = "_", remove = FALSE) %>%
  filter(kit != "001") %>% 
  filter(kit != "002") %>% 
  filter(kit != "003") %>% 
  filter(kit != "004") %>% 
  filter(kit != "006") %>% 
  filter(kit != "007") %>% 
  filter(kit != "008") %>% 
  filter(kit != "010") %>% 
  filter(kit != "015") %>% 
  filter(kit != "020") %>% 
  filter(kit != "028") %>% 
  filter(kit != "043") %>% 
  filter(kit != "050") %>% 
  filter(kit != "062") %>% 
  mutate(Percent_Mud = Percent_Clay + Percent_Silt) %>% 
  dplyr::select(-c(CM, an, Sample_ID)) 

grn$kit <- sub('.', '', grn$kit)

##Dg = exp(0.01 * sum(%mass * ln(mean particle size)))
#Coarse sand = , med sand = , fine sand = , silt = , clay =


# grn_long <- grn %>% 
#   dplyr::select(-c("CM", "kit", "an")) %>% 
#   pivot_longer(!Sample_ID, names_to = "size", values_to = "percent") %>% 
#   filter(size != "Percent_Tot_Sand")

# png(file = paste0("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/", as.character(Sys.Date()),"_GRN_stacked.png"), width = 8, height = 8, units = "in", res = 300)
# 
# ggplot(grn_long, aes(fill = size, y = percent, x = Sample_ID)) + 
#   geom_bar(position="fill", stat="identity")+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust =1))
# 
# dev.off()


###pH, SpC, temp data

#Use this for individual correlation matrix
chem_all = map %>% 
  dplyr::select(-`Disk Calibration date`, -Location, -`Time on`, -`Time off`, -Notes, -`Disk_ID`, -`map.file[i]`) %>% 
  filter(!grepl("Blank", Sample_ID))

#Use this for wet/dry correlation matrices
mean_chem <- chem_all %>% 
  separate(Sample_ID, c("ECA", "kit", "Analysis"), sep = "_", remove = FALSE)%>% 
  mutate(Treat = case_when(
    endsWith(Sample_ID, "W1") ~ "Wet",
    endsWith(Sample_ID, "W2") ~ "Wet",
    endsWith(Sample_ID, "W3") ~ "Wet",
    endsWith(Sample_ID, "W4") ~ "Wet", 
    endsWith(Sample_ID, "W5") ~ "Wet",
    endsWith(Sample_ID, "D1") ~ "Dry",
    endsWith(Sample_ID, "D2") ~ "Dry",
    endsWith(Sample_ID, "D3") ~ "Dry",
    endsWith(Sample_ID, "D4") ~ "Dry",
    endsWith(Sample_ID, "D5") ~ "Dry"
  )) %>% 
  unite(kit_treat, kit, Treat, remove = FALSE ) %>% 
  group_by(kit_treat) %>% 
  mutate(Mean_Sp_Conductivity = mean(SpC)) %>% 
  mutate(Mean_Temperature = mean(Temp)) %>% 
  mutate(Mean_pH = mean(pH)) %>% 
  dplyr::select(c(kit_treat, kit, Treat,Mean_Sp_Conductivity, Mean_Temperature, Mean_pH)) %>% 
  distinct(kit_treat, .keep_all = TRUE)

#Use this for effect size matrix
mean_chem_diff <- mean_chem %>% 
  group_by(kit) %>% 
  mutate(Sp_Conductivity_Difference = (Mean_Sp_Conductivity[Treat == "Wet"] - Mean_Sp_Conductivity[Treat == "Dry"])) %>% 
  mutate(Temp_Difference = (Mean_Temperature[Treat == "Wet"] - Mean_Temperature[Treat == "Dry"])) %>% 
  mutate(pH_Difference = (Mean_pH[Treat == "Wet"] - Mean_pH[Treat == "Dry"])) %>% 
  dplyr::select(-c(Mean_Sp_Conductivity,Mean_Temperature, Mean_pH,kit_treat,Treat)) %>% 
  distinct(kit, .keep_all = TRUE)  

###Gravimetric Moisture 

mean_wet_wt <- wet_wt %>% 
  separate(col = sample_name, into = c("Project", "kit", "analysis"), sep = "_") %>%  
  separate(col = analysis, into = c("Analysis", "Replicate"), sep = "-") %>% 
  group_by(kit) %>% 
  summarise(mean_wet_grav = (mean(percent_water_content_wet)/100))

grav <- grav_inc %>% 
  group_by(Sample_Name) %>% 
  mutate(Sample_weight_initial_g = first(Sample_weight_g)) %>% 
  mutate(Sample_weight_final_g = last(Sample_weight_g)) %>% 
  ungroup() %>% 
  separate(col = Sample_Name, into = c("Project", "kit", "analysis"), sep = "_") %>%  
  separate(col = analysis, into = c("Analysis", "Replicate"), sep = "-")


final_grav <- merge(grav, mean_wet_wt, by = "kit")

all_grav <- final_grav %>% 
  mutate(mass_sed = (Sample_weight_initial_g - (Sample_weight_initial_g * mean_wet_grav))) %>% 
  mutate(mass_water_initial = (Sample_weight_initial_g * mean_wet_grav) + Water_added_initial_g) %>%  
  mutate(grav_dry_initial = mass_water_initial/mass_sed) %>% 
  mutate(mass_water_final_g = (Sample_weight_final_g - mass_sed)+Water_added_initial_g) %>% 
  mutate(grav_dry_final = mass_water_final_g/mass_sed) %>% 
  mutate(lost_grav_perc = grav_dry_initial - grav_dry_final)

#Use this for individual correlation matrix
all_grav_final <- all_grav %>% 
  mutate(Treat = case_when(
    endsWith(Replicate, "W1") ~ "Wet",
    endsWith(Replicate, "W2") ~ "Wet",
    endsWith(Replicate, "W3") ~ "Wet",
    endsWith(Replicate, "W4") ~ "Wet", 
    endsWith(Replicate, "W5") ~ "Wet",
    endsWith(Replicate, "D1") ~ "Dry",
    endsWith(Replicate, "D2") ~ "Dry",
    endsWith(Replicate, "D3") ~ "Dry",
    endsWith(Replicate, "D4") ~ "Dry",
    endsWith(Replicate, "D5") ~ "Dry"
  )) %>% 
  unite(kit_treat, kit, Treat, remove = FALSE) %>% 
  unite(Sample, Project, kit, Analysis, sep = "_", remove = FALSE) %>% 
 unite(Sample_ID, Sample, Replicate, sep = "-") %>% 
  dplyr::select(-c(Date, Tare_weight_g, Sample_weight_g,Project, Analysis,Water_added_initial_g,Sample_weight_final_g)) 

all_grav_ind <- all_grav_final %>% 
  distinct(Sample_ID, .keep_all = TRUE) %>% 
  dplyr::select(-c(Sample_weight_initial_g,mean_wet_grav,mass_sed,mass_water_initial,mass_water_final_g))

#Use this for wet/dry correlation matrix
average_grav <- all_grav_final %>% 
  group_by(kit_treat) %>% 
  mutate(average_grav_intial = mean(grav_dry_initial)) %>% 
  mutate(average_grav_final = mean(grav_dry_final)) %>% 
  mutate(average_grav_lost_subt = average_grav_intial - average_grav_final) %>% 
  mutate(average_grav_lost = mean(lost_grav_perc)) %>% 
  distinct(kit_treat, .keep_all = TRUE) %>% 
  dplyr::select(c(kit_treat,average_grav_intial,average_grav_final,average_grav_lost))

#use this for effect difference correlation matrix
average_grav_lost <- average_grav %>% 
  separate(col = kit_treat, into = c("kit", "Treat"), sep = "_") %>% 
  group_by(kit) %>% 
  mutate(Final_Gravimetric_Moisture_Difference = (average_grav_final[Treat == "Wet"] - average_grav_final[Treat == "Dry"])) %>% 
  dplyr::select(c(kit, Final_Gravimetric_Moisture_Difference)) %>% 
  distinct(kit, .keep_all = TRUE)

##Individual Samples Correlation Matrix

all_list <- list(fe_all, resp, chem_all, all_grav_ind)

all_samples <- all_list %>% 
  reduce(merge, by = "Sample_ID")%>% 
  dplyr::select(-c(kit_treat.x,kit.x,Treat.x,kit_treat.y, kit.y, Treat.y, Log_Fe_mg_kg, ECA,rep, cv_before_removal, cv_after_removal,flag, rate_mg_per_L_per_min, rate_mg_per_L_per_h, Mean_Slope_All, Slope_Removed_Mean, slope_of_the_regression))

all_samples_grn <- merge(all_samples, grn, by = "kit")

all_samples_clean <- all_samples_grn %>% 
  na.omit()  %>% 
  remove_rownames %>% 
  column_to_rownames(var = "Sample_ID") %>% 
  rename(`Rate (mg/L)` = rate_mg) %>% 
  rename(`Initial Gravimetric Water` = grav_dry_initial) %>% 
  rename(`Final Gravimetric Water` = grav_dry_final) %>% 
  rename(`Lost Gravimetric Water` = lost_grav_perc) %>% 
  rename(`Fe (II) (mg/kg)`  = Fe_mg_kg) %>% 
  rename(`Sp. Conducitivity` = SpC) %>% 
  rename(`pH` = pH) %>% 
  rename(`% Fine Sand` = Percent_Fine_Sand) %>% 
  rename(`% Med. Sand` = Percent_Med_Sand) %>% 
  rename(`% Coarse Sand` = Percent_Coarse_Sand) %>% 
  rename(`% Mud` = Percent_Mud) %>% 
  separate(kit_treat, into = c("kit", "Treat"), sep = "_") 

%>% 
  dplyr::select(-c("Percent_Tot_Sand", "Percent_Silt", "Percent_Clay","kit","Treat","Temp"))

"#D55E00" "#0072B2"
###EGU figures

all_samples_clean$Treat <- as.factor(all_samples_clean$Treat)

color_pallete <- colorRampPalette(colors = c("#D55E00", "#0072B2"))

num_colors <- nlevels(all_samples_clean$Treat)

samples_colors <- color_pallete(num_colors)

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/ESS-PI_EGU/", as.character(Sys.Date()),"_Respiration_vs_Mud_lowess.png"), width = 8, height = 8, units = "in", res = 300)

par(mar = c(5 ,6 , 6, 6), xpd = FALSE)

plot(all_samples_clean$`% Mud`, all_samples_clean$`Rate (mg/L)`, cex= 1.8, cex.lab = 2.2, cex.axis = 1.8 ,xlab = "% Mud" , ylab = "Rate (mg/L)", lwd = 2, col = samples_colors[all_samples_clean$Treat], pch = c(16,17)[as.numeric(all_samples_clean$Treat)])

lines(lowess(all_samples_clean$`% Mud`[all_samples_clean$Treat=="Dry"], all_samples_clean$`Rate (mg/L)`[all_samples_clean$Treat=="Dry"]), col = "#D55E00", lwd = 3)


lines(lowess(all_samples_clean$`% Mud`[all_samples_clean$Treat=="Wet"], all_samples_clean$`Rate (mg/L)`[all_samples_clean$Treat=="Wet"]), col = "#0072B2", lwd = 3)

par(xpd = TRUE)

legend(x = "topright", inset = c(-0.15,0.45), legend = paste(levels(all_samples_clean$Treat)), col = samples_colors, pch = c(16,17), cex = 1)

dev.off()


##Fe (mg/L) vs. Respiration
png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/ESS-PI_EGU/", as.character(Sys.Date()),"_Respiration_vs_Fe_lowess.png"), width = 8, height = 8, units = "in", res = 300)

par(mar = c(5 ,6 , 6, 6), xpd = FALSE)

plot(all_samples_clean$`Mean_Fe_mg_per_L`, all_samples_clean$`Rate (mg/L)`, cex= 1.8, cex.lab = 2.2, cex.axis = 1.8 ,xlab = "Fe (II) (mg/L)" , ylab = "Rate (mg/L)", lwd = 2, col = samples_colors[all_samples_clean$Treat], pch = c(16,17)[as.numeric(all_samples_clean$Treat)])

lines(lowess(all_samples_clean$`Mean_Fe_mg_per_L`[all_samples_clean$Treat=="Dry"], all_samples_clean$`Rate (mg/L)`[all_samples_clean$Treat=="Dry"]), col = "#D55E00", lwd = 3)


lines(lowess(all_samples_clean$`Mean_Fe_mg_per_L`[all_samples_clean$Treat=="Wet"], all_samples_clean$`Rate (mg/L)`[all_samples_clean$Treat=="Wet"]), col = "#0072B2", lwd = 3)

par(xpd = TRUE)

legend(x = "topright", inset = c(-0.15,0.45), legend = paste(levels(all_samples_clean$Treat)), col = samples_colors,  pch = c(16,17), cex = 1)

dev.off()

##Fe (mg/kg) vs. Respiration

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/ESS-PI_EGU/", as.character(Sys.Date()),"_Respiration_vs_Fe_mg_kg_lowess.png"), width = 8, height = 8, units = "in", res = 300)

par(mar = c(5 ,6 , 6, 6), xpd = FALSE)

plot(all_samples_clean$`Fe (II) (mg/kg)`, all_samples_clean$`Rate (mg/L)`, cex= 1.8, cex.lab = 2.2, cex.axis = 1.8 ,xlab = "Fe (II) (mg/kg)" , ylab = "Rate (mg/L)", lwd = 2, col = samples_colors[all_samples_clean$Treat], pch = 20)

lines(lowess(all_samples_clean$`Mean_Fe_mg_per_L`, all_samples_clean$`Rate (mg/L)`), col = "blue", lwd = 3)

par(xpd = TRUE)

legend(x = "topright", inset = c(-0.15,0.45), legend = paste(levels(all_samples_clean$Treat)), col = samples_colors, pch = 19, cex = 1)

dev.off()

all_samples_corr <- cor(all_samples_clean, method = "spearman")

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/ESS-PI_EGU/", as.character(Sys.Date()),"_All_Samples_Correlation_Matrix.png"), width = 8, height = 8, units = "in", res = 300)

corrplot(all_samples_corr, title = "All Samples Correlation")

dev.off()

all_samples_dry <- all_samples_grn %>%
  na.omit()  %>% 
  filter(!grepl("Wet", kit_treat)) %>% 
  remove_rownames %>% 
  column_to_rownames(var = "Sample_ID") %>% 
  dplyr::select(-c("Percent_Tot_Sand", "Percent_Silt", "Percent_Clay","kit","Treat", "kit_treat"))

all_samples_dry_corr <- cor(all_samples_dry, method = "spearman")

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/ESS-PI_EGU/", as.character(Sys.Date()),"_All_Dry_Samples_Correlation_Matrix.png"), width = 8, height = 8, units = "in", res = 300)

corrplot(all_samples_dry_corr, title = "All Dry Samples Correlation")

dev.off()

all_samples_wet <- all_samples_grn %>%
  na.omit()  %>% 
  filter(!grepl("Dry", kit_treat)) %>% 
  remove_rownames %>% 
  column_to_rownames(var = "Sample_ID") %>% 
  dplyr::select(-c("Percent_Tot_Sand", "Percent_Silt", "Percent_Clay","kit","Treat", "kit_treat"))

all_samples_wet_corr <- cor(all_samples_wet, method = "spearman")

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/ESS-PI_EGU/", as.character(Sys.Date()),"_All_Wet_Samples_Correlation_Matrix.png"), width = 8, height = 8, units = "in", res = 300)

corrplot(all_samples_wet_corr, title = "All Wet Samples Correlation")

dev.off()

##Wet/Dry Correlation Matrices
wd_list <- list(effect_all, mean_fe_treat, mean_chem, average_grav)

#merge all data frames in list
wet_dry <- wd_list %>% 
  reduce(merge, by = "kit_treat") %>% 
  dplyr::select(-c(kit.x,Treat.x, kit.y, Treat.y))

mean_wet_dry <- merge(wet_dry, grn, by = "kit")

mean_wet <- mean_wet_dry %>%
  na.omit()  %>% 
  filter(!grepl("Dry", kit_treat)) %>% 
  remove_rownames %>% 
  column_to_rownames(var = "kit_treat") %>% 
  rename(`Effect Size` = pos_effect) %>% 
  rename(`Mean Rate (mg/L)` = Average_Rate) %>% 
  rename(`Initial Gravimetric Water` = average_grav_intial) %>% 
  rename(`Final Gravimetric Water` = average_grav_final) %>% 
  rename(`Lost Gravimetric Water` = average_grav_lost) %>% 
  rename(`Fe (II) (mg/kg)`  = Mean_Fe_mg_kg) %>% 
  rename(`Sp. Conducitivity` = Mean_Sp_Conductivity) %>% 
  rename(`pH` = Mean_pH) %>% 
  rename(`% Fine Sand` = Percent_Fine_Sand) %>% 
  rename(`% Med. Sand` = Percent_Med_Sand) %>% 
  rename(`% Coarse Sand` = Percent_Coarse_Sand) %>% 
  rename(`% Mud` = Percent_Mud) %>% 
  dplyr::select(-c("Percent_Tot_Sand", "Percent_Silt", "Percent_Clay","kit","Treat","Mean_Temperature"))
 
mean_wet_corr <- cor(mean_wet, method = "spearman")

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/ESS-PI_EGU/", as.character(Sys.Date()),"_Mean_Wet_Treatment_Correlation_Matrix.png"), width = 8, height = 8, units = "in", res = 300)

corrplot(mean_wet_corr, title = "Wet Treatment Correlation", type = "upper", tl.col = "black", tl.cex = 1.6, cl.cex = 1.25)

dev.off()

mean_dry <- mean_wet_dry %>%
  na.omit()  %>% 
  filter(!grepl("Wet", kit_treat)) %>% 
  remove_rownames %>% 
  column_to_rownames(var = "kit_treat") %>% 
  rename(`Effect Size` = pos_effect) %>% 
  rename(`Mean Rate (mg/L)` = Average_Rate) %>% 
  rename(`Initial Gravimetric Water` = average_grav_intial) %>% 
  rename(`Final Gravimetric Water` = average_grav_final) %>% 
  rename(`Lost Gravimetric Water` = average_grav_lost) %>% 
  rename(`Fe (II) (mg/kg)`  = Mean_Fe_mg_kg) %>% 
  rename(`Sp. Conducitivity` = Mean_Sp_Conductivity) %>% 
  rename(`pH` = Mean_pH) %>% 
  rename(`% Fine Sand` = Percent_Fine_Sand) %>% 
  rename(`% Med. Sand` = Percent_Med_Sand) %>% 
  rename(`% Coarse Sand` = Percent_Coarse_Sand) %>% 
  rename(`% Mud` = Percent_Mud) %>% 
  dplyr::select(-c("Percent_Tot_Sand", "Percent_Silt", "Percent_Clay","kit","Treat","Mean_Temperature"))

mean_dry_corr <- cor(mean_dry, method = "spearman")

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/ESS-PI_EGU/", as.character(Sys.Date()),"_Mean_Dry_Treatment_Correlation_Matrix.png"), width = 8, height = 8, units = "in", res = 300)

corrplot(mean_dry_corr, title = "Dry Treatment Correlation", type = "upper", tl.col = "black", tl.cex = 1.6, cl.cex = 1.25)

dev.off()

##Effect Differences Correlation Matrix

effect_list <- list(effect_diff, mean_fe_diff, mean_chem_diff, average_grav_lost,grn)

#merge all data frames in list
effect <- effect_list %>% 
  reduce(merge, by = "kit") %>% 
  rename(`Effect Size` = pos_effect) %>% 
  rename(`Fe (II) Diff.`  = Fe_Difference) %>% 
  rename(`SpC Diff.` = Sp_Conductivity_Difference) %>% 
  rename(`Temp. Diff.` = Temp_Difference) %>% 
  rename(`pH Diff.` = pH_Difference) %>% 
  rename(`Moisture Diff.` = Final_Gravimetric_Moisture_Difference) %>% 
  rename(`% Fine Sand` = Percent_Fine_Sand) %>% 
  rename(`% Med. Sand` = Percent_Med_Sand) %>% 
  rename(`% Coarse Sand` = Percent_Coarse_Sand) %>% 
  rename(`% Mud` = Percent_Mud) %>% 
  na.omit()  %>% 
  remove_rownames %>% 
  column_to_rownames(var = "kit") %>% 
  dplyr::select(-c("Percent_Tot_Sand", "Percent_Silt", "Percent_Clay"))
  

effect_corr <- cor(effect, method = "spearman")

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/ESS-PI_EGU/", as.character(Sys.Date()),"_Effect_Difference_Correlation_Matrix.png"), width = 8, height = 8, units = "in", res = 300)

corrplot(effect_corr, type = 'upper', tl.col = "black", tl.cex = 1.6, cl.cex = 1.25)

dev.off()

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/ESS-PI_EGU/", as.character(Sys.Date()),"_Effect_vs_Mud_Scatter.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(effect, aes(x = Percent_Fine_Sand, y = Positive_Effect_Size))+
  geom_point()+
  geom_smooth(method = lm)

dev.off()

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/ESS-PI_EGU/", as.character(Sys.Date()),"_Effect_vs_Mud_Lowess.png"), width = 8, height = 8, units = "in", res = 300)

par(mar = c(5 ,6 , 4, 1))

plot(effect$`% Mud`, effect$`Effect Size`, cex= 1.8, cex.lab = 2.2, cex.axis = 1.8 ,xlab = "% Mud" , ylab = "Effect Size", lwd = 2)

lines(lowess(effect$`% Mud`, effect$`Effect Size`), col = "blue", lwd = 3)

dev.off()

###PCA/RDA

mean_wet_pca <- princomp(mean_wet_corr)
summary(mean_wet_pca)

mean_wet_pca$loadings[, 1:2]

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/ESS-PI_EGU/", as.character(Sys.Date()),"_Wet_Mean_PCA.png"), width = 8, height = 8, units = "in", res = 300)

fviz_pca_var(mean_wet_pca, col.var = "black")

dev.off()


mean_dry_pca <- princomp(mean_dry_corr)
summary(mean_dry_pca)

mean_dry_pca$loadings[, 1:2]

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/ESS-PI_EGU/", as.character(Sys.Date()),"_Dry_Mean_PCA.png"), width = 8, height = 8, units = "in", res = 300)

fviz_pca_var(mean_dry_pca, col.var = "black")

dev.off()

mean_wet_rda <- rda(mean_wet$`Mean Rate (mg/L)` ~., data = mean_wet)


summary(mean_wet_rda)

fwd_wet <- ordiR2step(rda(mean_wet$`Mean Rate (mg/L)` ~ 1, data = mean_wet), 
          scope = formula(mean_wet_rda), 
          direction = "forward", 
          R2scope = TRUE,
          pstep = 1000,
          trace = FALSE)

fwd_wet$call

mean_wet_signif <- rda(formula = mean_wet$`Mean Rate (mg/L)` ~ `% Coarse Sand` + `Initial Gravimetric Water`, data = mean_wet)

RsquareAdj(mean_wet_signif)

ordiplot(mean_wet_signif, scaling = 1, type = "text")

ordiplot(mean_wet_rda, scaling = 1, type = "text")
