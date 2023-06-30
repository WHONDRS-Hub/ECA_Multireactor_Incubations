##ECA Data Analysis 

#### Read in Data ####

#read in libraries

library(lubridate);library(writexl);library(raster);library(tidyverse);library(devtools)
library(readxl)
library(corrplot)
library(corrr)
library(vegan)
library(FactoMineR)
library(factoextra)

rm(list=ls());graphics.off()

# Set working directory to data file
#Example:
pnnl.user = 'laan208'

# choose file dates to read in 

effect.date = '2023-06-28'
respiration.date = '2023-06-28'
grav.date = '2023-06-27'

#Read in all data
setwd(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/"))

#effect size - change date to most recent
effect_size <- read_csv(paste0("Optode multi reactor/Optode_multi_reactor_incubation/effect size/Effect_Size_Data/Effect_Size_merged_by_laan208_on_",effect.date,".csv"))


#Respiration rates with removals from dist matrix to calculate effect size 
respiration <- read.csv(paste0("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/rates/Plots/All_Respiration_Rates/removed_respiration_merged_by_laan208_on_",respiration.date,".csv"))

#ECA Iron
iron <- read_csv(paste0("FE/03_ProcessedData/EC_SFE_ReadyForBoye_06-29-2023.csv"))


#ICON Grain Size
grain <- read_csv("C:/GitHub/ECA_Multireactor_Incubations/Data/v2_CM_SSS_Sediment_Grain_Size.csv")


#All incubation pH, SpC, temp
chemistry <- paste0("Optode multi reactor/Optode_multi_reactor_incubation/")

import_data = function(chemistry){ 
  
  #pull a list of files in target folder with correct pattern
  #read all files and combine
  
  map.file <-  list.files(chemistry, recursive = T, pattern = "\\.xlsx$",full.names = T)
  
  map.file <- map.file[!grepl("red|Red|EC_01_|EC_02_|EC_03_|EC_06_07|EC_10_15|EC_04_08|fast|combined|SPC|IC|QA", map.file)]
  
  mapping <- lapply(map.file, read_xlsx)
  
  for (i in 1:length(mapping)){mapping[[i]] <- cbind(mapping[[i]], map.file[i])}
  
  all.map <- 
    do.call(rbind,mapping)
}

map = import_data(chemistry)

#Gravimetric Moisture

grav_inc <- read.csv(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/INC/03_ProcessedData/ECA_Drying_Masses_merged_by_laan208_on_",grav.date,".csv"))

wet_wt <- read.csv(paste0("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/MOI/03_ProcessedData/EC_Moisture_Content_2022.csv"))



#### EFFECT SIZE/RESPIRATION #### 

#Use this for individual samples correlation matrix

resp <- respiration %>% 
  separate(Sample_ID, into = c("EC", "kit", "INC"), sep = "_", remove = FALSE) %>% 
  mutate(Treat = case_when(grepl("W",INC)~"Wet",
                           grepl("D", INC) ~"Dry")) %>% 
  relocate(Treat, .after = kit) %>% 
  dplyr::select(-c(EC, INC))

#Use this for wet/dry correlation matrix
effect_all <- effect_size %>%
  mutate(Average_Rate = abs(Mean_Slope_Removed)) %>%
  rename(kit = Site) %>% 
  dplyr::select(c(kit, Treat,Average_Rate,effect,log_effect))
 
#Use this for effect size correlation matrix
effect_diff <- effect_all %>% 
  distinct(kit, .keep_all = TRUE) %>% 
  dplyr::select(-c(Average_Rate,Treat))


#### IRON ####

#Individual Fe Samples after averaging for analytical reps

#calculate mean Fe for kit/treatment from analytical replicates

fe_all <- iron %>% 
  separate(Sample_Name, into = c("Sample_ID", "rep"), sep = -1, convert = TRUE) %>% 
  group_by(Sample_ID) %>% 
  mutate(Mean_Rep_Fe_mg_per_L = mean(Fe_mg_per_L)) %>% 
  #mutate(Mean_Rep_Fe_mg_kg= mean(Fe_mg_per_kg_sediment)) %>% 
  separate(Sample_ID, c("ECA", "kit", "rep"), sep = "_") %>% 
  mutate(Treat = case_when(grepl("W",rep)~"Wet",
                           grepl("D", rep) ~"Dry")) %>% 
  relocate(Treat, .after = rep) %>% 
 # mutate(Log_Mean_Rep_Fe_mg_kg = log10(Mean_Rep_Fe_mg_kg + 1)) %>% 
  mutate(Log_Mean_Rep_Fe_mg_L = log10(Mean_Rep_Fe_mg_per_L + 1)) %>% 
  mutate(rep = str_replace(rep, "SFE", "INC")) %>% 
  na.omit()

for (i in 1:nrow(fe_all)){
  
  if (str_count(fe_all$kit[i], "[0-9]") <= 2){
    
    fe_all$kit[i] = paste0("0", fe_all$kit[i])
    
  }
  
  else {
    
    fe_all$kit[i] = fe_all$kit[i]
  }
  
}


fe_all <- fe_all %>% 
  unite(Sample_ID, c("ECA", "kit", "rep"), sep = "_", remove = FALSE) %>% 
  mutate(Sample_ID = str_replace(Sample_ID, "SFE", "INC")) %>% 
  distinct(Sample_ID, .keep_all = TRUE) %>% 
  dplyr::select(-c(Fe_mg_per_L,#Fe_mg_per_kg_sediment, 
                   rep, ECA, Methods_Deviation))


##Check CV for Mean_Rep_Fe_mg_kg for all reps (W/D)
mean_fe_check <- fe_all %>% 
  group_by(kit, Treat) %>% 
  mutate(CV = (sd(Mean_Rep_Fe_mg_per_L)/mean(Mean_Rep_Fe_mg_per_L))*100)
  #mutate(CV = (sd(Mean_Rep_Fe_mg_kg)/mean(Mean_Rep_Fe_mg_kg))*100)

mean_fe_check$Mean_Rep_Fe_mg_kg <- format(mean_fe_check$Mean_Rep_Fe_mg_kg, scientific = FALSE)

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

# ggplot(mean_fe, aes(y = Fe_mg_kg, x = kit_treat))+
#     geom_boxplot()+
#     theme_bw()+
#     theme(axis.title.x = element_text(size = 18),
#           axis.title.y = element_text(size = 18),
#           axis.text.x = element_text(size = 16),
#           axis.text.y = element_text(size = 16))+
#     xlab("\n Fe (II) (mg/L)")+
#     ylab("Count\n")
  

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


#Fe Means for use in wet/dry correlation matrix
mean_fe_treat <- fe_all %>% 
  group_by(kit, Treat) %>%
  mutate(Mean_Treat_Fe_mg_L = mean(Mean_Rep_Fe_mg_per_L)) %>% 
 # mutate(Mean_Treat_Fe_mg_kg = mean(Mean_Rep_Fe_mg_kg)) %>% 
  mutate(Log_Mean_Treat_Fe_mg_L = log10(Mean_Treat_Fe_mg_L)) %>% 
 # mutate(Log_Mean_Treat_Fe_mg_kg = log10(Mean_Treat_Fe_mg_kg)) %>% 
  dplyr::select(c(kit,Treat,Mean_Treat_Fe_mg_L,#Mean_Treat_Fe_mg_kg,Log_Mean_Treat_Fe_mg_kg
                  )) %>% 
  distinct(.keep_all = TRUE)

#Differences for Effect Size Correlation Matrix
mean_fe_diff <- mean_fe_treat %>% 
  group_by(kit) %>% 
  mutate(Fe_Difference_mg_L = (Mean_Treat_Fe_mg_L[Treat == "Wet"] - Mean_Treat_Fe_mg_L[Treat == "Dry"])) %>% 
  # mutate(Fe_Difference_mg_kg = (Mean_Treat_Fe_mg_kg[Treat == "Wet"] - Mean_Treat_Fe_mg_kg[Treat == "Dry"])) %>% 
  dplyr::select(-c(Mean_Treat_Fe_mg_L, #Mean_Treat_Fe_mg_kg, Log_Mean_Treat_Fe_mg_kg, 
                   Treat)) %>% 
  distinct(kit, .keep_all = TRUE)
 

#### GRAIN SIZE ####

#this is the only data type that is not vial specific

grain2 <- grain[-c(1,3:13, 133), -1]
names(grain2) <- grain2[1,]
grain2 <- grain2[-1,]


grn <- grain2 %>% 
  filter(!grepl("SSS", Sample_Name)) %>% 
  separate(Sample_Name, c("CM", "kit", "an"), sep = "_", remove = FALSE) %>%
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
  filter(kit != "025") %>% 
  filter(kit != "026") %>% 
  filter(kit != "029") %>% 
  filter(kit != "030") %>% 
  filter(kit != "028") %>% 
  filter(kit != "043") %>% 
  filter(kit != "050") %>% 
  filter(kit != "062") %>% 
  mutate_at(c("Percent_Coarse_Sand", "Percent_Med_Sand", "Percent_Fine_Sand", "Percent_Silt", "Percent_Clay", "Percent_Tot_Sand"), as.numeric) %>% 
  mutate(Percent_Mud = Percent_Clay + Percent_Silt) %>%
  dplyr::select(-c(CM, an, Sample_Name, Methods_Deviation, Material)) 

#grn$kit <- sub('.', '', grn$kit)

grn <- grn %>% 
  mutate(geom_rusle = exp(0.01 * ((Percent_Coarse_Sand* log(1.25)) + (Percent_Med_Sand* log(0.375)) + (Percent_Fine_Sand * log(0.1)) + (Percent_Silt * log(0.026)) + (Percent_Clay * log(0.01))))) %>% 
  mutate(geom = ((Percent_Coarse_Sand/100)*1.25) + ((Percent_Med_Sand/100) * 0.375) + ((Percent_Fine_Sand/100) *0.1) + ((Percent_Silt/100) * 0.026) + ((Percent_Clay/100 * 0.01))) 

grn_all <- grn %>% 
  dplyr::select(-c(geom_rusle, geom, Percent_Tot_Sand, Percent_Mud)) %>% 
  mutate(Percent_Coarse_Sand_Finer = Percent_Coarse_Sand + Percent_Med_Sand + Percent_Fine_Sand + Percent_Silt + Percent_Clay) %>% 
  mutate(Percent_Med_Sand_Finer = Percent_Med_Sand + Percent_Fine_Sand + Percent_Silt + Percent_Clay) %>% 
  mutate(Percent_Fine_Sand_Finer = Percent_Fine_Sand + Percent_Silt + Percent_Clay) %>% 
  mutate(Percent_Silt_Finer = Percent_Silt + Percent_Clay) %>% 
  mutate(Percent_Clay_Finer = Percent_Clay) %>% 
  dplyr::select(-c(Percent_Coarse_Sand, Percent_Med_Sand, Percent_Fine_Sand, Percent_Silt, Percent_Clay)) %>% 
  pivot_longer(cols = c("Percent_Coarse_Sand_Finer", "Percent_Med_Sand_Finer", "Percent_Fine_Sand_Finer", "Percent_Silt_Finer", "Percent_Clay_Finer"), names_to = "Category", values_to = "fraction") %>% 
  mutate(size = case_when(
    grepl("Coarse", Category) ~ 1.25,
    grepl("Med", Category) ~ 0.375,
    grepl("Fine_Sand", Category) ~ 0.1,
    grepl("Silt", Category) ~ 0.026, 
    grepl("Clay", Category) ~ 0.01)) %>% 
  mutate(y_int = 50) %>% 
  group_by(kit) %>% 
  mutate(slope = if_else(fraction[Category == "Percent_Med_Sand_Finer"] < 50, ((fraction[Category == "Percent_Coarse_Sand_Finer"] - fraction[Category == "Percent_Med_Sand_Finer"]) / (size[Category == "Percent_Coarse_Sand_Finer"] - size[Category == "Percent_Med_Sand_Finer"])), 
                if_else(fraction[Category == "Percent_Fine_Sand_Finer"] < 50, ((fraction[Category == "Percent_Med_Sand_Finer"] - fraction[Category == "Percent_Fine_Sand_Finer"]) / (size[Category == "Percent_Med_Sand_Finer"] - size[Category == "Percent_Fine_Sand_Finer"])), 
                if_else(fraction[Category == "Percent_Silt_Finer"] < 50, ((fraction[Category == "Percent_Fine_Sand_Finer"] - fraction[Category == "Percent_Silt_Finer"]) / (size[Category == "Percent_Fine_Sand_Finer"] - size[Category == "Percent_Silt_Finer"])), ((fraction[Category == "Percent_Silt_Finer"] - fraction[Category == "Percent_Clay_Finer"]) / (size[Category == "Percent_Silt_Finer"] - size[Category == "Percent_Clay_Finer"]))
    )))) %>% 
  mutate(x_int = if_else(fraction[Category == "Percent_Med_Sand_Finer"] < 50, (((50 - fraction[Category == "Percent_Med_Sand_Finer"])/slope[Category == "Percent_Coarse_Sand_Finer"])+size[Category == "Percent_Med_Sand_Finer"]),
                  if_else(fraction[Category == "Percent_Fine_Sand_Finer"] < 50, (((50 - fraction[Category == "Percent_Fine_Sand_Finer"])/slope[Category == "Percent_Med_Sand_Finer"])+size[Category == "Percent_Fine_Sand_Finer"]),
                  if_else(fraction[Category == "Percent_Silt_Finer"] < 50,(((50 - fraction[Category == "Percent_Silt_Finer"])/slope[Category == "Percent_Fine_Sand_Finer"])+size[Category == "Percent_Silt_Finer"]), (((50 - fraction[Category == "Percent_Clay_Finer"])/slope[Category == "Percent_Silt_Finer"])+size[Category == "Percent_Clay_Finer"])))))
  
#pdf(file = "C:/Users/laan208/OneDrive - PNNL/Documents/d50_plots/my_plots.pdf")

#unique.samples = unique(grn_all$kit)

# for (i in 1:length(unique.samples)){
# 
#   site_subset = subset(grn_all, grn_all$kit == unique.samples[i])
#   
# kit <- ggplot(site_subset, aes(x = size, y = fraction))+
#         geom_point()+
#         geom_line()+
#         geom_hline(yintercept = 50, color = "red", linetype = "dashed")+
#         geom_vline(xintercept = site_subset$x_int, color = "red", linetype = "dashed") +
#   ggtitle(site_subset$kit[1])
# 
# print(kit)
# 
# }
# 
# dev.off()

d50 <- grn_all %>% 
  distinct(kit, .keep_all = TRUE) %>% 
  dplyr::select(c(kit, x_int)) %>% 
  rename(d50 = x_int)

grn <- merge(d50, grn)

grn <- grn %>% 
  relocate(d50, .after = geom)

# geom <- grn %>% 
#   select(c(kit, geom_rusle,geom)) %>% 
#   rename(log_geom_mean = geom_rusle) %>% 
#   rename(geom_mean = geom)
# 
# all_geom <- left_join(d50, geom, by = c("kit"))
# 
# geom_vs_d50 <- ggplot(all_geom, aes(x = geom_mean, y = d50))+
#   geom_point()+
#   geom_smooth()
# 
# geom_vs_rusle <- ggplot(all_geom, aes(x = geom_mean, y = log_geom_mean))+
#   geom_point()+
#   geom_smooth()
# 
# d50_vs_rusle <- ggplot(all_geom, aes(x = d50, y = log_geom_mean))+
#   geom_point()+
#   geom_smooth()


#Coarse sand = 0.5 - 2 mm (1.25 mm), med sand = 0.25 - 0.499 (0.375) , fine sand = 0.05 - 0.25 (0.1), silt = , clay =


# png(file = paste0("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/", as.character(Sys.Date()),"_GRN_stacked.png"), width = 8, height = 8, units = "in", res = 300)
# 
# ggplot(grn_long, aes(fill = size, y = percent, x = Sample_ID)) + 
#   geom_bar(position="fill", stat="identity")+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust =1))
# 
# dev.off()


#### pH, SpC, Temp ####

#Use this for individual correlation matrix
map_corr = map %>% 
  separate(Sample_ID, c("ECA"
                        , "kit", "Analysis"), sep = "_", remove = FALSE)%>% 
  filter(!grepl("Blank", Sample_ID))

for (i in 1:nrow(map_corr)){
  
  if (str_count(map_corr$kit[i], "[0-9]") <= 2){
    
    map_corr$kit[i] = paste0("0", map_corr$kit[i])
    
  }
  
  else {
    
   map_corr$kit[i] = map_corr$kit[i]
  }
  
}

map_corr <- map_corr %>% 
  unite(Sample_ID, c("ECA", "kit", "Analysis"), sep = "_")

chem_all = map_corr %>% 
  separate(Sample_ID, c("ECA"
, "kit", "Analysis"), sep = "_", remove = FALSE) %>% 
  mutate(Treat = case_when(grepl("W",Analysis)~"Wet",
                           grepl("D", Analysis) ~"Dry")) %>% 
  dplyr::select(-`Disk Calibration date`, -Location, -`Time on`, -`Time off`, -Notes, -`Disk_ID`, -`map.file[i]`, -`ECA`, -`Analysis`) %>% 
  filter(!grepl("Blank", Sample_ID))

#Use this for wet/dry correlation matrices
mean_chem <- chem_all %>% 
  group_by(kit, Treat) %>% 
  mutate(Mean_Sp_Conductivity = mean(SpC)) %>% 
  mutate(Mean_Temperature = mean(Temp)) %>% 
  mutate(Mean_pH = mean(pH)) %>% 
  dplyr::select(c(kit, Treat,Mean_Sp_Conductivity, Mean_Temperature, Mean_pH)) %>% 
  distinct(.keep_all = TRUE)

#Use this for effect size matrix
mean_chem_diff <- mean_chem %>% 
  group_by(kit) %>% 
  mutate(Sp_Conductivity_Difference = (Mean_Sp_Conductivity[Treat == "Wet"] - Mean_Sp_Conductivity[Treat == "Dry"])) %>% 
  mutate(Temp_Difference = (Mean_Temperature[Treat == "Wet"] - Mean_Temperature[Treat == "Dry"])) %>% 
  mutate(pH_Difference = (Mean_pH[Treat == "Wet"] - Mean_pH[Treat == "Dry"])) %>% 
  dplyr::select(-c(Mean_Sp_Conductivity,Mean_Temperature, Mean_pH,Treat)) %>% 
  distinct(.keep_all = TRUE)  

#### Gravimetric Moisture ####

mean_wet_wt <- wet_wt %>% 
  separate(col = sample_name, into = c("Project", "kit", "analysis"), sep = "_") %>%  
  separate(col = analysis, into = c("Analysis", "Replicate"), sep = "-") %>% 
  group_by(kit) %>% 
  summarise(mean_wet_grav = (mean(percent_water_content_wet)/100))

for (i in 1:nrow(mean_wet_wt)){
  
  if (str_count(mean_wet_wt$kit[i], "[0-9]") <= 2){
    
    mean_wet_wt$kit[i] = paste0("0", mean_wet_wt$kit[i])
    
  }
  
  else {
    
    mean_wet_wt$kit[i] = mean_wet_wt$kit[i]
  }
  
}

grav <- grav_inc %>% 
  group_by(Sample_Name) %>% 
  mutate(Sample_weight_initial_g = first(Sample_weight_g)) %>% 
  mutate(Sample_weight_final_g = last(Sample_weight_g)) %>% 
  ungroup() %>% 
  separate(col = Sample_Name, into = c("Project", "kit", "analysis"), sep = "_", remove = FALSE) %>%  
  separate(col = analysis, into = c("Analysis", "Replicate"), sep = "-") %>%
  mutate(Treat = case_when(grepl("W",Replicate)~"Wet",
                           grepl("D", Replicate) ~"Dry")) %>% 
  relocate(Treat, .after = Sample_Name) 

for (i in 1:nrow(grav)){
  
  if (str_count(grav$kit[i], "[0-9]") <= 2){
    
    grav$kit[i] = paste0("0", grav$kit[i])
    
  }
  
  else {
    
    grav$kit[i] = grav$kit[i]
  }
  
}

grav <- grav %>% 
  unite(Sample_Name, c("Project", "kit", "Analysis"), sep = "_", remove = FALSE) %>% 
  unite(Sample_Name, c("Sample_Name", "Replicate"), sep = "-") %>% 
  select(-c(Project, Analysis))
  

final_grav <- merge(grav, mean_wet_wt, by = "kit")

all_grav <- final_grav %>% 
  mutate(mass_sed = (Sample_weight_initial_g - (Sample_weight_initial_g * mean_wet_grav))) %>% 
  mutate(mass_water_initial = (Sample_weight_initial_g * mean_wet_grav) + Water_added_initial_g) %>%  
  mutate(grav_dry_initial = mass_water_initial/mass_sed) %>% 
  mutate(mass_water_final_g = (Sample_weight_final_g - mass_sed)+Water_added_initial_g) %>% 
  mutate(grav_dry_final = mass_water_final_g/mass_sed) %>% 
  mutate(lost_grav_perc = grav_dry_initial - grav_dry_final)

#Use this for individual correlation matrix
all_grav_ind <- all_grav %>% 
  rename(Sample_ID = Sample_Name) %>% 
  distinct(Sample_ID, .keep_all = TRUE) %>% 
  dplyr::select(-c(Sample_weight_initial_g,mean_wet_grav,mass_sed,mass_water_initial,mass_water_final_g,Date, Sample_weight_g,Water_added_initial_g,Sample_weight_final_g))

#Use this for wet/dry correlation matrix
average_grav <- all_grav_ind %>% 
  group_by(kit, Treat) %>% 
  mutate(average_grav_intial = mean(grav_dry_initial)) %>% 
  mutate(average_grav_final = mean(grav_dry_final)) %>% 
  mutate(average_grav_lost_subt = average_grav_intial - average_grav_final) %>% 
  mutate(average_grav_lost = mean(lost_grav_perc)) %>% 
  dplyr::select(c(kit, Treat,average_grav_intial,average_grav_final,average_grav_lost)) %>% 
  distinct(.keep_all = TRUE)

#use this for effect difference correlation matrix
average_grav_lost <- average_grav %>% 
  group_by(kit) %>% 
  mutate(Final_Gravimetric_Moisture_Difference = (average_grav_final[Treat == "Wet"] - average_grav_final[Treat == "Dry"])) %>% 
  dplyr::select(c(kit, Final_Gravimetric_Moisture_Difference)) %>% 
  distinct(.keep_all = TRUE)

#### Individual Samples Correlation Matrix ####

all_list <- list(fe_all, resp, chem_all, all_grav_ind)

all_samples <- all_list %>% 
  reduce(merge, by = c("Sample_ID", "kit", "Treat"), all = TRUE)%>% 
  dplyr::select(-c(#Mean_Rep_Fe_mg_per_L,
    Sample_weight_Fill_g,
    rate_mg_per_L_per_h,
    #Log_Mean_Rep_Fe_mg_kg
                   ))
 
all_samples_grn <- merge(all_samples, grn, by = "kit", all = TRUE)

all_samples_clean <- all_samples_grn %>% 
  na.omit()  %>% 
  #filter(!is.na(Sample_ID)) %>% 
  remove_rownames %>% 
  column_to_rownames(var = "Sample_ID") %>% 
  rename(`Rate (mg/L)` = rate_mg_per_L_per_min) %>% 
  rename(`Initial Gravimetric Water` = grav_dry_initial) %>% 
  rename(`Final Gravimetric Water` = grav_dry_final) %>% 
  rename(`Lost Gravimetric Water` = lost_grav_perc) %>% 
 # rename(`Fe (II) (mg/kg)`  = Mean_Rep_Fe_mg_kg) %>% 
  rename(`Fe (II) (mg/L)` = Mean_Rep_Fe_mg_per_L) %>% 
  rename(`Sp. Conducitivity` = SpC) %>% 
  rename(`pH` = pH) %>% 
  rename(`% Fine Sand` = Percent_Fine_Sand) %>% 
  rename(`% Med. Sand` = Percent_Med_Sand) %>% 
  rename(`% Coarse Sand` = Percent_Coarse_Sand) %>% 
  rename(`% Mud` = Percent_Mud) %>% 
  rename(`Geometric Mean` = geom) %>% 
  rename(`RUSLE Geometric Mean` = geom_rusle) %>% 
  rename(`D50` = d50)


#### EGU figures ####

all_samples_clean$Treat <- as.factor(all_samples_clean$Treat)

color_pallete <- colorRampPalette(colors = c("#D55E00", "#0072B2"))

num_colors <- nlevels(all_samples_clean$Treat)

samples_colors <- color_pallete(num_colors)

## Respiration vs. Mud ####
png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/ESS-PI_EGU/", as.character(Sys.Date()),"_Respiration_vs_Mud_lowess.png"), width = 8, height = 8, units = "in", res = 300)

par(mar = c(5 ,6 , 6, 6), xpd = FALSE)

plot(all_samples_clean$`% Mud`, all_samples_clean$`Rate (mg/L)`, cex= 1.8, cex.lab = 2.2, cex.axis = 1.8 ,xlab = "% Mud" , ylab = expression("Respiration Rate (mg O"[2]*" L"^- 1*" min"^-1*")"), lwd = 2, col = samples_colors[all_samples_clean$Treat], pch = c(16,17)[as.numeric(all_samples_clean$Treat)])

#lines(lowess(all_samples_clean$`% Mud`[all_samples_clean$Treat=="Dry"], all_samples_clean$`Rate (mg/L)`[all_samples_clean$Treat=="Dry"]), col = "#D55E00", lwd = 3)


#lines(lowess(all_samples_clean$`% Mud`[all_samples_clean$Treat=="Wet"], all_samples_clean$`Rate (mg/L)`[all_samples_clean$Treat=="Wet"]), col = "#0072B2", lwd = 3)

par(xpd = TRUE)

legend(x = "topright", inset = c(-0.15,0.45), legend = paste(levels(all_samples_clean$Treat)), col = samples_colors, pch = c(16,17), cex = 1)

dev.off()

## Fe (mg/kg) vs. Respiration ####

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/ESS-PI_EGU/", as.character(Sys.Date()),"_Respiration_vs_Fe_mg_kg_lowess.png"), width = 8, height = 8, units = "in", res = 300)

par(mar = c(5 ,6 , 6, 6), xpd = FALSE)

plot(all_samples_clean$`Fe (II) (mg/L)`, all_samples_clean$`Rate (mg/L)`, cex= 1.8, cex.lab = 2.2, cex.axis = 1.8 ,xlab = expression("Fe (II) (mg kg"^-1*")"), ylab = expression("Respiration Rate (mg O"[2]*" L"^- 1*" min"^-1*")"), lwd = 2, col = samples_colors[all_samples_clean$Treat], pch = c(16,17)[as.numeric(all_samples_clean$Treat)])

#lines(lowess(all_samples_clean$`Mean_Fe_mg_per_L`, all_samples_clean$`Rate (mg/L)`), col = "blue", lwd = 3)

par(xpd = TRUE)

legend(x = "topright", inset = c(-0.15,0.45), legend = paste(levels(all_samples_clean$Treat)), col = samples_colors,  pch = c(16,17), cex = 1)

dev.off()

## Correlation Ind ####

all_samples_clean_corr <- all_samples_clean %>% 
  dplyr::select(-c(kit, Treat, ID, rep, kit_treat))

all_samples_corr <- cor(all_samples_clean_corr, method = "spearman")

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/", as.character(Sys.Date()),"_All_Samples_Correlation_Matrix.png"), width = 8, height = 8, units = "in", res = 300)

corrplot(all_samples_corr,type = "upper", tl.col = "black", tl.cex = 1.6, cl.cex = 1.25,  title = "All Samples Correlation")

dev.off()

## Dry Corr Ind ####

all_samples_dry <- all_samples_clean %>%
  na.omit()  %>% 
  filter(!grepl("Wet", Treat)) %>% 
  dplyr::select(-c("Percent_Tot_Sand", "Percent_Silt", "Percent_Clay", "kit", "Treat", "% Fine Sand", "% Med. Sand", "% Coarse Sand", "% Mud", "RUSLE Geometric Mean", "Geometric Mean", "ID", "rep", "kit_treat"))

all_samples_dry_corr <- cor(all_samples_dry,method = "spearman")

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/", as.character(Sys.Date()),"_All_Dry_Samples_Correlation_Matrix.png"), width = 8, height = 8, units = "in", res = 300)

corrplot(all_samples_dry_corr, type = "upper", tl.col = "black", tl.cex = 1.6, cl.cex = 1.25, title = "All Dry Samples Correlation")

dev.off()

## Wet Correlation Ind ####
all_samples_wet <- all_samples_clean %>%
  na.omit()  %>% 
  filter(!grepl("Dry", Treat)) %>% 
  dplyr::select(-c("Percent_Tot_Sand", "Percent_Silt", "Percent_Clay","kit", "Treat", "% Fine Sand", "% Med. Sand", "% Coarse Sand", "% Mud", "RUSLE Geometric Mean", "Geometric Mean", "ID", "rep", "kit_treat"))

all_samples_wet_corr <- cor(all_samples_wet, method = "spearman")

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/", as.character(Sys.Date()),"_All_Wet_Samples_Correlation_Matrix.png"), width = 8, height = 8, units = "in", res = 300)

corrplot(all_samples_wet_corr,type = "upper", tl.col = "black", tl.cex = 1.6, cl.cex = 1.25, title = "All Wet Samples Correlation")

dev.off()

#### Wet/Dry Correlation Matrices ####
wd_list <- list(effect_all, mean_fe_treat, mean_chem, average_grav)

#merge all data frames in list
wet_dry <- wd_list %>% 
  reduce(merge, by = c("kit", "Treat"), remove = FALSE)

mean_wet_dry <- merge(wet_dry, grn, by = "kit")

mean_wet_dry_clean <- mean_wet_dry %>% 
  rename(`Average Rate (mg/L)` = Average_Rate) %>% 
  rename(`Effect Size` = effect) %>% 
    rename(`Mean Initial Gravimetric Water` = average_grav_intial) %>% 
  rename(`Mean Final Gravimetric Water` = average_grav_final) %>% 
  rename(`Mean Lost Gravimetric Water` = average_grav_lost) %>% 
 # rename(`Mean Fe (II) (mg/kg)`  = Mean_Treat_Fe_mg_kg) %>% 
  rename(`Mean Sp. Conducitivity` = Mean_Sp_Conductivity) %>% 
  rename(`Mean pH` = Mean_pH) %>% 
  rename(`Mean Temp.` = Mean_Temperature) %>% 
  rename(`% Fine Sand` = Percent_Fine_Sand) %>% 
  rename(`% Med. Sand` = Percent_Med_Sand) %>% 
  rename(`% Coarse Sand` = Percent_Coarse_Sand) %>% 
  rename(`% Mud` = Percent_Mud) %>% 
  rename(`Geometric Mean` = geom) %>% 
  rename(`RUSLE Geometric Mean` = geom_rusle) %>% 
  rename(`D50` = d50) %>% 
  dplyr::select(-c(#Mean_Treat_Fe_mg_L,
    #Log_Mean_Treat_Fe_mg_kg,
    `RUSLE Geometric Mean`, `Geometric Mean`, log_effect, Percent_Clay, `% Fine Sand`,`% Med. Sand`, `% Coarse Sand`,`% Mud`, Percent_Silt, Percent_Tot_Sand)) %>% 
  na.omit()

mean_wet_dry_clean_corr <- mean_wet_dry_clean %>% 
  unite(kit_Treat, kit, Treat, sep = "_") %>% 
  na.omit() %>% 
  remove_rownames %>% 
  column_to_rownames(var = c("kit_Treat")) 

mean_samples_corr <- cor(mean_wet_dry_clean_corr, method = "spearman")

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/", as.character(Sys.Date()),"_Mean_Wet_Dry_Samples_Correlation_Matrix.png"), width = 8, height = 8, units = "in", res = 300)

corrplot(mean_samples_corr,type = "upper", tl.col = "black", tl.cex = 1.6, cl.cex = 1.25,  title = "Mean Samples Correlation")

dev.off()

## Mean Wet Corr ####
mean_wet <- mean_wet_dry_clean %>%
  filter(!grepl("Dry", Treat)) %>% 
  remove_rownames %>% 
  column_to_rownames(var = "kit") %>% 
  dplyr::select(-c(Treat))
 
mean_wet_corr <- cor(mean_wet, method = "spearman")

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/", as.character(Sys.Date()),"_Mean_Wet_Treatment_Correlation_Matrix.png"), width = 8, height = 8, units = "in", res = 300)

corrplot(mean_wet_corr, title = "Wet Treatment Correlation", type = "upper", tl.col = "black", tl.cex = 1.6, cl.cex = 1.25)

dev.off()

## Mean Dry Corr ####
mean_dry <-  mean_wet_dry_clean %>%
  filter(!grepl("Wet", Treat)) %>% 
  remove_rownames %>% 
  column_to_rownames(var = "kit") %>% 
  dplyr::select(-c(Treat))

mean_dry_corr <- cor(mean_dry, method = "spearman")

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/", as.character(Sys.Date()),"_Mean_Dry_Treatment_Correlation_Matrix.png"), width = 8, height = 8, units = "in", res = 300)

corrplot(mean_dry_corr, title = "Dry Treatment Correlation", type = "upper", tl.col = "black", tl.cex = 1.6, cl.cex = 1.25)

dev.off()

#### Effect Differences Correlation Matrix ####

effect_list <- list(effect_diff, mean_fe_diff, mean_chem_diff, average_grav_lost,grn)

#merge all data frames in list
effect <- effect_list %>% 
  reduce(merge, by = "kit") %>% 
  rename(`Effect Size` = effect) %>% 
  #rename(`Fe (II) Diff.`  = Fe_Difference_mg_kg) %>% 
  rename(`SpC Diff.` = Sp_Conductivity_Difference) %>% 
  rename(`Temp. Diff.` = Temp_Difference) %>% 
  rename(`pH Diff.` = pH_Difference) %>% 
  rename(`Moisture Diff.` = Final_Gravimetric_Moisture_Difference) %>% 
  rename(`% Fine Sand` = Percent_Fine_Sand) %>% 
  rename(`% Med. Sand` = Percent_Med_Sand) %>% 
  rename(`% Coarse Sand` = Percent_Coarse_Sand) %>% 
  rename(`% Mud` = Percent_Mud) %>% 
  rename(D50 = d50) %>% 
  na.omit()  %>% 
  remove_rownames %>% 
  column_to_rownames(var = "kit") %>% 
  dplyr::select(-c("Percent_Tot_Sand", "Percent_Silt", "Percent_Clay", "log_effect"))# %>% 
 # dplyr::select(-c("% Fine Sand", "% Med. Sand", "% Coarse Sand", "% Mud", "geom_rusle", "geom"))
  
effect_corr <- cor(effect, method = "spearman")

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/", as.character(Sys.Date()),"_Effect_Difference_Correlation_Matrix.png"), width = 8, height = 8, units = "in", res = 300)

corrplot(effect_corr, type = 'upper', tl.col = "black", tl.cex = 1.6, cl.cex = 1.25)

dev.off()

## Fine Sand vs. Effect Size ####

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/ESS-PI_EGU/", as.character(Sys.Date()),"_Effect_vs_Mud_Scatter.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(effect, aes(x = `% Fine Sand`, y = `Effect Size`))+
  geom_point()+
  geom_smooth(method = lm)

dev.off()

## Mud vs. Effect Size ####

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/ESS-PI_EGU/", as.character(Sys.Date()),"_Effect_vs_Mud.png"), width = 8, height = 8, units = "in", res = 300)

par(mar = c(5 ,6 , 4, 1))

plot(effect$`% Mud`, effect$`Effect Size`, cex= 1.8, cex.lab = 1.8, cex.axis = 1.8 ,xlab = "% Mud" , ylab = expression(paste("Effect Size (Wet - Dry Rate) (mg O"[2]*" L"^-1*" min"^-1*")")), lwd = 2)

#lines(lowess(effect$`% Mud`, effect$`Effect Size`), col = "blue", lwd = 3)

dev.off()

## D50 vs. Effect Size ####

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/ESS-PI_EGU/", as.character(Sys.Date()),"_Effect_vs_D50.png"), width = 8, height = 8, units = "in", res = 300)

par(mar = c(5 ,6 , 4, 1))

plot(effect$D50, effect$`Effect Size`,  cex= 1.8, cex.lab = 1.8, cex.axis = 1.8 , xlab = "D50" , ylab = expression(paste("Effect Size (Wet - Dry Rate) (mg O"[2]*" L"^-1*" min"^-1*")")), lwd = 2)

lines(lowess(effect$D50, effect$`Effect Size`), col = "blue", lwd = 3)

dev.off()

# James Correlation Matrix ####

mean_samples_corr_js <- mean_wet_dry_clean_corr %>%

  rename("Iron (mg/kg)" = "Mean_Fe_Treat") %>% 
  rename("Rate (mg/L)"  = "Average_Rate") %>% 
  rename("Effect Size" = "Effect_Size") %>% 
  rename("Sp. Conductivity" = "SpC_mean_chem") %>% 
  rename("pH" = "pH_mean_chem") %>% 
  rename("Initial Grav. Water" = "mean_grav_initial") %>% 
  rename("Final Grav. Water" = "mean_grav_final") %>% 
  rename("Lost Grav. Water" = "mean_grav_lost") %>% rename("Fine Sand (%)" = "Percent_Fine_Sand") %>% 
  rename("Medium Sand (%)" = "Percent_Med_Sand") %>% 
  rename("Coarse Sand (%)" = "Percent_Coarse_Sand") %>% 
  rename("Silt + Clay (%)" = 
           "percent_mud")

fe_do_chem_grn_wet <- fe_do_chem_grn_wet %>%
  rename("Initial Grav. Water" = "Initial Grav. Water Content") %>% 
  rename("Final Grav. Water" = "Final Grav. Water Content") %>% 
  rename("Lost Grav. Water" = "Lost Grav. Water Content")

panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r = (cor(x, y))
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

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/", as.character(Sys.Date()),"_James_Average_Correlation_Matrix.png"), width = 12, height = 12, units = "in", res = 300)

pairs(mean_wet_dry_clean_corr, lower.panel = panel.smooth,upper.panel = panel.cor, gap = 0, cex.labels = 1, cex = 1)

dev.off()

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/ESS-PI_EGU/", as.character(Sys.Date()),"_Dry_Correlation_Matrix.png"), width = 12, height = 12, units = "in", res = 300)

pairs(fe_do_chem_grn_dry, lower.panel = panel.smooth, upper.panel = panel.cor, gap = 0, cex.labels = 1, cex = 1)

dev.off()





#### PCA/RDA ####
## PCA ####
mean_wet_pca <- mean_wet %>% 
  dplyr::select(-c("Mean Initial Gravimetric Water", "Mean Lost Gravimetric Water"#, "log_effect",
                   #"Log_Mean_Treat_Fe_mg_kg",
                   #"% Coarse Sand", "% Med. Sand", "% Fine Sand", "% Mud", "geom", "geom_rusle", "Mean_Treat_Fe_mg_L"
                   ))

mean_wet_norm <- scale(mean_wet_pca)

mean_wet_norm_corr <- cor(mean_wet_norm)

corrplot(mean_wet_norm_corr, type = "upper", tl.col = "black", tl.cex = 1.6, cl.cex = 1.25)

mean_wet_pca <- princomp(mean_wet_norm_corr)
summary(mean_wet_pca)

mean_wet_pca$loadings[, 1:2]

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/ESS-PI_EGU/", as.character(Sys.Date()),"_Wet_Mean_PCA.png"), width = 11.5, height = 11.5, units = "in", res = 500)

fviz_pca_var(mean_wet_pca, col.var = "black")

dev.off()

mean_dry_pca <- mean_dry %>% 
  dplyr::select(-c("log_effect", "Log_Mean_Treat_Fe_mg_kg", "% Coarse Sand", "% Med. Sand", "% Fine Sand", "% Mud", "geom", "geom_rusle", "Mean_Treat_Fe_mg_L"))

mean_dry_norm <- scale(mean_dry_pca)

mean_dry_norm_corr <- cor(mean_dry_norm)
corrplot(mean_dry_norm_corr, type = "upper", tl.col = "black", tl.cex = 1, cl.cex = 1)

mean_dry_pca <- princomp(mean_dry_norm_corr)
summary(mean_dry_pca)

mean_dry_pca$loadings[, 1:2]

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/ESS-PI_EGU/", as.character(Sys.Date()),"_Dry_Mean_PCA.png"), width = 13, height = 13, units = "in", res = 500)

fviz_pca_var(mean_dry_pca, col.var = "black")

dev.off()

## RDA ####

mean_wet_rate <- mean_wet %>% 
  dplyr::select(c(`Average Rate (mg/L)`))

mean_wet_rate_hel <- decostand(mean_wet_rate, method = "hellinger")

mean_wet_pred <- mean_wet %>% 
  dplyr::select(-c(`Average Rate (mg/L)`))

mean_wet_pred_std <- decostand(mean_wet_pred, method = "standardize")

mean_wet_rda <- rda(mean_wet_rate ~., data = mean_wet_pred_std)

summary(mean_wet_rda)

ordiplot(mean_wet_rda, scaling = 1, type = "text")

fwd_wet <- ordiR2step(rda(mean_wet_rate$`Average Rate (mg/L)` ~ 1, data = mean_wet_pred_std), 
          scope = formula(mean_wet_rda), 
          direction = "forward", 
          R2scope = TRUE,
          pstep = 1000,
          trace = FALSE)

fwd_wet$call

mean_wet_signif <- rda(formula = mean_wet_rate$`Average Rate (mg/L)` ~ D50 + `Mean Initial Gravimetric Water`, data = mean_wet_pred_std)

RsquareAdj(mean_wet_signif)

anova.cca(mean_wet_signif, step = 1000)
anova.cca(mean_wet_signif, step = 1000, by = "term")


ordiplot(mean_wet_signif, scaling = 1, type = "text")

ordiplot(mean_wet_rda, scaling = 1, type = "text")

