#read in libraries

library(lubridate);library(writexl);library(raster);library(tidyverse);library(devtools);library(readxl)

##### Load data ######
rm(list=ls());graphics.off()

# Set working directory to data file
#Example:
pnnl.user = 'laan208'


#Read in all data
setwd(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/"))

#effect size - change date to most recent
effect.size <- paste0("Optode multi reactor/Optode_multi_reactor_incubation/rates/Effect_Size_merged_by_laan208_on_2023-04-10.csv")

effect <- read_csv(effect.size)

#Respiration rates
respiration <- paste0("Optode multi reactor/Optode_multi_reactor_incubation/rates/Plots/ECA_Sediment_Incubations_Respiration_Rates_merged_by_laan208_on_2023-03-07_pts_rem_res.csv")

resp <- read_csv(respiration)

#ECA Iron
iron <- paste0("FE/03_ProcessedData/20230110_Data_Processed_SFE_ECA_EC_1-270/20230110_Data_Processed_SFE_ECA_EC_1-270.csv")

fe <- read_csv(iron)

#ICON Grain Size
grain <- paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ICON_ModEx_SSS/09_Grain_Size/03_ProcessedData/20221230_Grain_Size_SBR_RC4_CM_1-42/20221230_Data_Processed_Grain_Size_SBR_RC4_CM_1-42.csv")

grn <- read_csv(grain)

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

all.chem = map %>% 
  dplyr::select(-`Disk Calibration date`, -Location, -`Time on`, -`Time off`, -Notes, -`Disk_ID`, -`map.file[i]`) %>% 
  filter(!grepl("Blank", Sample_ID))


###IRON DATA

#calculate mean Fe for kit/treatment
mean_fe <- fe %>% 
  separate(Sample_Name, into = c("Sample_ID", "rep"), sep = -1, convert = TRUE) %>% 
  group_by(Sample_ID) %>% 
  summarise_at(vars(Fe_mg_per_kg_sediment), list(mean = mean)) 

mean_fe <- mean_fe %>% 
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

##All Fe data
ggplot(mean_fe, aes(x = mean))+
  geom_histogram(binwidth = 0.1, fill = "cornflowerblue", col = "black")+
  theme_bw()+
  theme(axis.title.x = element_text(size = 18), 
        axis.title.y = element_text(size = 18), 
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16))+
  xlab("\n Fe (II) (mg/L)")+
  ylab("Count\n")

#All Fe data, faceted by wet vs dry
ggplot(mean_fe, aes(x = mean, fill = Treat))+
  geom_histogram(binwidth = 0.1)+
  facet_wrap(~Treat)+
  theme_bw()+
  theme(axis.title.x = element_text(size = 18), 
        axis.title.y = element_text(size = 18), 
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16))+
  xlab("\n Fe (II) (mg/L)")+
  ylab("Count\n")

#log Fe data, not faceted
mean_fe$log_iron_mean = log10(mean_fe$mean)

#png(file = paste0("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/ECA/EC 2022 Experiment/Optode multi reactor/Optode_multi_reactor_incubation/effect size/", as.character(Sys.Date()),"_log_fe_hist.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(mean_fe, aes(x = log_iron_mean))+
  geom_histogram(binwidth = 0.1, fill = "cornflowerblue", col = "black")+
  theme_bw()+
  theme(axis.title.x = element_text(size = 18), 
        axis.title.y = element_text(size = 18), 
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16))+
  xlab("\n Log of Fe (II) (mg/L)")+
  ylab("Count\n")

#dev.off()

#Log Fe, facetted by wet vs. dry
ggplot(mean_fe, aes(x = log_iron_mean, fill = Treat))+
  geom_histogram(binwidth = 0.1)+
  facet_wrap(~Treat)+
  theme_bw()+
  theme(axis.title.x = element_text(size = 18), 
        axis.title.y = element_text(size = 18), 
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16))+
  xlab("\n Log of Fe (II) (mg/L)")+
  ylab("Count\n")

#all means
mean_fe_all <- mean_fe %>% 
  group_by(kit) %>% 
  summarise_at(vars(mean), list(mean_all = mean))

#wet kit means
fe_wet <- mean_fe %>% 
  group_by(kit, Treat) %>% 
  summarise_at(vars(mean), list(mean_wet = mean)) %>% 
  filter(Treat != "Dry")

#dry kit means
fe_dry <- mean_fe %>% 
  group_by(kit, Treat) %>% 
  summarise_at(vars(mean), list(mean_dry = mean)) %>% 
  filter(Treat != "Wet")

#change names to merge with rates
fe.new <- mean_fe %>% 
  unite(kit_treat, kit, Treat, remove = FALSE) %>%
  mutate(rep = str_replace(rep, "SFE", "INC"))


###Grain Size Data

#this is the only data type that is not vial specific

grn <- grn %>% 
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
  filter(kit != "062")

grn$kit <- sub('.', '', grn$kit)

grn$percent_mud <- grn$Percent_Clay + grn$Percent_Silt

grn_long <- grn %>% 
  dplyr::select(-c("CM", "kit", "an")) %>% 
  pivot_longer(!Sample_ID, names_to = "size", values_to = "percent") %>% 
  filter(size != "Percent_Tot_Sand")

png(file = paste0("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/", as.character(Sys.Date()),"_GRN_stacked.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(grn_long, aes(fill = size, y = percent, x = Sample_ID)) + 
  geom_bar(position="fill", stat="identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust =1))

dev.off()


#pH, SpC, temp data
all.chem <- all.chem %>% 
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
  ))

#means, if needed
mean.chem <- all.chem %>% 
  group_by(kit, Treat) %>% 
  summarise_at(vars("SpC", "Temp", "pH"), mean) 

mean.chem <- mean.chem %>% 
  separate(Sample_ID, c("ECA", "kit", "Analysis"), sep = "_", remove = FALSE) 


###Iron vs. rates (FIX)
fe_do <- merge(fe.new, slope.new, by = c("rep", "kit_treat"), all = TRUE)

mean_fe$log_iron_mean = log10(mean_fe$mean)

lmFe = lm(mean~rate_mg_per_L_per_min, data = fe_do)

summary(lmFe)

png(file = paste0("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/ECA/EC 2022 Experiment/Optode multi reactor/Optode_multi_reactor_incubation/effect size/", as.character(Sys.Date()),"do_by_fe.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(fe_do, aes(x = mean, y = rate_mg_per_L_per_min))+
  geom_point()+
  stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth")+
  ylab(expression(paste("Rate (mg " ~L^-1 ~ min^-1~")"))) +
  xlab(expression(paste('Fe (II)  (mg' ~ L^-1~")")))+
  theme_bw()

dev.off()


###Grain Size vs. Iron
fe_grn <- merge (mean_fe_all, grn, by = "kit")

lmFeGrn = lm(mean_all~percent_mud, data = fe_grn)

summary(lmFeGrn)

ggplot(fe_grn, aes(x = percent_mud, y = mean_all))+
  geom_point()+
  stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth")+
  ylab(expression(paste('Fe (II)  (mg' ~ kg^-1~")"))) +
  xlab(expression(paste('Percent Mud')))+
  theme_bw()

fe_wet_grn <- merge(fe_wet, grn, by = "kit")

lmFeWetGrn = lm(mean_wet~percent_mud, data = fe_wet_grn)

summary(lmFeWetGrn)

ggplot(fe_wet_grn, aes(x = percent_mud, y = mean_wet))+
  geom_point()+
  stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth")+
  ylab(expression(paste('Fe (II) Wet (mg' ~ kg^-1~")"))) +
  xlab(expression(paste('Percent Mud')))+
  theme_bw()

fe_dry_grn <- merge(fe_dry, grn, by = "kit")

lmFeDryGrn = lm(mean_dry~percent_mud, data = fe_dry_grn)

summary(lmFeDryGrn)

ggplot(fe_dry_grn, aes(x = percent_mud, y = mean_dry))+
  geom_point()+
  stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth")+
  ylab(expression(paste('Fe (II) Dry (mg' ~ kg^-1~")"))) +
  xlab(expression(paste('Percent Mud')))+
  theme_bw()


### Fe vs. Grain Size vs. rates (FIX)

fe_do_grn <- merge(fe_do, grn, by = c("kit"), all = TRUE)

fe_do_grn <- fe_do_grn %>%
  filter(!is.na(slope_of_the_regression))

fe_do_grn_pca <- prcomp (fe_do_grn[,c(7, 10, 29:34)], center = TRUE, scale. = TRUE)

summary(fe_do_grn_pca)

ggbiplot(fe_do_grn_pca, labels = fe_do_grn$Sample_ID.x)

