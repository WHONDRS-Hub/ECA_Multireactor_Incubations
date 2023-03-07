library(tidyverse);library(dplyr);library(readxl)

rm(list=ls());graphics.off()

pnnl.user = 'laan208'

input.path = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/FE/")

setwd(input.path)

raw.data = ("01_Rawdata/")
formatted.data = ("02_FormattedData/")
processed.data = ("03_ProcessedData/")


#read in masses of sediment used in incubation

inc.masses <- paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/INC/01_RawData/")

#import_iron = function(input.path){
  
  # import map
  ferrozine_map = read_excel("01_Rawdata/20230110_Data_Raw_SFE_ECA_EC_1-270/20230110_Mapping_SFE_ECA_EC_1-270.xlsx", sheet = "map") %>% mutate_all(as.character) 
  # import data files (plate reader)
  filePaths_ferrozine <- list.files(path = raw.data, pattern = "xlsx", full.names = TRUE, recursive = TRUE)
  ferrozine_data <- do.call(bind_rows, lapply(filePaths_ferrozine, function(raw.data) {
    df <- read_xlsx(raw.data, skip = 24) %>% mutate_all(as.character) %>% janitor::clean_names()
    df = df %>% mutate(source = basename(raw.data))
    df}))
  
 # list(ferrozine_map = ferrozine_map,
 #     ferrozine_data = ferrozine_data)
  #}



# clean the map
map_processed = 
  ferrozine_map %>% 
  mutate(tray_number = parse_number(tray_number)) %>% 
  filter(!is.na(sample_name) & !is.na(tray_number)) %>% 
  rename(sample_label = sample_name)

#remove trays 4 and 6, figure out how to remove bad dilutions in 12 and 13

data_formatted = 
  ferrozine_data %>% 
  mutate_all(na_if,"") %>% 
  rename("x" = "x1") %>% 
  rename("x1" = "x1_2") %>% 
  dplyr::select(-x14) %>% 
  #fill(x_1) %>% 
  #filter(x_2 == "562") %>% 
  #dplyr::select(-x_2) %>% 
  pivot_longer(-c(source, x), values_to = "absorbance_562") %>% 
  mutate(name = str_remove(name, "x"),
         well_position = paste0(x, name),
         tray_number = str_extract(source, "Tray[0-9]+"),
         tray_number = parse_number(tray_number),
         absorbance_562 = as.numeric(absorbance_562)) %>% 
  dplyr::select(tray_number, well_position, absorbance_562) %>% 
  right_join(map_processed, by = c("tray_number", "well_position")) %>% 
  filter(!notes %in% "skip") %>%
  filter(tray_number != "4") %>% 
  filter(tray_number != "6") 
 

data_flag <- data_formatted %>%
  separate(col = sample_label, into = c("Project", "kit", "analysis"), sep = "_") %>% 
  separate(col = analysis, into = c("Analysis", "Replicate"), sep = "-") %>% 
  separate(Replicate, into = c("Replicate", "Technical"), sep = "(?<=\\d)(?=[a-z]?)") %>% 
  group_by(kit, Replicate) %>% 
  summarise(CV = ((sd(absorbance_562)/mean(absorbance_562))*100))
  
#### Formatted absorbance data ####
write.csv(data_formatted,"02_FormattedData/20230110_Data_Formatted_SFE_ECA_EC_1-270/20230110_Data_Formatted_SFE_ECA_EC_1-270.csv", row.names = F)
  

####Processed Data #####

#choose which standard curve to use for Fe 

calibrate_ferrozine_data = function(data_formatted){
  standards = 
    data_formatted %>% 
    filter(grepl("standard", sample_label)) %>% 
    dplyr::select(tray_number, absorbance_562, standard_ppm) %>% 
    mutate(standard_ppm = as.numeric(standard_ppm))
  
  standards %>% 
    ggplot(aes(x = standard_ppm, y = absorbance_562, color = as.character(tray_number)))+
    geom_point()+
    geom_smooth(method = "lm", se = F)+
    facet_wrap(~tray_number)
  
  calibration_coef = 
    standards %>% 
    dplyr::group_by(tray_number) %>% 
    dplyr::summarize(slope = lm(absorbance_562 ~ standard_ppm)$coefficients["standard_ppm"], 
                     intercept = lm(absorbance_562 ~ standard_ppm)$coefficients["(Intercept)"])
  
  # y = mx + c
  # abs = m*ppm + c
  # ppm = abs-c/m
  
  
 data_formatted = data_formatted %>% 
    left_join(calibration_coef) %>% 
    mutate(ppm_calculated = ((absorbance_562 - intercept) / slope))
  
}

samples = 
  calibrate_ferrozine_data(data_formatted) %>% 
  filter(grepl("EC", sample_label)) %>% 
  dplyr::select(sample_label, analysis, ppm_calculated) %>% 
  mutate(ppm_calculated = case_when(
  analysis  == "dilute" ~ ppm_calculated * 2, 
  analysis == "Fe2" ~ ppm_calculated)) %>% 
  rename("mg_Fe_per_L" = "ppm_calculated") %>% 
  mutate(mg_Fe_per_L = if_else(mg_Fe_per_L<0,0,mg_Fe_per_L))

#ggplot(samples) +
 # geom_boxplot(aes(x = sample_label, y = mg_Fe_per_L))

#largest sd 0.18, for diluted sample with differences in abs 0.2

blanks = 
  calibrate_ferrozine_data(data_formatted) %>% 
  filter(grepl("blank", sample_label)) 
  
# mean(blanks$ppm_calculated)
# sd(blanks$ppm_calculated)
# 
# mean(blanks$absorbance_562)
# sd(blanks$absorbance_562)
  
### pull in moisture data to correct to mg Fe/kg dry sediment

moisture <- read.csv(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/MOI/03_ProcessedData/EC_Moisture_Content_2022.csv"))

moisture$wet_g_by_dry_g <- moisture$wet_weight_g/moisture$true_dry_weight_g

mean.moisture <- moisture %>% 
  separate(sample_name, sep = "-", c("sample_name", "replicate")) %>% 
  group_by(sample_name) %>% 
  summarise_at(vars(wet_g_by_dry_g), funs(mean)) %>% 
  separate(sample_name, sep = "_", c("Study Code", "Site"))

#merge moisture and Fe samples

samples_2 <- samples %>% 
  separate(sample_label, sep = "_", c("Study Code", "Site", "Replicate")) 

merged <- merge(samples_2,mean.moisture, by = "Site")

merged_sep <- merged %>% 
  separate(Replicate, sep = "-", c("Iron",  "Treat"))
  
  merged_sep$`Technical Replicate` <- str_sub(merged_sep$Treat, 3, - 1) 
  
  merged_sep$Treat <- str_sub(merged_sep$Treat, 0,2)
###mapping files 

import_data = function(inc.masses){ 
  
  #pull a list of files in target folder with correct pattern
  #read all files and combine
  
  map.file <-  list.files(inc.masses, recursive = T, pattern = "\\.xlsx$",full.names = T)
  
  map.file <- map.file[!grepl("red|Red|EC_01_|EC_02_|EC_03_|EC_06_07|EC_10_15|EC_04_08|fast", map.file)]
  
  mapping <- lapply(map.file, read_xlsx)
  
  for (i in 1:length(mapping)){mapping[[i]] <- cbind(mapping[[i]], map.file[i])}
  
  all.map <- 
    do.call(rbind,mapping)
}

all.map = import_data(inc.masses)

all.map = all.map %>% 
  dplyr::select(-`...8`, -`...9`, -`...10`, -`map.file[i]`) %>% 
  filter(Tare_weight_g != -9999) %>% 
  filter(`Jars or 40 mL vials` != "Jars") %>% separate(Sample_Name, sep = "_", c("Study Code", "Site", "Inc")) %>% 
  separate(Inc, sep = "-", c("Inc", "Treat"))
  
all.map$wet_sediment_wt_g <- all.map$Sample_weight_g - all.map$Tare_weight_g


###merging Fe, moisture, and wet sediment mass in vial

final.merge <- merge(merged_sep, all.map, by = c("Site", "Treat"))


###assuming 50 mL water in 50 mL vial

final.merge$mg_Fe_per_kg_sediment <- final.merge$mg_Fe_per_L*(0.05)*(1/final.merge$wet_sediment_wt_g)*(final.merge$wet_g_by_dry_g)*1000

processed.data <- final.merge %>% 
  dplyr::select(-`Study Code.y`, -`Study Code`, -Date, -Inc, -`Sample_weight_Fill_g`, -Notes,-`Jars or 40 mL vials`, -wet_g_by_dry_g, -Tare_weight_g, -Sample_weight_g, - wet_sediment_wt_g) %>% 
  unite(Sample_Name, `Study Code.x`, Site,sep = "_") %>% 
  unite(Sample_Name, Sample_Name, Iron, sep = "_") %>% 
  unite(Sample_Name, Sample_Name, Treat, sep = "-") %>% 
  unite(Sample_Name, Sample_Name, `Technical Replicate`, sep = "")


write.csv(processed.data, "C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/FE/03_ProcessedData/20230110_Data_Processed_SFE_ECA_EC_1-270/20230110_Data_Processed_SFE_ECA_EC_1-270.csv")
 


####not finished
means = 
  samples %>% 
  group_by(sample_label) %>% 
  summarise_at(vars(ppm_calculated),
               funs(mean,sd)) %>% 
  rename(mean, mean_ppm_calculated = mean) %>% 
  mutate(mean_ppm_calculated = if_else(mean_ppm_calculated<0,0,mean_ppm_calculated))

png(file = paste0("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/ECA/EC 2022 Experiment/Optode multi reactor/Optode_multi_reactor_incubation/effect size/", as.character(Sys.Date()),"_all_iron_hist.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(means, aes(x = mean_ppm_calculated)) +
geom_histogram(binwidth = 0.02, fill = "cornflowerblue", col = "black")+
  theme_bw()+
  theme(axis.title.x = element_text(size = 18), 
        axis.title.y = element_text(size = 18), 
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16))+
  xlab("\nFe (II) (mg/L)")+
  ylab("Count\n")
 
dev.off()

all_treat <- means %>% 
  separate(sample_label, sep = "-", c("kit", "rep")) %>% 
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
  
mean_treat = 
  all_treat %>% 
  group_by(kit) %>% 
  mutate(mean_kit = mean(mean_ppm_calculated)) %>%
  mutate(sd_kit = sd(mean_ppm_calculated)) %>% 
  ungroup() %>% 
  group_by(kit, Treat) %>% 
  mutate(mean_kit_treat = mean(mean_ppm_calculated)) %>% 
  mutate(sd_kit_treat = sd(mean_ppm_calculated)) %>% 
  ungroup() %>% 
  mutate(cv_kit = sd_kit/mean_kit) %>% 
  mutate(cv_kit_treat = sd_kit_treat/mean_kit_treat)
  

slope.new <- slope.new %>% 
  group_by(kit_treat) %>% 
  mutate(Slope.All = mean(slope_of_the_regression)) %>% 
  ungroup()

#figure out how to use this - right now can't because of triplicates
%>%
  pivot_wider(names_from = "analysis", values_from = "ppm_calculated") %>%
  mutate(across(where(is.numeric), round, 2))

samples2 = 
  samples %>% 
  dplyr::select(sample_label, starts_with("Fe")) %>% 
  pivot_longer(cols = starts_with("Fe"), names_to = "species", values_to = "ppm") %>% 
  left_join(moisture_processed) %>% 
  left_join(subsampling %>% dplyr::select(sample_label, iron_g)) %>% 
  rename(fm_g = iron_g) %>% 
  mutate(ppm = as.numeric(ppm),
         od_g = fm_g/((gwc_perc/100)+1),
         soilwater_g = fm_g - od_g,
         ug_g = ppm * ((25 + soilwater_g)/od_g),
         ug_g = round(ug_g, 2)) %>% 
  dplyr::select(sample_label, species, ppm, ug_g) %>% 
  arrange(sample_label) %>% 
  pivot_longer(-c(sample_label, species)) %>% 
  mutate(name = paste0(species, "_", name)) %>% 
  dplyr::select(-species) %>% 
  filter(!grepl("blank", sample_label)) %>% 
  pivot_wider()%>% 
  mutate(analysis = "Ferrozine")

samples2

