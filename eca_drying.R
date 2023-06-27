library(readxl);library(corrplot);library(tidyverse)

rm(list=ls());graphics.off()

pnnl.user = 'laan208'

dry_wt <- read.csv(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/MOI/03_ProcessedData/EC_Moisture_Content_2022.csv"))

moisture <- paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/INC/01_RawData/2022_Data_Raw_INC_ECA_EC.xlsx")

all_moisture <- read_xlsx(moisture)

all_moisture$Date <- as.Date(all_moisture$Date)

all_moisture$Sample_weight_Fill_g <- as.numeric(all_moisture$Sample_weight_Fill_g)


#37 was incubated on 11/23 so needs to be removed from sheet, 72-5 were not incubated (not enough sediment), 12-D5 didn't have 20 g of sediment, 21 and 33 incubated on 9/28

all_moisture <- all_moisture %>% 
  filter(`Jars or 40 mL vials` != "Jars") %>% 
  filter(Sample_Name != "EC_37" & Date != "2022-11-23") %>% 
  filter(Sample_Name != "EC_72_INC-D5") %>% 
  filter(Sample_Name != "EC_72_INC-W5") %>% 
  filter(Sample_Name != "EC_12_INC-D5") %>% 
  filter(Sample_Name != "EC_21" & Date != "2022-09-28") %>% 
  filter(Sample_Name != "EC_33" & Date != "2022-09-28")

dry_wt$wet_g_by_dry_g <- dry_wt$wet_weight_g/dry_wt$true_dry_weight_g

mean_dry_wt <- dry_wt %>% 
  separate(sample_name, sep = "-", c("sample_name", "replicate")) %>% 
  group_by(sample_name) %>% 
  summarise_at(vars(wet_g_by_dry_g), funs(mean)) %>% 
  separate(sample_name, sep = "_", c("Study Code", "Site"))

all_moisture <- all_moisture %>% 
  separate(Sample_Name, into = c("EC", "Site", "INC"), sep = "_", remove = FALSE)

merged <- merge(all_moisture, mean_dry_wt, by = "Site")

merged_clean <- merged %>% 
  dplyr::select(-c(Site, EC, INC, Notes, `Jars or 40 mL vials`, ...8, ...9, ...10, `Study Code`))

##All Samples - wet and dry (all dates) ####

location = c("-W1", "-W2", "-W3", "-W4", "-W5", "-D1", "-D2", "-D3", "-D4", "-D5")

wet_dry <- as.data.frame(matrix(NA, ncol = 8, nrow =1))

colnames(wet_dry) = c("Sample_Name","Date",  "Tare_weight_g", "Sample_weight_g","Sample_weight_Fill_g", "dry_wt_sed_g","Water_added_initial_g", "Water_total_g")

all_dates = as.data.frame(matrix(NA, ncol = 8, nrow = length(unique(all_moisture$Sample_Name))))

colnames(all_dates) = c("Sample_Name","Date",  "Tare_weight_g", "Sample_weight_g","Sample_weight_Fill_g", "dry_wt_sed_g","Water_added_initial_g", "Water_total_g")

for (i in 1:length(location)){
  
  all_dates= as.data.frame(matrix(NA, ncol =8, nrow = length(unique(merged_clean$Sample_Name))))
  
  colnames(all_dates) = c("Sample_Name","Date",  "Tare_weight_g", "Sample_weight_g","Sample_weight_Fill_g", "dry_wt_sed_g","Water_added_initial_g", "Water_total_g")
  

  data_location_subset = merged_clean[grep(location[i],merged_clean$Sample_Name),]
  
  unique.incubations = unique(data_location_subset$Sample_Name)
  
  
  for (j in 1:length(unique.incubations)){
    
    data_site_subset = subset(data_location_subset, data_location_subset$Sample_Name == unique.incubations[j])
    
    data_site_subset <- data_site_subset[with(data_site_subset, order(Date)),]
    
    data_site_subset = data_site_subset %>%
      mutate(Tare_weight_g = first(Tare_weight_g))
    
    data_site_subset = data_site_subset %>% 
      mutate(dry_wt_sed_g = (first(Sample_weight_g) - Tare_weight_g)*(1/wet_g_by_dry_g)) 
  
  
    merge_dates= as.data.frame(matrix(NA, ncol = 8, nrow = length(unique(all_moisture$Sample_Name))))
    
    colnames(merge_dates) = c("Sample_Name","Date",  "Tare_weight_g", "Sample_weight_g", "Sample_weight_Fill_g", "dry_wt_sed_g", "Water_added_initial_g", "Water_total_g")
    
  for (k in 1:nrow(data_site_subset)){
    
    if (data_site_subset$Sample_weight_Fill_g[k] > -9999) {
     
      merge_dates$Sample_Name[k] = as.character(data_site_subset$Sample_Name[k])
    
    merge_dates$Date[k] = as.character(data_site_subset$Date[k]) 
    
    merge_dates$Tare_weight_g[k] = as.numeric(data_site_subset$Tare_weight_g[k]) 
    
    merge_dates$Sample_weight_g[k] = as.numeric(data_site_subset$Sample_weight_g[k] - data_site_subset$Tare_weight_g[k])
    
    merge_dates$Sample_weight_Fill_g[k] = as.numeric(data_site_subset$Sample_weight_Fill_g[k] - data_site_subset$Tare_weight_g[k])
    
    merge_dates$dry_wt_sed_g[k] = as.numeric(data_site_subset$dry_wt_sed_g[k])
    
    merge_dates$Water_added_initial_g[k] = as.numeric(data_site_subset$Sample_weight_Fill_g[1] - data_site_subset$Sample_weight_g[1])
    
  merge_dates$Water_total_g[k] = merge_dates$Sample_weight_Fill_g[k] - merge_dates$dry_wt_sed_g[k]
    
    }
    
    else {
      
      merge_dates$Sample_Name[k] = as.character(data_site_subset$Sample_Name[k])
      
      merge_dates$Date[k] = as.character(data_site_subset$Date[k]) 
      
      merge_dates$Tare_weight_g[k] = as.numeric(data_site_subset$Tare_weight_g[k]) 
      
      merge_dates$Sample_weight_g[k] = as.numeric(data_site_subset$Sample_weight_g[k] - data_site_subset$Tare_weight_g[k])
      
      merge_dates$Water_added_initial_g[k] = as.numeric(0)
      
      merge_dates$dry_wt_sed_g[k] = as.numeric(data_site_subset$dry_wt_sed_g[k])
      
      merge_dates$Water_total_g[k] = merge_dates$Sample_weight_g[k] - merge_dates$dry_wt_sed_g[k]
      
    }
  
    }
    
    #all_dates = all_dates %>%
     # mutate(Sample_weight_initial_g = first(Sample_weight_initial_g), round, 4)
    
   all_dates = rbind(all_dates, merge_dates)
    
  }
  
  
  
  all_dates = all_dates[!is.na(all_dates$Sample_Name),]
  
  wet_dry = rbind(wet_dry, all_dates)
  
}

wet_dry = wet_dry[-1,] 

wet_dry <- wet_dry %>% 
  dplyr::select(-c(Tare_weight_g))

write.csv(wet_dry,paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/INC/03_ProcessedData/ECA_Drying_Masses_merged_by_",pnnl.user,"_on_",Sys.Date(),".csv"), row.names = F)






wet_wt <- paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/MOI/03_ProcessedData/EC_Moisture_Content_2022.csv")

corr <- read_csv(wet_wt)

mean_wet_wt <- corr %>% 
  separate(col = sample_name, into = c("Project", "kit", "analysis"), sep = "_") %>%  
  separate(col = analysis, into = c("Analysis", "Replicate"), sep = "-") %>% 
  group_by(kit) %>% 
  summarise(mean_wet_grav = (mean(percent_water_content_wet)/100))

pivot = wet_dry %>% 
  pivot_longer(cols = ends_with("_g"),
               names_to = "Sample_weight_type",
               values_to = "Sample_weight_g")

grav <- wet_dry %>% 
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
  unite(kit_Treat, kit, Treat) %>% 
  select(-c(Project, Analysis, Date, Tare_weight_g, Sample_weight_g)) %>% 
  distinct(kit_Treat, Replicate, .keep_all = TRUE)

average_grav <- all_grav_final %>% 
  group_by(kit_Treat) %>% 
  mutate(average_grav_intial = mean(grav_dry_initial)) %>% 
  mutate(average_grav_final = mean(grav_dry_final)) %>% 
  mutate(average_grav_lost_subt = average_grav_intial - average_grav_final) %>% 
  mutate(average_grav_lost = mean(lost_grav_perc)) %>% 
  distinct(kit_Treat, .keep_all = TRUE) %>% 
  dplyr::select(c(kit_Treat,average_grav_intial,average_grav_final,average_grav_lost))

average_grav_lost <- average_grav %>% 
  separate(col = kit_Treat, into = c("kit", "Treat"), sep = "_") %>% 
  group_by(kit) %>% 
mutate(grav_final_diff = (average_grav_final[Treat == "Wet"] - average_grav_final[Treat == "Dry"])) %>% 
  dplyr::select(c(kit, grav_final_diff)) %>% 
  distinct(kit, .keep_all = TRUE)

ess_ex <- all_grav %>% 
  filter(kit == "05") %>% 
  filter(Replicate == "W1"| Replicate == "D1") %>% 
  mutate(mass_sed = (Sample_weight_initial_g - (Sample_weight_initial_g*mean_wet_grav)))%>% 
  mutate(grav_dry = (((Sample_weight_g - mass_sed) + Water_added_initial_g)/mass_sed))

ess_ex$Date <- as.Date(ess_ex$Date, format = "%Y-%m-%d")

cbPalette <- c("#D55E00","#0072B2","#999999","#E69F00", "#56B4E9","#009E73","#F0E442","#CC79A7")

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/ESS-PI_EGU/", as.character(Sys.Date()),"05_R1_Grav_Water.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(ess_ex, aes(x = Date, y = grav_dry))+
  geom_point(aes(color = Replicate), size = 2.5)+
  geom_line(aes(color = Replicate), linewidth = 1.5)+
  theme_bw()+
  ylab("Gravimetric water content \n")+
  xlab("\n Date")+
  scale_x_date(date_labels = "%Y - %m - %d", #breaks = unique(ess_ex$Date)
               )+
  theme(axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text.x = element_text(size = 16, angle = 90, hjust = 1, vjust = 1),
        axis.text.y = element_text(size =18),
        title = element_text(size = 24),
        legend.text = element_text(size = 18))+
  scale_colour_manual(values=cbPalette, labels = c("Dry", "Wet"))

dev.off()