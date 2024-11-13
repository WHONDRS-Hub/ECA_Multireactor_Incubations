library(tidyverse)

cinc = read.csv("C:/GitHub/Cincinnati_Multireactor_Respiration/Data/RS_Sediment_Incubations_Respiration_Rates_ReadyForBoye_2024-08-28.csv") %>% 
  select(c(Sample_Name, Respiration_Rate_mg_DO_per_kg_per_H))

ev = read.csv("C:/Users/laan208/OneDrive - PNNL/Shared Documents - Core Richland and Sequim Lab-Field Team/Data Generation and Files/ECA/INC/03_ProcessedData/EV_Sediment_Incubations_Respiration_Rates_ReadyForBoye_2024-11-01.csv") %>% 
  filter(Respiration_Rate_mg_DO_per_kg_per_H != "-9999") %>% 
  mutate(Respiration_Rate_mg_DO_per_kg_per_H = as.numeric(Respiration_Rate_mg_DO_per_kg_per_H))%>% 
  select(c(Sample_Name, Respiration_Rate_mg_DO_per_kg_per_H, Methods_Deviation))

eca = read.csv("C:/GitHub/ECA_Multireactor_Incubations/Data/EC_Sediment_SpC_pH_Temp_Respiration.csv", skip = 2) %>% 
  filter(grepl("EC", Sample_Name)) %>% 
  filter(Respiration_Rate_mg_DO_per_kg_per_H != "-9999") %>% 
  mutate(Respiration_Rate_mg_DO_per_kg_per_H = as.numeric(Respiration_Rate_mg_DO_per_kg_per_H))%>% 
  select(c(Sample_Name, Respiration_Rate_mg_DO_per_kg_per_H, Methods_Deviation))

mel = read.csv("C:/Users/laan208/OneDrive - PNNL/Shared Documents - Core Richland and Sequim Lab-Field Team/Data Generation and Files/ECA/INC/03_ProcessedData/EL_Sediment_Incubations_Respiration_Rates_ReadyForBoye_2024-10-24.csv")%>% 
  filter(Respiration_Rate_mg_DO_per_kg_per_H != "-9999") %>% 
  mutate(Respiration_Rate_mg_DO_per_kg_per_H = as.numeric(Respiration_Rate_mg_DO_per_kg_per_H))%>% 
  select(c(Sample_Name, Respiration_Rate_mg_DO_per_kg_per_H, Methods_Deviation))


all_data = bind_rows(cinc, eca, ev, mel) %>% 
  mutate(Data = if_else(grepl("EC", Sample_Name), "CONUS Sediment", if_else(grepl("EV", Sample_Name), "Texas Sediment", if_else(grepl("EL", Sample_Name), "CONUS Soil", if_else(grepl("R0", Sample_Name), "Cincinnati Sediment", "Cincinnati Soil"))))) 

ggplot(all_data, aes(x = Respiration_Rate_mg_DO_per_kg_per_H, fill = Data)) + 
  geom_histogram()+ 
  theme_bw()+
  theme(axis.title.x = element_text(size = 15))

real_data = all_data %>% 
  filter(!grepl("RATE_005", Methods_Deviation))

real_data %>% 
  filter(Respiration_Rate_mg_DO_per_kg_per_H >-20) %>% 
ggplot(aes(x = Respiration_Rate_mg_DO_per_kg_per_H, fill = Data)) + 
  geom_histogram()+ 
  theme_bw()+
  theme(axis.title.x = element_text(size = 15))


                        