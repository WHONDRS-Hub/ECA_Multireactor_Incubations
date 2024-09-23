library(tidyverse)

grain = read.csv("C:/GitHub/ECA_Multireactor_Incubations/Data/v3_CM_SSS_Sediment_Sample_Data_Summary.csv", skip = 2) %>% 
  filter(grepl("CM", Sample_Name)) %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "CM", "EC")) %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "Sediment", "all")) %>% 
  select(c(Sample_Name, Percent_Tot_Sand, Percent_Coarse_Sand, Percent_Med_Sand, Percent_Fine_Sand, Percent_Silt, Percent_Clay, Mean_Specific_Surface_Area_m2_per_g))


grain_range = grain %>%
  select(-c(Percent_Tot_Sand, Mean_Specific_Surface_Area_m2_per_g)) %>%
  pivot_longer(!Sample_Name, names_to = "Fraction", values_to = "Percent") %>%
  mutate(Percent = as.numeric(Percent)) %>% 
  mutate(Size = ifelse(Fraction == "Percent_Coarse_Sand", 0.500, ifelse(Fraction == "Percent_Med_Sand", 0.250, ifelse(Fraction == "Percent_Fine_Sand", 0.053, ifelse(Fraction == "Percent_Silt", 0.002, 0.000))))) %>% 
  group_by(Sample_Name) %>%
  mutate(P_Finer = case_when(
    Fraction == "Percent_Coarse_Sand" ~ Percent_Coarse_Sand + Percent_Med_Sand + Percent_Fine_Sand + Percent_Silt + Percent_Clay,
    Fraction == "Percent_Med_Sand" ~ Percent_Med_Sand + Percent_Fine_Sand + Percent_Silt + Percent_Clay, 
    Fraction == "Percent_Fine_Sand" ~ Percent_Fine_Sand + Percent_Silt + Percent_Clay,
    Fraction == "Percent_Silt" ~ Percent_Silt + Percent_Clay,
    Fraction == "Percent_Clay" ~ Percent_Clay))


 