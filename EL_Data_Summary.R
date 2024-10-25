## ECA Data Package Summary File
# 5/29/2024 M.Laan

library(tidyverse); library(readxl)

rm(list=ls());graphics.off()

# Set working directory to data file
#Example:
pnnl.user = 'laan208'

#Set wd to Boye Files
setwd("Z:/00_Cross-SFA_ESSDIVE-Data-Package-Upload/01_Study-Data-Package-Folders/ECA_Data_Package/EC_Data_Package/Sample_Data/")


# Respiration -------------------------------------------------------------

all_respiration <- read.csv("C:/Users/laan208/OneDrive - PNNL/Shared Documents - Core Richland and Sequim Lab-Field Team/Data Generation and Files/ECA/INC/03_ProcessedData/EL_Sediment_Incubations_Respiration_Rates_ReadyForBoye_2024-10-24.csv") %>%  
  dplyr::select(c(Sample_Name, SpC, pH, Temp, Respiration_Rate_mg_DO_per_L_per_H, Respiration_Rate_mg_DO_per_kg_per_H, Methods_Deviation)) %>% 
  mutate(across(c(SpC:Respiration_Rate_mg_DO_per_kg_per_H), as.numeric))


median_respiration = all_respiration %>% 
  mutate(Respiration_Rate_mg_DO_per_L_per_H = ifelse(grepl("INC_Method_001|INC_Method_002|INC_QA_004", Methods_Deviation), NA, Respiration_Rate_mg_DO_per_L_per_H)) %>% 
  #missing replicates (EC_072-W5/D5),  overexposed samples (EC_027, EC_013, EC_014), less sediment in sample (EC_012-D5)
  mutate(Respiration_Rate_mg_DO_per_kg_per_H = ifelse(grepl("INC_Method_001|INC_Method_002|INC_QA_004", Methods_Deviation), NA, Respiration_Rate_mg_DO_per_kg_per_H)) %>% 
  mutate(SpC = ifelse(grepl("INC_Method_001|INC_Method_002", Methods_Deviation), NA, SpC)) %>% 
  # mutate(SpC_microsiemens_per_cm = ifelse(grepl("INC_Method_001|INC_Method_002", Methods_Deviation), NA, SpC_microsiemens_per_cm)) %>% 
  mutate(pH = ifelse(grepl("INC_Method_001|INC_Method_002", Methods_Deviation), NA, pH)) %>% 
  mutate(Temp = ifelse(grepl("INC_Method_001|INC_Method_002", Methods_Deviation), NA, Temp)) %>%
  # mutate(Temperature_degC = ifelse(grepl("INC_Method_001|INC_Method_002", Methods_Deviation), NA, Temperature_degC)) %>% 
  mutate(Respiration_Rate_mg_DO_per_kg_per_H = ifelse(Respiration_Rate_mg_DO_per_kg_per_H == "-9999", NA, Respiration_Rate_mg_DO_per_kg_per_H)) %>% 
  separate(Sample_Name, c("Sample_ID", "Rep"), sep = "-") %>% 
  mutate(Rep = if_else(grepl("D", Rep), "D", "W")) %>%
  group_by(Sample_ID, Rep) %>%
  summarise(across(where(is.numeric),
                   list(Median = ~median(.x, na.rm = TRUE),
                        cv = ~sd(.x, na.rm = TRUE)/mean(.x, na.rm =TRUE),
                        n = ~sum(!is.na(.x))), 
                   .names = "{.fn}_{.col}")) %>% 
  ungroup() %>% 
  group_by(Sample_ID, Rep) %>%
  mutate(Remove = ifelse(all(c_across(starts_with("n_")) == 5), "FALSE", "TRUE")) %>% ## Check CV's, then remove 
  dplyr::select(c(Sample_ID, Rep, Median_SpC, Median_pH, Median_Temp, Median_Respiration_Rate_mg_DO_per_L_per_H, Median_Respiration_Rate_mg_DO_per_kg_per_H, Remove)) %>% 
  # select(c(Sample_ID, Rep, Median_SpC_microsiemens_per_cm, Median_pH, Median_Temperature_degC, Median_Respiration_Rate_mg_DO_per_L_per_H, Median_Respiration_Rate_mg_DO_per_kg_per_H, Remove)) %>% 
  unite(Sample_Name, c("Sample_ID", "Rep")) %>% 
  mutate(Median_Respiration_Rate_mg_DO_per_L_per_H = Median_Respiration_Rate_mg_DO_per_L_per_H * -1) %>% 
  mutate(Median_Respiration_Rate_mg_DO_per_kg_per_H = Median_Respiration_Rate_mg_DO_per_kg_per_H * -1) %>% 
  separate(Sample_Name, c("EL", "Site", "IN", "Treat"), sep = "_") %>% 
  unite(Sample_Name, c("EL", "Site", "IN"), sep = "_", remove = F) %>% 
  ungroup()

## Effect Size ####

effect_data <- median_respiration %>% 
  group_by(Sample_Name) %>% 
  mutate(across(where(is.numeric), ~. [Treat == "W"] - .[Treat  == "D"])) %>% 
  rename_with(~ str_replace_all(., "Median_", "Effect_Size_")) %>% 
  distinct(Sample_Name, .keep_all = TRUE) %>% 
  select(-c(Treat)) %>% 
  ungroup()
## Aggregates ####

agg = read_xlsx("C:/GitHub/ECA_Multireactor_Incubations/Data/EL/water_stable_aggregates_10_24_2024.xlsx") %>% 
  mutate(Sample_Name = str_replace(sample_name, "AG", "IN")) %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "-1", "")) %>% 
  select(-c(sample_name)) %>%
  relocate(Sample_Name, .before=`53um_aggregates_wt%`) %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "MEL", "EL"))

## Figures ####

resp_agg = left_join(agg, median_respiration, by = "Sample_Name") %>% 
  select(-c(Remove, Median_SpC, Median_Temp, Median_pH)) %>% 
  relocate(Treat, .after = "Sample_Name") %>% 
  pivot_longer(cols = starts_with("53um_aggregates_wt%"):starts_with("total%_aggregates"), 
               names_to = "agg_size", values_to = "perc")

effect_agg = left_join(agg, effect_data, by = "Sample_Name") %>% 
  select(-c(Remove, Effect_Size_SpC, Effect_Size_Temp, Effect_Size_pH)) %>% 
  pivot_longer(cols = starts_with("53um_aggregates_wt%"):starts_with("total%_aggregates"), 
               names_to = "agg_size", values_to = "perc")

ggplot(resp_agg, aes(x = perc)) +
  geom_histogram()+ 
  facet_wrap(~agg_size)

ggplot(resp_agg, aes(x = Median_Respiration_Rate_mg_DO_per_kg_per_H)) + 
  geom_histogram()

ggplot(resp_agg, aes(x = Median_Respiration_Rate_mg_DO_per_L_per_H)) + 
  geom_histogram()

ggplot(resp_agg, aes(x = perc, y = Median_Respiration_Rate_mg_DO_per_L_per_H)) +
  geom_point(aes(color = Treat, shape = IN)) + 
  facet_wrap(~agg_size)

ggplot(resp_agg, aes(x = perc, y = Median_Respiration_Rate_mg_DO_per_kg_per_H)) +
  geom_point(aes(color = Treat, shape = IN)) + 
  facet_wrap(~agg_size)

ggplot(effect_agg, aes(x = perc, y = Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H)) +
  geom_point(aes(color = IN))+
  facet_wrap(~agg_size)

ggplot(effect_agg, aes(x = perc, y = Effect_Size_Respiration_Rate_mg_DO_per_L_per_H)) +
  geom_point(aes(color = IN))+
  facet_wrap(~agg_size)

## PCA? ####

effect_agg_wide = left_join(agg, effect_data, by = "Sample_Name")


effect_agg_pca = effect_agg_wide %>% 
  select(-c(Remove, Effect_Size_SpC, Effect_Size_Temp, Effect_Size_pH, Site, EL, IN, #Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H, 
            Effect_Size_Respiration_Rate_mg_DO_per_L_per_H)) %>% column_to_rownames("Sample_Name")

effect_agg_scale = scale(effect_agg_pca)

pca_result <- prcomp(effect_agg_scale, center = TRUE)

# View PCA summary
summary(pca_result)

# Get PCA loadings
pca_result$rotation

# Get scores (PC coordinates for each observation)
pca_scores = as.data.frame(pca_result$x)
pca_scores$Group = effect_agg_wide$IN


library(ggbiplot)

ggbiplot(pca_result, obs.scale = 1, var.scale = 1, groups = effect_agg_wide$IN, ellipse = T) + 
  #scale_color_gradient(low = "blue", high = "red", name = "Effect Size") + 
  theme_minimal()+
  theme(plot.title = element_text(size = 16, hjust = 0.5),         # Title text size
    axis.title = element_text(size = 14),                      # Axis title text size
    axis.text = element_text(size = 12),                       # Axis tick text size
    legend.title = element_text(size = 12),                    # Legend title text size
    legend.text = element_text(size = 10)                     # Legend item text size
  )

  