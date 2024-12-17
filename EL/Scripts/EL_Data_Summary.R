## ECA Data Package Summary File
# 5/29/2024 M.Laan

#pull out 12D-H5, 12L-W5, 16L-D2

library(tidyverse); library(readxl); library(ggpmisc); library(corrplot); library(Hmisc)

rm(list=ls());graphics.off()

# Set working directory to data file
#Example:
pnnl.user = 'laan208'



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
  unite(Site, c("EL", "Site"), sep = "_") %>% 
  ungroup()

ggplot(median_respiration, aes(x = Median_Respiration_Rate_mg_DO_per_kg_per_H)) + 
  geom_histogram() + 
  facet_grid(~Treat, scales = "free")

## Effect Size ####

effect_data <- median_respiration %>% 
  group_by(Sample_Name) %>% 
  mutate(across(where(is.numeric), ~. [Treat == "W"] - .[Treat  == "D"])) %>% 
  rename_with(~ str_replace_all(., "Median_", "Effect_Size_")) %>% 
  distinct(Sample_Name, .keep_all = TRUE) %>% 
  select(-c(Treat)) %>% 
  ungroup()

## Aggregates ####

agg = read_xlsx("C:/GitHub/ECA_Multireactor_Incubations/EL/Data/water_stable_aggregates_10_24_2024.xlsx") %>% 
  mutate(Sample_Name = str_replace(sample_name, "AG", "IN")) %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "-1", "")) %>% 
  select(-c(sample_name)) %>%
  relocate(Sample_Name, .before=`53um_aggregates_wt%`) %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "MEL", "EL")) %>% 
  separate(Sample_Name, c("Site", "Level"), sep = "_IN", remove = F)

## Figures ####

resp_agg = left_join(agg, median_respiration, by = c("Sample_Name", "Site")) %>% 
  select(-c(Remove, Median_SpC, Median_Temp, Median_pH)) %>% 
  relocate(Treat, .after = "Sample_Name") %>% 
  pivot_longer(cols = starts_with("53um_aggregates_wt%"):starts_with("total%_aggregates"), 
               names_to = "agg_size", values_to = "perc")

effect_agg = left_join(agg, effect_data, by = c("Sample_Name", "Site")) %>% 
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


## PCA? ####

effect_agg_wide = left_join(agg, effect_data, by = c("Sample_Name", "Site"))

effect_real = effect_agg_wide %>% 
  filter(Effect_Size_Respiration_Rate_mg_DO_per_L_per_H < 100)

ggplot(effect_real, aes(x = reorder(Sample_Name, Effect_Size_Respiration_Rate_mg_DO_per_L_per_H), y = Effect_Size_Respiration_Rate_mg_DO_per_L_per_H)) +
  geom_col()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 


effect_agg_pca = effect_agg_wide %>% 
  select(-c(Remove, Effect_Size_SpC, Effect_Size_Temp, Effect_Size_pH, Site, IN, Level, Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H, 
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


## Effect Size Histogram

## EFFECT CUBE ROOT


effect_limits <- c(-1600, 1600)


ggplot(effect_data, aes(x = Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H))+
  # geom_histogram(binwidth = 0.15, fill = "#009E73")+
  geom_histogram(binwidth = 20, aes(fill = after_stat(x))) +
  scale_fill_gradient2(name = "Effect Size", limits = effect_limits, low = "firebrick2", mid = "goldenrod2",
                       high = "dodgerblue2", midpoint = (max(effect_limits)+min(effect_limits))/2) +
  theme_bw()+
  #theme(axis.title.x = element_text(size = 4),
  #  axis.title.y = element_text(size = 4),
  #  axis.text.x = element_text(size = 4),
  #  axis.text.y = element_text(size =4))+
  xlim(c(-1600, 1600))+
  ylab("Count\n")+
  theme(legend.position = c(0.8, 0.8),
        legend.key.size = unit(0.15, "in"), 
        legend.title = element_text(size = 8),
        axis.title.x = element_text(size = 10)) + 
  xlab(expression(atop("\n Effect Size", "(Median Wet - Median Dry Rate; mg O"[2]*" kg"^-1*" H"^-1*")"))) 


ggsave("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/MEL_Effect_Histogram.png", width = 30, height = 30)


ggplot(effect_data, aes(x = Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H)) + 
  geom_histogram(aes(fill = after_stat(x))) + 
  facet_grid(~IN, scales = "free") + 
  scale_fill_gradient2(name = "Effect Size", limits = effect_limits, low = "firebrick2", mid = "goldenrod2",
                       high = "dodgerblue2", midpoint = (max(effect_limits)+min(effect_limits))/2) + 
  theme_bw()

ggplot(agg, aes(x = `total%_aggregates`)) + 
  geom_histogram(bins = 20) + 
  facet_grid(~Level, scales = "free") +
  theme_bw()


## Correlation Matrices ####

corr_data = effect_agg_wide %>% 
  select(c(Sample_Name, Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H,`53um_aggregates_wt%`, `125um_aggregates_wt%`, `250um_aggregates_wt%`, `2mm_aggregates_wt%`, `total%_aggregates`)) %>% 
  dplyr::rename(`53 um aggregates %`  = `53um_aggregates_wt%`) %>% 
  dplyr::rename(`125 um aggregates %` = `125um_aggregates_wt%`) %>% 
  dplyr::rename(`250 um aggregates %` = `250um_aggregates_wt%`) %>% 
  dplyr::rename(`2 mm aggregates %` = `2mm_aggregates_wt%`) %>% 
  dplyr::rename(`Total % Aggregates` = `total%_aggregates`) %>% 
  dplyr:: rename("Effect Size Respiration Rates (Wet - Dry)" = Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H) %>% 
  column_to_rownames("Sample_Name")

low_corr_data = effect_agg_wide %>% 
  filter(IN == "INL") %>% 
  select(c(Sample_Name, Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H,`53um_aggregates_wt%`, `125um_aggregates_wt%`, `250um_aggregates_wt%`, `2mm_aggregates_wt%`, `total%_aggregates`)) %>% 
  dplyr::rename(`53 um aggregates %`  = `53um_aggregates_wt%`) %>% 
  dplyr::rename(`125 um aggregates %` = `125um_aggregates_wt%`) %>% 
  dplyr::rename(`250 um aggregates %` = `250um_aggregates_wt%`) %>% 
  dplyr::rename(`2 mm aggregates %` = `2mm_aggregates_wt%`) %>% 
  dplyr::rename(`Total % Aggregates` = `total%_aggregates`) %>% 
  dplyr:: rename("Effect Size Respiration Rates (Wet - Dry)" = Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H)%>% 
  column_to_rownames("Sample_Name")
  
high_corr_data = effect_agg_wide %>% 
  filter(IN == "INH") %>% 
  select(c(Sample_Name, Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H,`53um_aggregates_wt%`, `125um_aggregates_wt%`, `250um_aggregates_wt%`, `2mm_aggregates_wt%`, `total%_aggregates`)) %>% 
  dplyr::rename(`53 um aggregates %`  = `53um_aggregates_wt%`) %>% 
  dplyr::rename(`125 um aggregates %` = `125um_aggregates_wt%`) %>% 
  dplyr::rename(`250 um aggregates %` = `250um_aggregates_wt%`) %>% 
  dplyr::rename(`2 mm aggregates %` = `2mm_aggregates_wt%`) %>% 
  dplyr::rename(`Total % Aggregates` = `total%_aggregates`) %>% 
  dplyr:: rename("Effect Size Respiration Rates (Wet - Dry)" = Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H)%>% 
  column_to_rownames("Sample_Name")

## Cube Root Transformations ####

cube_root <- function(x) sign(x) * (abs(x))^(1/3)

cube_corr_data = corr_data %>%
  dplyr::mutate(across(where(is.numeric), cube_root))

cube_high_data = high_corr_data %>%
  dplyr::mutate(across(where(is.numeric), cube_root))

cube_low_data = low_corr_data %>%
  dplyr::mutate(across(where(is.numeric), cube_root))

## Pearson Correlations ####
pearson = cor(corr_data, method = "pearson", use = "complete.obs")

pear_corr_all = corrplot(pearson, type = "upper", method = "number", tl.col = "black", tl.cex = 1.6, cl.cex = 1.25,  title = "Pearson Correlation")

## Surface Pearson Correlations ##

pearson_high <- cor(high_corr_data, method = "pearson", use = "complete.obs")

# calculate p-values for pearson correlation matrix, this should be the same as linear regression?
pear_high_p = as.data.frame(matrix(((rcorr(as.matrix(high_corr_data)))$P)[1, ], nrow = 1)) 
colnames(pear_high_p) = colnames(pearson_high)
rownames(pear_high_p) = c("Pearson_p_value")
pear_high_p = pear_high_p[,-1]

# Pull out Pearson correlation values
pear_effect_high = matrix(pearson_high[1, ], nrow = 1)
colnames(pear_effect_high) = colnames(pearson_high)
rownames(pear_effect_high) = c("Pearson Correlation Coefficient")

## Subsurface Pearson Correlations ##

pearson_low <- cor(low_corr_data, method = "pearson", use = "complete.obs")

# calculate p-values for pearson correlation matrix, this should be the same as linear regression?
pear_low_p = as.data.frame(matrix(((rcorr(as.matrix(low_corr_data)))$P)[1, ], nrow = 1)) 
colnames(pear_low_p) = colnames(pearson_low)
rownames(pear_low_p) = c("Pearson_p_value")
pear_low_p = pear_low_p[,-1]

# pull out Pearson correlation values 
pear_effect_low = matrix(pearson_low[1, ], nrow = 1)
colnames(pear_effect_low) = colnames(pearson_low)
rownames(pear_effect_low) = c("Pearson Correlation Coefficient")

pear_effect_high_df = as.data.frame(pear_effect_high) %>% 
  select(-c(`Effect Size Respiration Rates (Wet - Dry)`)) %>% 
  rbind(pear_high_p) %>% 
  rownames_to_column("Type") %>% 
  pivot_longer(-Type, names_to = "Aggregate Size", values_to = "Value") %>% 
  pivot_wider(
    names_from = Type,
    values_from = Value
  ) %>% 
  mutate(y = "Surface")

pear_effect_low_df = as.data.frame(pear_effect_low) %>% 
  select(-c(`Effect Size Respiration Rates (Wet - Dry)`)) %>% 
  rbind(pear_low_p) %>% 
  rownames_to_column("Type") %>% 
  pivot_longer(-Type, names_to = "Aggregate Size", values_to = "Value") %>% 
  pivot_wider(
    names_from = Type,
    values_from = Value
  ) %>% 
  mutate(y = "Subsurface")

color_palette_s = colorRampPalette(c("#B2182B", "#F7F7F7", "#2166AC"))(200)

agg_order = c("53 um aggregates %", "125 um aggregates %", "250 um aggregates %", "2 mm aggregates %", "Total % Aggregates")

max_p = 1
min_p = -1

pear_matrix = ggplot() +
  geom_tile(pear_effect_high_df, fill = "white", color = "black", mapping = aes(factor(`Aggregate Size`, level = agg_order), y)) +
  geom_text(pear_effect_high_df, mapping = aes(x = factor(`Aggregate Size`, level = agg_order), y = y, label = round(`Pearson Correlation Coefficient`, 2), color = `Pearson Correlation Coefficient`), size = 5, fontface = "bold") + 
  geom_tile(pear_effect_low_df, fill = "white", color = "black", mapping = aes(factor(`Aggregate Size`, level = agg_order), y)) +
  geom_text(pear_effect_low_df, mapping = aes(x = factor(`Aggregate Size`, level = agg_order), y = y, label = round(`Pearson Correlation Coefficient`, 2), color = `Pearson Correlation Coefficient`), size = 5, fontface = "bold") + 
  scale_color_gradientn(colors = color_palette_s, 
                        limit = c(min_p, max_p),
                        guide = "none") + 
  theme_bw() + 
  theme(aspect.ratio = 0.45, 
        axis.text.x = element_text(angle = 90, hjust = 0, size = 15), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 15))+#, 
  #axis.text.y = element_blank()) +
  scale_x_discrete(position = "top")

ggsave("C:/Github/ECA_Multireactor_Incubations/EL/Figures/Pearson_Correlation_Depths.png", plot = pear_matrix, width = 12, height = 6)

## Pearson Cube Root Correlations ##

## Cube Root Transformations ##

cube_corr_data = corr_data %>%
  dplyr::mutate(across(where(is.numeric), cube_root))

cube_high_data = high_corr_data %>%
  dplyr::mutate(across(where(is.numeric), cube_root))

cube_low_data = low_corr_data %>%
  dplyr::mutate(across(where(is.numeric), cube_root))

## Pearson Correlations ####
cube_pearson = cor(cube_corr_data, method = "pearson", use = "complete.obs")

cube_pear_corr_all = corrplot(cube_pearson, type = "upper", method = "number", tl.col = "black", tl.cex = 1.6, cl.cex = 1.25,  title = "Cube Root Pearson Correlation")

## Surface Pearson Correlations ##

cube_pearson_high <- cor(cube_high_data, method = "pearson", use = "complete.obs")

# calculate p-values for pearson correlation matrix, this should be the same as linear regression?
cube_pear_high_p = as.data.frame(matrix(((rcorr(as.matrix(cube_high_data)))$P)[1, ], nrow = 1)) 
colnames(cube_pear_high_p) = colnames(cube_pearson_high)
rownames(cube_pear_high_p) = c("Pearson_p_value")
cube_pear_high_p = cube_pear_high_p[,-1]

# Pull out Pearson correlation values
cube_pear_effect_high = matrix(cube_pearson_high[1, ], nrow = 1)
colnames(cube_pear_effect_high) = colnames(cube_pearson_high)
rownames(cube_pear_effect_high) = c("Pearson Correlation Coefficient")

## Subsurface Pearson Correlations ##

cube_pearson_low <- cor(cube_low_data, method = "pearson", use = "complete.obs")

# calculate p-values for pearson correlation matrix, this should be the same as linear regression?
cube_pear_low_p = as.data.frame(matrix(((rcorr(as.matrix(cube_low_data)))$P)[1, ], nrow = 1)) 
colnames(cube_pear_low_p) = colnames(cube_pearson_low)
rownames(cube_pear_low_p) = c("Pearson_p_value")
cube_pear_low_p = cube_pear_low_p[,-1]

# pull out Pearson correlation values 
cube_pear_effect_low = matrix(cube_pearson_low[1, ], nrow = 1)
colnames(cube_pear_effect_low) = colnames(cube_pearson_low)
rownames(cube_pear_effect_low) = c("Pearson Correlation Coefficient")

cube_pear_effect_high_df = as.data.frame(cube_pear_effect_high) %>% 
  select(-c(`Effect Size Respiration Rates (Wet - Dry)`)) %>% 
  rbind(cube_pear_high_p) %>% 
  rownames_to_column("Type") %>% 
  pivot_longer(-Type, names_to = "Aggregate Size", values_to = "Value") %>% 
  pivot_wider(
    names_from = Type,
    values_from = Value
  ) %>% 
  mutate(y = "Surface")

cube_pear_effect_low_df = as.data.frame(cube_pear_effect_low) %>% 
  select(-c(`Effect Size Respiration Rates (Wet - Dry)`)) %>% 
  rbind(cube_pear_low_p) %>% 
  rownames_to_column("Type") %>% 
  pivot_longer(-Type, names_to = "Aggregate Size", values_to = "Value") %>% 
  pivot_wider(
    names_from = Type,
    values_from = Value
  ) %>% 
  mutate(y = "Subsurface")

color_palette_s = colorRampPalette(c("#B2182B", "#F7F7F7", "#2166AC"))(200)

agg_order = c("53 um aggregates %", "125 um aggregates %", "250 um aggregates %", "2 mm aggregates %", "Total % Aggregates")

max_p = 1
min_p = -1

cube_pear_matrix = ggplot() +
  geom_tile(cube_pear_effect_high_df, fill = "white", color = "black", mapping = aes(factor(`Aggregate Size`, level = agg_order), y)) +
  geom_text(cube_pear_effect_high_df, mapping = aes(x = factor(`Aggregate Size`, level = agg_order), y = y, label = round(`Pearson Correlation Coefficient`, 2), color = `Pearson Correlation Coefficient`), size = 5, fontface = "bold") + 
  geom_tile(cube_pear_effect_low_df, fill = "white", color = "black", mapping = aes(factor(`Aggregate Size`, level = agg_order), y)) +
  geom_text(cube_pear_effect_low_df, mapping = aes(x = factor(`Aggregate Size`, level = agg_order), y = y, label = round(`Pearson Correlation Coefficient`, 2), color = `Pearson Correlation Coefficient`), size = 5, fontface = "bold") + 
  scale_color_gradientn(colors = color_palette_s, 
                        limit = c(min_p, max_p),
                        guide = "none") + 
  theme_bw() + 
  theme(aspect.ratio = 0.45, 
        axis.text.x = element_text(angle = 90, hjust = 0, size = 15), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 15))+#, 
  #axis.text.y = element_blank()) +
  scale_x_discrete(position = "top")

ggsave("C:/Github/ECA_Multireactor_Incubations/EL/Figures/Cube_Pearson_Correlation_Depths.png", plot = cube_pear_matrix, width = 12, height = 6)

## Spearman Correlations ####
spearman <- cor(corr_data, method = "spearman", use = "complete.obs")

corrplot(spearman,type = "upper", method = "number", tl.col = "black", tl.cex = 1.6, cl.cex = 1.25,  title = "Spearman Correlation")


## Surface Spearman Correlations ##

spearman_high <- cor(high_corr_data, method = "spearman", use = "complete.obs")

# calculate p-values for pearson correlation matrix, this should be the same as linear regression?
spear_high_p = as.data.frame(matrix(((rcorr(as.matrix(high_corr_data), type = "spearman"))$P)[1, ], nrow = 1))
colnames(spear_high_p) = colnames(spearman_high)
rownames(spear_high_p) = c("Spearman_p_value")
spear_high_p = spear_high_p[,-1]

# Spearman Correlation Coefficients
corr_effect_high = matrix(spearman_high[1, ], nrow = 1)
colnames(corr_effect_high) = colnames(spearman_high)
rownames(corr_effect_high) = c("Spearman Correlation Coefficient")

## Subsurface Spearman Correlations
spearman_low <- cor(low_corr_data, method = "spearman", use = "complete.obs")

# P values
spear_low_p = as.data.frame(matrix(((rcorr(as.matrix(low_corr_data), type = "spearman"))$P)[1, ], nrow = 1))
colnames(spear_low_p) = colnames(spearman_low)
rownames(spear_low_p) = c("Spearman_p_value")
spear_low_p = spear_low_p[,-1]

# Spearman Correlation Coefficients
corr_effect_low = matrix(spearman_low[1, ], nrow = 1)
colnames(corr_effect_low) = colnames(spearman_low)
rownames(corr_effect_low) = c("Spearman Correlation Coefficient")

corr_effect_high_df = as.data.frame(corr_effect_high) %>% 
  select(-c(`Effect Size Respiration Rates (Wet - Dry)`)) %>% 
  rbind(spear_high_p) %>% 
  rownames_to_column("Type") %>% 
  pivot_longer(-Type, names_to = "Aggregate Size", values_to = "Value") %>% 
  pivot_wider(
    names_from = Type,
    values_from = Value
  ) %>% 
  mutate(y = "Surface")

corr_effect_low_df = as.data.frame(corr_effect_low) %>% 
  select(-c(`Effect Size Respiration Rates (Wet - Dry)`)) %>% 
  rbind(spear_low_p) %>% 
  rownames_to_column("Type") %>% 
  pivot_longer(-Type, names_to = "Aggregate Size", values_to = "Value") %>% 
  pivot_wider(
    names_from = Type,
    values_from = Value
  ) %>% 
  mutate(y = "Subsurface") 


color_palette_s = colorRampPalette(c("#B2182B", "#F7F7F7", "#2166AC"))(200)

max_p = 1
min_p = -1

spear_matrix = ggplot() +
  geom_tile(corr_effect_high_df, fill = "white", color = "black", mapping = aes(factor(`Aggregate Size`, level = agg_order), y)) +
  geom_text(corr_effect_high_df, mapping = aes(x = factor(`Aggregate Size`, level = agg_order), y = y, label = round(`Spearman Correlation Coefficient`, 2), color = `Spearman Correlation Coefficient`), size = 5, fontface = "bold") + 
  geom_tile(corr_effect_low_df, fill = "white", color = "black", mapping = aes(factor(`Aggregate Size`, level = agg_order), y)) +
  geom_text(corr_effect_low_df, mapping = aes(x = factor(`Aggregate Size`, level = agg_order), y = y, label = round(`Spearman Correlation Coefficient`, 2), color = `Spearman Correlation Coefficient`), size = 5, fontface = "bold") + 
  scale_color_gradientn(colors = color_palette_s, 
                        limit = c(min_p, max_p),
                        guide = "none") + 
  theme_bw() + 
  theme(aspect.ratio = 0.45, 
        axis.text.x = element_text(angle = 90, hjust = 0, size = 15), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())+#, 
        #axis.text.y = element_blank()) +
  scale_x_discrete(position = "top")

ggsave("C:/Github/ECA_Multireactor_Incubations/EL/Figures/Spearman_Correlation_Depths.png", plot = spear_matrix, width = 12, height = 6)

## Aggregate Scatter Plots  ####

## Untransformed Data ##

effect_agg$agg_size = factor(effect_agg$agg_size, levels = c('53um_aggregates_wt%', "125um_aggregates_wt%", "250um_aggregates_wt%", "2mm_aggregates_wt%", "total%_aggregates"))

agg_depth_scatter = ggplot(effect_agg, mapping = aes(x = perc, y = Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H, color = IN)) + 
  geom_point() + 
  theme_bw() +
  facet_wrap(agg_size~., scales = "free", nrow = 1) + 
  stat_poly_line(data = effect_agg, se = FALSE) + 
  stat_correlation(method = "pearson", label.y = "middle", label.x = "right", size = 3, use_label(c("R2", "P")))

ggsave("C:/Github/ECA_Multireactor_Incubations/EL/Figures/EL_Depth_Agg_Scatter_Plots_No_Transform.png", plot = agg_depth_scatter, width = 20, height = 10)

## Cube Root Transformed Data ##

cube_effect_agg = effect_agg %>% 
  dplyr::mutate(across(where(is.numeric), cube_root))

cube_effect_agg$agg_size = factor(cube_effect_agg$agg_size, levels = c('53um_aggregates_wt%', "125um_aggregates_wt%", "250um_aggregates_wt%", "2mm_aggregates_wt%", "total%_aggregates"))

cube_agg_depth_scatter = ggplot(cube_effect_agg, mapping = aes(x = perc, y = Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H, color = IN)) + 
  geom_point() + 
  theme_bw() +
  facet_wrap(agg_size~., scales = "free", nrow = 1) + 
  stat_poly_line(data = cube_effect_agg, se = FALSE) + 
  stat_correlation(method = "pearson", label.y = "middle", label.x = "right", size = 3, use_label(c("R2", "P")))

ggsave("C:/Github/ECA_Multireactor_Incubations/EL/Figures/EL_Depth_Agg_Scatter_Plots_Cube_Transform.png", plot = cube_agg_depth_scatter, width = 20, height = 10)

## Cube Effect vs. Untransformed Aggregate %'s ##

cube_effect = effect_agg %>% 
  mutate(cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H = cube_root(Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H)) 

cube_effect$agg_size = factor(cube_effect$agg_size, levels = c('53um_aggregates_wt%', "125um_aggregates_wt%", "250um_aggregates_wt%", "2mm_aggregates_wt%", "total%_aggregates"))


cube_effect_agg_depth = ggplot(cube_effect, mapping = aes(x = perc, y = cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H, color = IN)) + 
  geom_point() + 
  theme_bw() +
  facet_wrap(agg_size~., scales = "free", nrow = 1) + 
  stat_poly_line(data = cube_effect, se = FALSE) + 
  stat_correlation(method = "pearson", label.y = "middle", label.x = "right", size = 3, use_label(c("R2", "P")))


ggsave("C:/Github/ECA_Multireactor_Incubations/EL/Figures/EL_Depth_Agg_Scatter_Plots.png", plot = cube_effect_agg_depth, width = 20, height = 10)

