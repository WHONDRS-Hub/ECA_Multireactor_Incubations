#### ECA LASSO ####
library(tidyverse)
library(glmnet)

rm(list=ls());graphics.off()

## Read in all data
# remove NEON sites - 052, 053, 057
# remove samples with no respiration rate
# Total: 512 samples

all_data = read.csv("C:/GitHub/ECA_Multireactor_Incubations/Data/Cleaned Data/All_ECA_Data_03-06-2024.csv")  %>% 
  filter(!grepl("052", Sample_Name)) %>% 
  filter(!grepl("053", Sample_Name)) %>% 
  filter(!grepl("057", Sample_Name)) %>% 
  filter(Respiration_Rate_mg_DO_per_L_per_H != -9999)

#remove site 023 because no moisture data
all_data_mg_L = all_data %>% 
  select(-c(Respiration_Rate_mg_DO_per_kg_per_H, Fe_mg_per_kg, Sample_ID, INC)) %>%
  filter(!grepl("023", Sample_Name)) %>% 
  mutate(Fe_mg_per_L = if_else(grepl("SFE_Below", Fe_mg_per_L), "0.002", Fe_mg_per_L)) %>% 
  mutate(Fe_mg_per_L = as.numeric(Fe_mg_per_L)) %>% 
  mutate(Respiration_Rate_mg_DO_per_L_per_H = abs(Respiration_Rate_mg_DO_per_L_per_H))

ggplot(all_data_mg_L, aes(x = Respiration_Rate_mg_DO_per_L_per_H)) + 
  geom_histogram()

### TRANSFORM DATA

# Write functions for log/cube transformations 
half_min_values = function(x) {
  min_val <- min(x)
  if (min_val == 0) {
    min_val <- sort(x[x != 0])[1]
  }
  return(0.5 * min_val)
}

cube_root <- function(x) sign(x) * (abs(x))^(1/3)

#Log transform all data
log_lasso = all_data_mg_L %>% 
  select(-c(Lost_Gravimetric_Water)) %>% 
  mutate(across(where(is.numeric), ~. + half_min_values(.))) %>% 
  mutate(across(where(is.numeric), ~log10(.))) %>% 
  column_to_rownames("Sample_Name")

z_log_lasso = log_lasso %>% 
  mutate(across(where(is.numeric), function(x) (x - mean(x)) / sd(x)))

z_log_resp <- z_log_lasso$Respiration_Rate_mg_DO_per_L_per_H

#Try all data
#z_log_pred <- data.matrix(z_log_lasso[, c('Fe_mg_per_L', 'Percent_Fine_Sand', 'Percent_Coarse_Sand', 'Percent_Tot_Sand', 'Percent_Silt', 'Percent_Clay', 'mean_ssa', 'SpC', 'Temp', 'pH', 'Initial_Gravimetric_Water', 'Final_Gravimetric_Water')])

## Log FS, SSA, Fe, Resp, Moi

z_log_pred <- data.matrix(z_log_lasso[, c("Fe_mg_per_L", "ATP_nanomol_per_L",  "Percent_Fine_Sand", "mean_ssa", "SpC", "Temp", "pH",  "Final_Gravimetric_Water")])

cv_model <- cv.glmnet(z_log_pred, z_log_resp, alpha = 1)

best_lambda <- cv_model$lambda.min
best_lambda

plot(cv_model)

best_model <- glmnet(z_log_pred, z_log_resp, alpha = 1, lambda = best_lambda)
coef(best_model)

resp_predict <- predict(best_model, s = best_lambda, newx = z_log_pred)

sst <- sum((z_log_resp - mean(z_log_resp))^2)
sse <- sum((resp_predict - z_log_resp)^2)

rsq = 1 - sse/sst

#rsq. 0.63
rsq

# Wet Log Lasso
wet_log_lasso = log_lasso %>% 
  rownames_to_column("Sample_Name") %>% 
  filter(grepl("W", Sample_Name))%>% 
  mutate(across(where(is.numeric), function(x) (x - mean(x)) / sd(x)))

z_log_resp_wet <- wet_log_lasso$Respiration_Rate_mg_DO_per_L_per_H

#Try all data
#z_log_pred <- data.matrix(z_log_lasso[, c('Fe_mg_per_L', 'Percent_Fine_Sand', 'Percent_Coarse_Sand', 'Percent_Tot_Sand', 'Percent_Silt', 'Percent_Clay', 'mean_ssa', 'SpC', 'Temp', 'pH', 'Initial_Gravimetric_Water', 'Final_Gravimetric_Water')])

## Log FS, SSA, Fe, Resp, Moi

z_log_pred_wet <- data.matrix(wet_log_lasso[, c("Fe_mg_per_L", "ATP_nanomol_per_L",  "Percent_Fine_Sand", "mean_ssa", "SpC", "Temp", "pH",  "Final_Gravimetric_Water")])

cv_model <- cv.glmnet(z_log_pred_wet, z_log_resp_wet, alpha = 1)

best_lambda <- cv_model$lambda.min
best_lambda

plot(cv_model)

best_model <- glmnet(z_log_pred_wet, z_log_resp_wet, alpha = 1, lambda = best_lambda)
coef(best_model)

resp_predict <- predict(best_model, s = best_lambda, newx = z_log_pred_wet)

sst <- sum((z_log_resp_wet - mean(z_log_resp_wet))^2)
sse <- sum((resp_predict - z_log_resp_wet)^2)

rsq = 1 - sse/sst

#rsq. 0.62
rsq

## Dry log lasso
dry_log_lasso = log_lasso %>% 
  rownames_to_column("Sample_Name") %>% 
  filter(grepl("D", Sample_Name))%>% 
  mutate(across(where(is.numeric), function(x) (x - mean(x)) / sd(x)))

z_log_resp_dry <- dry_log_lasso$Respiration_Rate_mg_DO_per_L_per_H

#Try all data
#z_log_pred <- data.matrix(z_log_lasso[, c('Fe_mg_per_L', 'Percent_Fine_Sand', 'Percent_Coarse_Sand', 'Percent_Tot_Sand', 'Percent_Silt', 'Percent_Clay', 'mean_ssa', 'SpC', 'Temp', 'pH', 'Initial_Gravimetric_Water', 'Final_Gravimetric_Water')])

## Log FS, SSA, Fe, Resp, Moi

z_log_pred_dry <- data.matrix(dry_log_lasso[, c("Fe_mg_per_L", "ATP_nanomol_per_L",  "Percent_Fine_Sand", "mean_ssa", "SpC", "Temp", "pH",  "Final_Gravimetric_Water")])

cv_model <- cv.glmnet(z_log_pred_dry, z_log_resp_dry, alpha = 1)

best_lambda <- cv_model$lambda.min
best_lambda

plot(cv_model)

best_model <- glmnet(z_log_pred_dry, z_log_resp_dry, alpha = 1, lambda = best_lambda)
coef(best_model)

resp_predict <- predict(best_model, s = best_lambda, newx = z_log_pred_dry)

sst <- sum((z_log_resp_dry - mean(z_log_resp_dry))^2)
sse <- sum((resp_predict - z_log_resp_dry)^2)

rsq = 1 - sse/sst

#rsq. 0.56
rsq

## Cube Root LASSO
cube_lasso = all_data_mg_L %>% 
  mutate(across(where(is.numeric), cube_root))%>% 
  column_to_rownames("Sample_Name")

z_cube_lasso = cube_lasso %>% 
  mutate(across(where(is.numeric), function(x) ((x - mean(x)) / sd(x))))

z_cube_resp <- z_cube_lasso$Respiration_Rate_mg_DO_per_L_per_H

#Try all data
z_cube_pred <- data.matrix(z_cube_lasso[, c('Fe_mg_per_L', 'Percent_Fine_Sand', 'Percent_Coarse_Sand', 'Percent_Tot_Sand', 'Percent_Silt', 'Percent_Clay', 'mean_ssa', 'SpC', 'Temp', 'pH', 'Initial_Gravimetric_Water', 'Final_Gravimetric_Water')])

## Log FS, SSA, Fe, Resp, Moi

#z_log_pred <- data.matrix(z_log_lasso[, c("Fe_mg_per_L", "ATP_nanomol_per_L",  "Percent_Fine_Sand", "mean_ssa", "SpC", "Temp", "pH",  "Final_Gravimetric_Water")])

cv_model <- cv.glmnet(z_cube_pred, z_cube_resp, alpha = 1)

best_lambda <- cv_model$lambda.min
best_lambda

plot(cv_model)

best_model <- glmnet(z_cube_pred, z_cube_resp, alpha = 1, lambda = best_lambda)
coef(best_model)

resp_predict <- predict(best_model, s = best_lambda, newx = z_cube_pred)

sst <- sum((z_cube_resp - mean(z_cube_resp))^2)
sse <- sum((resp_predict - z_cube_resp)^2)

rsq = 1 - sse/sst

#rsq. 0.63
rsq

# Wet Cube Lasso
wet_cube_lasso = cube_lasso %>% 
  rownames_to_column("Sample_Name") %>% 
  filter(grepl("W", Sample_Name))%>% 
  mutate(across(where(is.numeric), function(x) (x - mean(x)) / sd(x)))

z_cube_resp_wet <- wet_cube_lasso$Respiration_Rate_mg_DO_per_L_per_H

#Try all data
#z_log_pred <- data.matrix(z_log_lasso[, c('Fe_mg_per_L', 'Percent_Fine_Sand', 'Percent_Coarse_Sand', 'Percent_Tot_Sand', 'Percent_Silt', 'Percent_Clay', 'mean_ssa', 'SpC', 'Temp', 'pH', 'Initial_Gravimetric_Water', 'Final_Gravimetric_Water')])

## Log FS, SSA, Fe, Resp, Moi

z_cube_pred_wet <- data.matrix(wet_cube_lasso[, c("Fe_mg_per_L", "ATP_nanomol_per_L",  "Percent_Fine_Sand", "mean_ssa", "SpC", "Temp", "pH",  "Final_Gravimetric_Water")])

cv_model <- cv.glmnet(z_cube_pred_wet, z_cube_resp_wet, alpha = 1)

best_lambda <- cv_model$lambda.min
best_lambda

plot(cv_model)

best_model <- glmnet(z_cube_pred_wet, z_cube_resp_wet, alpha = 1, lambda = best_lambda)
coef(best_model)

resp_predict <- predict(best_model, s = best_lambda, newx = z_cube_pred_wet)

sst <- sum((z_cube_resp_wet - mean(z_cube_resp_wet))^2)
sse <- sum((resp_predict - z_cube_resp_wet)^2)

rsq = 1 - sse/sst

#rsq. 0.65
rsq

## Dry log lasso
dry_cube_lasso = cube_lasso %>% 
  rownames_to_column("Sample_Name") %>% 
  filter(grepl("D", Sample_Name))%>% 
  mutate(across(where(is.numeric), function(x) (x - mean(x)) / sd(x)))

z_cube_resp_dry <- dry_cube_lasso$Respiration_Rate_mg_DO_per_L_per_H

#Try all data
#z_log_pred <- data.matrix(z_log_lasso[, c('Fe_mg_per_L', 'Percent_Fine_Sand', 'Percent_Coarse_Sand', 'Percent_Tot_Sand', 'Percent_Silt', 'Percent_Clay', 'mean_ssa', 'SpC', 'Temp', 'pH', 'Initial_Gravimetric_Water', 'Final_Gravimetric_Water')])

## Log FS, SSA, Fe, Resp, Moi

z_cube_pred_dry <- data.matrix(dry_cube_lasso[, c("Fe_mg_per_L", "ATP_nanomol_per_L",  "Percent_Fine_Sand", "mean_ssa", "SpC", "Temp", "pH",  "Final_Gravimetric_Water")])

cv_model <- cv.glmnet(z_cube_pred_dry, z_cube_resp_dry, alpha = 1)

best_lambda <- cv_model$lambda.min
best_lambda

plot(cv_model)

best_model <- glmnet(z_cube_pred_dry, z_cube_resp_dry, alpha = 1, lambda = best_lambda)
coef(best_model)

resp_predict <- predict(best_model, s = best_lambda, newx = z_cube_pred_dry)

sst <- sum((z_cube_resp_dry - mean(z_cube_resp_dry))^2)
sse <- sum((resp_predict - z_cube_resp_dry)^2)

rsq = 1 - sse/sst

#rsq. 0.58
rsq

## Plot all log transformations vs. non log transformed
# 
# ggplot(gather(all_data_mg_L, cols, value), aes(x = value)) + 
#   geom_histogram() + 
#   facet_wrap(.~cols, scales = 'free_x')
# 
# ggplot(gather(log_lasso, cols, value), aes(x = value)) + 
#   geom_histogram() + 
#   facet_wrap(.~cols, scales = 'free_x')
# 
# ggplot(gather(cube_lasso, cols, value), aes(x = value)) + 
#   geom_histogram() + 
#   facet_wrap(.~cols, scales = 'free_x')



## Read in means for effect size

atp_means = read.csv("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/ATP/03_ProcessedData/EC_ATP_Summary_ReadyForBoye_03-05-2024.csv") %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "ATP", "INC"))

fe_means = read.csv("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Fe/03_ProcessedData/EC_SFE_Summary_ReadyForBoye_03-05-2024.csv") %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "SFE", "INC"))

resp_means = read.csv("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/INC/03_ProcessedData/ECA_Sediment_Incubations_Respiration_Rates_Summary_ReadyForBoye_2024-03-05.csv")

moi_means = read.csv("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/INC/03_ProcessedData/ECA_Drying_Masses_Summary_ReadyForBoye_01-29-2024.csv") %>% 
  separate(Sample_Name, c("Sample_Name", "Rep"), sep = -1) %>% 
  mutate(Initial_Gravimetric_Water = Initial_Water_Mass_g/Dry_Sediment_Mass_g) %>% 
  mutate(Final_Gravimetric_Water = Final_Water_Mass_g/Dry_Sediment_Mass_g) %>% 
  mutate(Lost_Gravimetric_Water = Initial_Gravimetric_Water - Final_Gravimetric_Water) %>% 
  group_by(Sample_Name) %>% 
  summarise(across(c(Initial_Gravimetric_Water:Lost_Gravimetric_Water), mean))

grain <- read.csv("C:/Github/ECA_Multireactor_Incubations/Data/v3_CM_SSS_Sediment_Grain_Size.csv", skip = 2, header = TRUE)

grain_all <- grain %>% 
  filter(!row_number() %in% c(1:11)) %>% 
  dplyr::select(-c(Field_Name, Material)) %>% 
  filter(!grepl("SSS", Sample_Name)) %>% 
  filter(Sample_Name != "") %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "CM", "EC")) %>% 
  mutate(Sample_ID = str_remove(Sample_Name, "_GRN")) %>% 
  mutate_at(c("Percent_Fine_Sand", "Percent_Med_Sand", "Percent_Coarse_Sand", "Percent_Tot_Sand", "Percent_Silt", "Percent_Clay"), as.numeric) %>% 
  select(-c(IGSN, Methods_Deviation, Sample_Name))

ssa <- read.csv("C:/GitHub/ECA_Multireactor_Incubations/Data/CM_SSS_Sediment_Specific_Surface_Area.csv", skip = 2, header = TRUE)

ssa_clean = ssa %>% 
  filter(!row_number() %in% c(1:11)) %>% 
  dplyr::select(-c(Field_Name, Material, IGSN)) %>% 
  filter(!grepl("SSS", Sample_Name)) %>% 
  filter(Sample_Name != "") %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "CM", "EC")) %>% 
  separate(Sample_Name, c("Sample_ID", "Rep"), sep = -6)

mean_ssa <- ssa_clean %>% 
  filter(Specific_Surface_Area_m2_per_g != -9999) %>% 
  filter(!grepl("Negative", Specific_Surface_Area_m2_per_g)) %>%
  mutate(Specific_Surface_Area_m2_per_g = as.numeric(Specific_Surface_Area_m2_per_g)) %>%
  group_by(Sample_ID) %>%
  summarise(mean_ssa = mean(Specific_Surface_Area_m2_per_g, na.rm = TRUE))

means = left_join(resp_means, moi_means, by = "Sample_Name") %>% 
  left_join(atp_means) %>% 
  left_join(fe_means) %>% 
  select(-contains("SD")) %>% 
  separate(Sample_Name, c("Sample_ID", "Rep"), sep = -6, remove = FALSE) %>% 
  separate(Rep, c("INC", "Treat"), sep = -1) %>% 
  left_join(grain_all, by = "Sample_ID") %>% 
  left_join(mean_ssa, by = "Sample_ID") %>% 
  filter(!grepl("052", Sample_ID)) %>% 
  filter(!grepl("053", Sample_ID)) %>% 
  filter(!grepl("057", Sample_ID))

best_means = means %>% 
  select(-contains("WithOutliers")) %>% 
  select(-c(Material, INC, Sample_Name)) %>% 
  select(-contains("kg")) %>% 
  mutate(Mean_OutliersRemoved_Respiration_Rate_mg_DO_per_L_per_H = abs(Mean_OutliersRemoved_Respiration_Rate_mg_DO_per_L_per_H))

best_effect = best_means %>% 
  filter(!grepl("011", Sample_ID)) %>% 
  filter(!grepl("012", Sample_ID)) %>% 
group_by(Sample_ID) %>% 
  mutate(across(c(Mean_SpC:Mean_Fe_mg_per_L), ~. [Treat == "W"] - .[Treat  == "D"])) %>% 
  distinct(Sample_ID, .keep_all = TRUE) 
 
cube_best_effect = best_effect %>% 
  mutate(across(where(is.numeric), cube_root)) %>% 
  select(-c(Treat)) %>% 
  column_to_rownames("Sample_ID")

z_cube_best_effect = cube_best_effect %>% 
  mutate(across(where(is.numeric), function(x) ((x - mean(x)) / sd(x))))


## Effect LASSO 
eff <- z_cube_best_effect$Mean_OutliersRemoved_Respiration_Rate_mg_DO_per_L_per_H

#Try all data
#z_cube_effect_pred <- data.matrix(z_cube_best_effect[, c('Mean_Fe_mg_per_L', 'Percent_Fine_Sand', 'Percent_Coarse_Sand', 'Percent_Tot_Sand', 'Percent_Silt', 'Percent_Clay', 'mean_ssa', 'Mean_SpC', 'Mean_Temp', 'Mean_pH', 'Initial_Gravimetric_Water', 'Final_Gravimetric_Water')])

## Log FS, SSA, Fe, Resp, Moi

z_cube_effect_pred <- data.matrix(z_cube_best_effect[, c("Mean_Fe_mg_per_L", "Mean_ATP_nanomol_per_L",  "Percent_Fine_Sand", "mean_ssa", "Mean_SpC", "Mean_Temp", "Mean_pH",  "Final_Gravimetric_Water")])

cv_model <- cv.glmnet(z_cube_effect_pred, eff, alpha = 1)

best_lambda <- cv_model$lambda.min
best_lambda

plot(cv_model)

best_model <- glmnet(z_cube_effect_pred, eff, alpha = 1, lambda = best_lambda)
coef(best_model)

eff_predict <- predict(best_model, s = best_lambda, newx = z_cube_effect_pred)

sst <- sum((eff - mean(eff))^2)
sse <- sum((eff_predict - eff)^2)

rsq = 1 - sse/sst

#RSq 0.54
rsq



