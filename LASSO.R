#### ECA LASSO ####
library(tidyverse)
library(glmnet)
library(factoextra)
library(janitor)

rm(list=ls());graphics.off()

## Read in and clean all data ####
# remove NEON sites - 052, 053, 057
# remove samples with no respiration rate
# Total: 512 samples

icon_resp = read.csv("C:/GitHub/ECA_Multireactor_Incubations/Data/v2_CM_SSS_Sediment_Sample_Data_Summary.csv") %>% 
  filter(!row_number() %in% c(1, 3:13)) %>% 
  janitor::row_to_names(row_number = 1) %>% 
  dplyr::select(c(Sample_Name,Mean_Normalized_Respiration_Rate_mg_DO_per_H_per_L_sediment)) %>% 
  filter(!grepl("SSS", Sample_Name)) %>% 
  filter(Sample_Name != "") %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "CM", "EC")) %>% 
  mutate(Sample_ID = str_remove(Sample_Name, "_Sediment")) %>% 
  filter(Mean_Normalized_Respiration_Rate_mg_DO_per_H_per_L_sediment != -9999) %>% 
  mutate(Mean_Normalized_Respiration_Rate_mg_DO_per_H_per_L_sediment = abs(as.numeric(Mean_Normalized_Respiration_Rate_mg_DO_per_H_per_L_sediment))) %>% 
  select(-c(Sample_Name))

all_data = read.csv("C:/GitHub/ECA_Multireactor_Incubations/Data/Cleaned Data/All_ECA_Data_03-06-2024.csv")  %>% 
  filter(!grepl("052", Sample_Name)) %>% 
  filter(!grepl("053", Sample_Name)) %>% 
  filter(!grepl("057", Sample_Name)) %>% 
  filter(Respiration_Rate_mg_DO_per_L_per_H != -9999) %>% 
  left_join(icon_resp, by = "Sample_ID")

#remove site 023 because no moisture data
all_data_mg_L = all_data %>% 
  select(-c(Sample_ID, INC)) %>%
  filter(!grepl("023", Sample_Name)) %>% 
  mutate(Fe_mg_per_L = if_else(grepl("SFE_Below", Fe_mg_per_L), "0.002", Fe_mg_per_L)) %>% 
  mutate(Fe_mg_per_kg = if_else(grepl("SFE_Below", Fe_mg_per_kg), "0.006", Fe_mg_per_kg)) %>% 
  mutate(Fe_mg_per_L = as.numeric(Fe_mg_per_L)) %>% 
  mutate(Respiration_Rate_mg_DO_per_L_per_H = abs(Respiration_Rate_mg_DO_per_L_per_H)) %>% 
    mutate(Fe_mg_per_kg = as.numeric(Fe_mg_per_kg)) %>% 
    mutate(Respiration_Rate_mg_DO_per_kg_per_H = abs(Respiration_Rate_mg_DO_per_kg_per_H)) %>% 
  rename(ICON_resp = Mean_Normalized_Respiration_Rate_mg_DO_per_H_per_L_sediment) %>% 
  mutate(Treat = ifelse(grepl("W", Sample_Name), "wet", "dry")) #%>% 
  #mutate(type = ifelse(Respiration_Rate_mg_DO_per_L_per_H > 90, "theoretical", "real"))

ggplot(all_data_mg_L, aes(x = ATP_picomol_per_g, y = Respiration_Rate_mg_DO_per_kg_per_H)) + 
  geom_point(aes(color = treat)) 

# all_data_th = all_data_mg_L %>% 
#   mutate(treat = ifelse(grepl("W", Sample_Name), "wet", "dry")) %>% 
#   mutate(type = ifelse(Respiration_Rate_mg_DO_per_L_per_H > 90, "theoretical", "real"))

ggplot(all_data_th, aes(x = Final_Gravimetric_Water, y = Respiration_Rate_mg_DO_per_L_per_H)) +
  geom_point(aes(color = treat))

ggplot(all_data_th, aes(y = Respiration_Rate_mg_DO_per_L_per_H, x = Final_Gravimetric_Water, color = treat)) +
  geom_point()

#ggplot(all_data_mg_L, aes(x = Respiration_Rate_mg_DO_per_L_per_H)) + 
 #geom_histogram()

### TRANSFORM DATA ####

# Write functions for log/cube transformations 
half_min_values = function(x) {
  min_val <- min(x)
  if (min_val == 0) {
    min_val <- sort(x[x != 0])[1]
  }
  return(0.5 * min_val)
}

cube_root <- function(x) sign(x) * (abs(x))^(1/3)

## CUBE TRANSFORMATIONS ####

## All mg/L values 
## Change to get w/d lasso
## theoretical vs. real: mg/L > 4.5


cube_lasso = all_data_mg_L %>% 
  mutate(across(where(is.numeric), cube_root)) %>%
  #filter(Respiration_Rate_mg_DO_per_L_per_H < 4.5) %>% 
  #filter(grepl("W", Sample_Name)) #%>% 
  #mutate(Treat = ifelse(grepl("W", Sample_Name), "Wet", "Dry")) %>% 
  column_to_rownames("Sample_Name") 

ggplot(cube_lasso, aes(x = Final_Gravimetric_Water, y = Respiration_Rate_mg_DO_per_kg_per_H)) +
  geom_point(aes(color = Treat))+
  geom_smooth()

ggplot(all_data_mg_L, aes(x = Final_Gravimetric_Water, y = Respiration_Rate_mg_DO_per_kg_per_H)) +
  geom_point(aes(color = Treat))+
  geom_smooth()




## Scale cube root values 
z_cube_lasso = cube_lasso %>% 
  mutate(across(where(is.numeric), function(x) ((x - mean(x)) / sd(x))))

## ALL DATA LASSO ####
## Set response variable
z_cube_resp <- z_cube_lasso$Respiration_Rate_mg_DO_per_kg_per_H

## Set predictor variables 

#Try all data
#z_cube_pred <- data.matrix(z_cube_lasso[, c('Fe_mg_per_L', 'Percent_Fine_Sand', 'Percent_Med_Sand', 'Percent_Coarse_Sand', 'Percent_Tot_Sand', 'Percent_Silt', 'Percent_Clay', 'mean_ssa', 'SpC', 'Temp', 'pH', 'Initial_Gravimetric_Water', 'Final_Gravimetric_Water')])

## Selected variables: FS, SSA, Fe, ICON Resp, Moi, ATP, SpC, Temp, pH

z_cube_pred <- data.matrix(z_cube_lasso[, c("Fe_mg_per_kg",  "ATP_picomol_per_g",  "Percent_Fine_Sand", "mean_ssa", "SpC", "Temp", "pH",  "Final_Gravimetric_Water", "ICON_resp")])

## Start LASSO
#alpha = 1 is for LASSO regression
cv_model <- cv.glmnet(z_cube_pred, z_cube_resp, alpha = 1)

#Find lambda for lowest MSE using k-fold cross-validation
best_lambda <- cv_model$lambda.min
best_lambda

plot(cv_model)

#Run model with the best lambda value and get coefficient estimates from model 
best_model <- glmnet(z_cube_pred, z_cube_resp, alpha = 1, lambda = best_lambda)
coef(best_model)

#Make predictions to calculate R sq.
resp_predict <- predict(best_model, s = best_lambda, newx = z_cube_pred)

sst <- sum((z_cube_resp - mean(z_cube_resp))^2)
sse <- sum((resp_predict - z_cube_resp)^2)

rsq = 1 - sse/sst

rsq

## Look at model residuals
residuals = z_cube_resp - resp_predict

par(mfrow=c(2,2)) # Set up a 2x2 grid of plots
hist(residuals, main="Histogram of Residuals")
qqnorm(residuals, main="QQ Plot of Residuals")
qqline(residuals)
plot(density(residuals), main="Density Plot of Residuals")
plot(residuals ~ resp_predict)

dev.off()

## EFFECT SIZE DATA ####

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
  filter(!grepl("023", Sample_Name)) %>% 
  group_by(Sample_Name) %>% 
  summarise(across(c(Initial_Gravimetric_Water:Lost_Gravimetric_Water), mean))

grain <- read.csv("C:/Github/ECA_Multireactor_Incubations/Data/v3_CM_SSS_Sediment_Grain_Size.csv", skip = 2, header = TRUE) %>% 
  filter(!row_number() %in% c(1:11)) %>% 
  dplyr::select(-c(Field_Name, Material)) %>% 
  filter(!grepl("SSS", Sample_Name)) %>% 
  filter(Sample_Name != "") %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "CM", "EC")) %>% 
  mutate(Sample_ID = str_remove(Sample_Name, "_GRN")) %>% 
  mutate_at(c("Percent_Fine_Sand", "Percent_Med_Sand", "Percent_Coarse_Sand", "Percent_Tot_Sand", "Percent_Silt", "Percent_Clay"), as.numeric) %>% 
  select(-c(IGSN, Methods_Deviation, Sample_Name))

ssa <- read.csv("C:/GitHub/ECA_Multireactor_Incubations/Data/CM_SSS_Sediment_Specific_Surface_Area.csv", skip = 2, header = TRUE) %>% 
  filter(!row_number() %in% c(1:11)) %>% 
  dplyr::select(-c(Field_Name, Material, IGSN)) %>% 
  filter(!grepl("SSS", Sample_Name)) %>% 
  filter(Sample_Name != "") %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "CM", "EC")) %>% 
  separate(Sample_Name, c("Sample_ID", "Rep"), sep = -6) %>% 
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
  left_join(grain, by = "Sample_ID") %>% 
  left_join(ssa, by = "Sample_ID") %>% 
  filter(!grepl("052", Sample_ID)) %>% 
  filter(!grepl("053", Sample_ID)) %>% 
  filter(!grepl("057", Sample_ID)) %>% 
  filter(!grepl("023", Sample_ID)) %>% 
  select(-c(Material, INC, Sample_Name)) %>% 
  mutate(Mean_OutliersRemoved_Respiration_Rate_mg_DO_per_L_per_H = abs(Mean_OutliersRemoved_Respiration_Rate_mg_DO_per_L_per_H)) %>% 
  mutate(Mean_OutliersRemoved_Respiration_Rate_mg_DO_per_kg_per_H = abs(Mean_OutliersRemoved_Respiration_Rate_mg_DO_per_kg_per_H)) %>% 
  mutate(Mean_WithOutliers_Respiration_Rate_mg_DO_per_L_per_H = abs(Mean_WithOutliers_Respiration_Rate_mg_DO_per_L_per_H))%>% 
  mutate(Mean_WithOutliers_Respiration_Rate_mg_DO_per_kg_per_H = abs(Mean_WithOutliers_Respiration_Rate_mg_DO_per_kg_per_H))

## EFFECT SIZE CLEANING ####

## Calculate Effect Size 
best_effect = means %>% 
  filter(!grepl("011", Sample_ID)) %>% 
  filter(!grepl("012", Sample_ID)) %>% 
  group_by(Sample_ID) %>% 
  mutate(across(c(Mean_SpC:Mean_Fe_mg_per_L), ~. [Treat == "W"] - .[Treat  == "D"])) %>% 
  rename_with(.cols = c(Mean_SpC:Mean_Fe_mg_per_L), .fn = ~ paste0("diff_", .x)) %>% 
  distinct(Sample_ID, .keep_all = TRUE) %>% 
  left_join(means, by = c("Sample_ID", "Treat")) %>% 
  mutate(grav_cat = ifelse(Final_Gravimetric_Water > 0.4, "VM", ifelse(Final_Gravimetric_Water > 0.25, "M", ifelse(Final_Gravimetric_Water > 0.15, "SM", ifelse(Final_Gravimetric_Water > 0.05, "D", "VD")))))


limits = c(-235, 77)

ggplot(best_effect, aes(x = Percent_Fine_Sand.x, y =diff_Mean_OutliersRemoved_Respiration_Rate_mg_DO_per_kg_per_H)) + 
  geom_point(aes(color = Percent_Fine_Sand.x)) +
  #geom_point(aes(color = diff_Mean_ATP_picomol_per_g)) +
  #geom_smooth() +
  scale_color_gradient2(limits = limits, low = "firebrick2", mid = "goldenrod2", high = "dodgerblue2", midpoint = (max(limits)+min(limits))/2) +
  theme_bw() +
  xlab("% Fine Sand") +
  ylab("Effect Size") +
  #labs(color = "ATP")


## Cube Effect Size
# To get theoretical, effect size > 3.5
cube_best_effect = best_effect %>% 
  mutate(across(where(is.numeric), cube_root)) %>% 
  rename_with(.cols = c(diff_Mean_SpC:mean_ssa), .fn = ~ paste0("cube_", .x)) %>% 
  select(-c(Treat)) %>% 
  #filter(Mean_WithOutliers_Respiration_Rate_mg_DO_per_L_per_H < 3.5) %>% 
  column_to_rownames("Sample_ID") #%>% 
  #mutate(Th = ifelse(Mean_OutliersRemoved_Respiration_Rate_mg_DO_per_L_per_H > 3.5, "theoretical", "real"))
#limits = c(0.9, 4.35)
limits = c(-2, 2)

ggplot(cube_best_effect, aes(x = Final_Gravimetric_Water, y = Mean_OutliersRemoved_Respiration_Rate_mg_DO_per_kg_per_H)) + 
  #geom_point() +
  geom_point(aes(color = Mean_Fe_mg_per_kg)) +
  #geom_smooth() +
scale_color_gradient2(limits = limits, low = "firebrick2", mid = "goldenrod2", high = "dodgerblue2", midpoint = (max(limits)+min(limits))/2) +
  theme_bw()

## Scale cube effect size
z_cube_best_effect = cube_best_effect %>% 
  mutate(across(where(is.numeric), function(x) ((x - mean(x)) / sd(x))))

## EFFECT SIZE LASSO ####
## Set response variable
eff <- z_cube_best_effect$Mean_OutliersRemoved_Respiration_Rate_mg_DO_per_kg_per_H

## Set predictor variables
#Try all data
#z_cube_effect_pred <- data.matrix(z_cube_best_effect[, c('Mean_Fe_mg_per_L', 'Mean_ATP_nanomol_per_L', 'Percent_Fine_Sand', 'Percent_Coarse_Sand', 'Percent_Tot_Sand', 'Percent_Silt', 'Percent_Clay', 'mean_ssa', 'Mean_SpC', 'Mean_Temp', 'Mean_pH', 'Initial_Gravimetric_Water', 'Final_Gravimetric_Water')])

## Chosen variables: FS, SSA, Fe, Moi, ATP, SpC, Temp, pH
z_cube_effect_pred <- data.matrix(z_cube_best_effect[, c("Mean_Fe_mg_per_kg", "Mean_ATP_picomol_per_g",  "Percent_Fine_Sand", "mean_ssa", "Mean_SpC", "Mean_Temp", "Mean_pH",  "Final_Gravimetric_Water")])

## Start LASSO
#alpha = 1 is for LASSO regression
cv_model <- cv.glmnet(z_cube_effect_pred, eff, alpha = 1)

#Find lambda for lowest MSE using k-fold cross-validation
best_lambda <- cv_model$lambda.min
best_lambda

plot(cv_model)

# Rerun LASSO with best lambda and get coefficients
best_model <- glmnet(z_cube_effect_pred, eff, alpha = 1, lambda = best_lambda)
coef(best_model)

## Get R Sq. of best model
eff_predict <- predict(best_model, s = best_lambda, newx = z_cube_effect_pred)

sst <- sum((eff - mean(eff))^2)
sse <- sum((eff_predict - eff)^2)

rsq = 1 - sse/sst
rsq

# Get model residuals
residuals = eff - eff_predict

par(mfrow=c(2,2)) # Set up a 2x2 grid of plots
hist(residuals, main="Histogram of Residuals")
qqnorm(residuals, main="QQ Plot of Residuals")
qqline(residuals)
plot(density(residuals), main="Density Plot of Residuals")
plot(residuals ~ eff_predict)

dev.off()


## START PCA ####
z_pca = z_cube_theor_effect %>% 
  select(c(Mean_OutliersRemoved_Respiration_Rate_mg_DO_per_L_per_H,Mean_Fe_mg_per_L, Mean_ATP_nanomol_per_L, Percent_Fine_Sand, mean_ssa,Mean_SpC,Mean_Temp, Mean_pH, Final_Gravimetric_Water)) 

effect_pca <- prcomp(z_pca, scale = TRUE,                   center = TRUE, retx = T)

# Summary
summary(effect_pca)

# See the principal components
dim(effect_pca$x)
effect_pca$x

limits = c(0,
           60)

ind <- get_pca_ind(effect_pca)
ind

#png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/", as.character(Sys.Date()),"_Effect_PCA.png"), width = 10, height = 10, units = "in", res = 500)

fviz_pca_biplot(effect_pca, col.var = "black",geom = "point"
)+
  geom_point(aes(color = theor_effect$Mean_OutliersRemoved_Respiration_Rate_mg_DO_per_L_per_H), size = 3.5)+
  scale_color_gradient2(limits = limits, low = "firebrick2", mid = "goldenrod2",
                        high = "dodgerblue2", midpoint = (max(limits)+min(limits))/2) +
  labs(color = paste0("Wet - Dry Rate"))

dev.off()


#####
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

## LOG TRANSFORMATIONS ####
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
#z_log_pred <- data.matrix(z_log_lasso[, c('Fe_mg_per_L', 'Percent_Fine_Sand', 'Percent_Med_Sand', 'Percent_Coarse_Sand', 'Percent_Tot_Sand', 'Percent_Silt', 'Percent_Clay', 'mean_ssa', 'SpC', 'Temp', 'pH', 'Initial_Gravimetric_Water', 'Final_Gravimetric_Water')])

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
#z_log_pred_wet <- data.matrix(wet_log_lasso[, c('Fe_mg_per_L', 'Percent_Fine_Sand', 'Percent_Coarse_Sand', 'Percent_Tot_Sand', 'Percent_Silt', 'Percent_Clay', 'mean_ssa', 'SpC', 'Temp', 'pH', 'Initial_Gravimetric_Water', 'Final_Gravimetric_Water')])

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

## All Samples PCA ####

z_all_pca = cube_lasso %>% 
  select(c(Respiration_Rate_mg_DO_per_L_per_H,Fe_mg_per_L, ATP_nanomol_per_L, Percent_Fine_Sand, mean_ssa,SpC,Temp, pH, Final_Gravimetric_Water)) 

all_pca <- prcomp(z_all_pca, scale = TRUE,                   center = TRUE, retx = T)

# Summary
summary(all_pca)

# See the principal components
dim(all_pca$x)
all_pca$x

limits = c(0,
           300)

ind <- get_pca_ind(all_pca)
ind

#png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/", as.character(Sys.Date()),"_Effect_PCA.png"), width = 10, height = 10, units = "in", res = 500)

fviz_pca_biplot(all_pca, col.var = "black",geom = "point"
)+
  geom_point(aes(color = all_data_mg_L$Respiration_Rate_mg_DO_per_L_per_H, shape = all_data_mg_L$treat), size = 2)+
  scale_color_gradient2(limits = limits, low = "firebrick2", mid = "goldenrod2",
                        high = "dodgerblue2", midpoint = (max(limits)+min(limits))/2) +
  labs(color = paste0("Respiration Rate"))

dev.off()





