#### ECA LASSO ####
library(tidyverse)
library(glmnet)
library(factoextra)
library(janitor)

rm(list=ls());graphics.off()

## EFFECT SIZE DATA ####
effect_data <- read.csv("C:/Github/ECA_Multireactor_Incubations/Data/Cleaned Data/Effect_Median_ECA_Data.csv",header = TRUE) %>% 
  select(-c(X))

## Cube Effect Size
# Take out Fe outlier
cube_best_effect = effect_data %>% 
  mutate(across(where(is.numeric), cube_root)) %>% 
  rename_with(.cols = c(diff_median_SpC:median_mean_ssa), .fn = ~ paste0("cube_", .x)) %>% 
  select(-c(Rep, cube_diff_median_Respiration_Rate_mg_DO_per_L_per_H, cube_diff_median_Fe_mg_per_L, cube_diff_median_ATP_nanomol_per_L, cube_diff_median_Initial_Gravimetric_Moisture, cube_diff_median_Lost_Gravimetric_Water, cube_median_Percent_Med_Sand, cube_median_Percent_Coarse_Sand, cube_median_Percent_Tot_Sand, cube_median_Percent_Silt, cube_median_Percent_Clay)) %>% 
  column_to_rownames("Sample_ID") %>% 
  filter(cube_diff_median_Fe_mg_per_kg >-1)#%>% 

## Scale cube effect size before LASSO
z_cube_best_effect = cube_best_effect %>% 
  mutate(across(where(is.numeric), function(x) ((x - mean(x)) / sd(x)))) 

# mean(z_cube_best_effect$cube_diff_median_SpC)
# sd(z_cube_best_effect$cube_diff_median_SpC)

## EFFECT SIZE LASSO #############################

## Set response variable
yvar <- data.matrix(z_cube_best_effect$cube_diff_median_Respiration_Rate_mg_DO_per_kg_per_H)

## Set predictor variables
#Try with all data
#z_cube_effect_pred <- data.matrix(z_cube_best_effect[, c('Mean_Fe_mg_per_L', 'Mean_ATP_nanomol_per_L', 'Percent_Fine_Sand', 'Percent_Coarse_Sand', 'Percent_Tot_Sand', 'Percent_Silt', 'Percent_Clay', 'mean_ssa', 'Mean_SpC', 'Mean_Temp', 'Mean_pH', 'Initial_Gravimetric_Water', 'Final_Gravimetric_Water')])

## Try with chosen variables: FS, SSA, Fe, Moi, ATP, SpC, Temp, pH
xvar <- data.matrix(z_cube_best_effect[, c("cube_diff_median_Fe_mg_per_kg", "cube_diff_median_ATP_picomol_per_g",  "cube_median_Percent_Fine_Sand", "cube_median_mean_ssa", "cube_diff_median_SpC", "cube_diff_median_Temp", "cube_diff_median_pH",  "cube_diff_median_Final_Gravimetric_Moisture")])

## Start LASSO
#set.seed(42)

#alpha = 1 is for LASSO regression
cv_model <- cv.glmnet(xvar, yvar, alpha = 1, standardize = FALSE,  standardize.response = FALSE, intercept = FALSE)

#Find lambda for lowest MSE using k-fold cross-validation
best_lambda <- cv_model$lambda.min
best_lambda

plot(cv_model)

# Rerun LASSO with best lambda and get coefficients
best_model <- glmnet(xvar, yvar, alpha = 1, lambda = best_lambda, family = 'gaussian', standardize = FALSE,  standardize.response = FALSE, intercept = FALSE)

coef(best_model)

## Get R Sq. of best model
yvar_predict <- predict(best_model, s = best_lambda, newx = xvar)

sst <- sum((yvar - mean(yvar_predict))^2)
sse <- sum((yvar_predict - yvar)^2)

rsq = 1 - sse/sst
rsq

# Get model residuals
residuals = yvar - yvar_predict

par(mfrow=c(2,2)) # Set up a 2x2 grid of plots
hist(residuals, main="Histogram of Residuals")
qqnorm(residuals, main="QQ Plot of Residuals")
qqline(residuals)
plot(density(residuals), main="Density Plot of Residuals")
plot(residuals ~ eff_predict)

dev.off()


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

all_data = read.csv("C:/Github/ECA_Multireactor_Incubations/Data/Cleaned Data/All_ECA_Data_05-08-2024.csv")  %>% 
  filter(!grepl("052", Sample_Name)) %>% 
  filter(!grepl("053", Sample_Name)) %>% 
  filter(!grepl("057", Sample_Name)) %>% 
  filter(Respiration_Rate_mg_DO_per_L_per_H != -9999) %>% 
  left_join(icon_resp, by = "Sample_ID")

#remove site 023 because no moisture data
all_data_mg_L = all_data %>% 
  select(-c(Sample_ID, INC)) %>%
  filter(!grepl("023", Sample_Name)) %>% 
  mutate(Respiration_Rate_mg_DO_per_L_per_H = abs(Respiration_Rate_mg_DO_per_L_per_H)) %>% 
    mutate(Respiration_Rate_mg_DO_per_kg_per_H = abs(Respiration_Rate_mg_DO_per_kg_per_H)) %>% 
  rename(ICON_resp = Mean_Normalized_Respiration_Rate_mg_DO_per_H_per_L_sediment) %>% 
  mutate(Treat = ifelse(grepl("W", Sample_Name), "wet", "dry"))

### TRANSFORM DATA ####

cube_root <- function(x) sign(x) * (abs(x))^(1/3)

## CUBE TRANSFORMATIONS ####

## All mg/L values 
## theoretical vs. real: mg/L > 4.5

cube_lasso = all_data_mg_L %>% 
  mutate(across(where(is.numeric), cube_root)) %>%
  #filter(Respiration_Rate_mg_DO_per_L_per_H < 4.5) %>% 
  #filter(grepl("W", Sample_Name)) #%>% 
  #mutate(Treat = ifelse(grepl("W", Sample_Name), "Wet", "Dry")) %>% 
  column_to_rownames("Sample_Name") 

## Scale cube root values 
#z_cube_lasso = cube_lasso %>% 
  # mutate(across(where(is.numeric), function(x) ((x - mean(x)) / sd(x))))

## ALL DATA LASSO ####
## Set response variable
z_cube_resp <- cube_lasso$Respiration_Rate_mg_DO_per_kg_per_H

## Set predictor variables 

#Try all data
#z_cube_pred <- data.matrix(z_cube_lasso[, c('Fe_mg_per_L', 'Percent_Fine_Sand', 'Percent_Med_Sand', 'Percent_Coarse_Sand', 'Percent_Tot_Sand', 'Percent_Silt', 'Percent_Clay', 'mean_ssa', 'SpC', 'Temp', 'pH', 'Initial_Gravimetric_Water', 'Final_Gravimetric_Water')])

## Selected variables: FS, SSA, Fe, ICON Resp, Moi, ATP, SpC, Temp, pH

z_cube_pred <- data.matrix(cube_lasso[, c("Fe_mg_per_kg",  "ATP_picomol_per_g",  "Percent_Fine_Sand", "mean_ssa", "SpC", "Temp", "pH",  "Final_Gravimetric_Moisture", "ICON_resp")])

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

