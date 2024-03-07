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

### LASSO REGRESSION
#Log transform all data

half_min_values = function(x) {
  min_val <- min(x)
  if (min_val == 0) {
    min_val <- sort(x[x != 0])[1]
  }
  return(0.5 * min_val)
}

cube_root <- function(x) sign(x) * (abs(x))^(1/3)

log_lasso = all_data_mg_L %>% 
  select(-c(Lost_Gravimetric_Water)) %>% 
  mutate(across(where(is.numeric), ~. + half_min_values(.))) %>% 
  mutate(across(where(is.numeric), ~log10(.))) %>% 
  column_to_rownames("Sample_Name")

cube_lasso = all_data_mg_L %>% 
  mutate(across(where(is.numeric), cube_root))%>% 
  column_to_rownames("Sample_Name")

cube_rows = cube_lasso %>% 
  rownames_to_column("Sample_Name")

cube_pivot = left_join(cube_rows, all_data_mg_L, by = "Sample_Name")

ggplot(cube_pivot, aes(y = Lost_Gravimetric_Water.x, x = Lost_Gravimetric_Water.y)) + 
  geom_point() + 
  ylab("cube")

## Plot all log transformations vs. non log transformed

ggplot(gather(all_data_mg_L, cols, value), aes(x = value)) + 
  geom_histogram() + 
  facet_wrap(.~cols, scales = 'free_x')

ggplot(gather(log_lasso, cols, value), aes(x = value)) + 
  geom_histogram() + 
  facet_wrap(.~cols, scales = 'free_x')

ggplot(gather(cube_lasso, cols, value), aes(x = value)) + 
  geom_histogram() + 
  facet_wrap(.~cols, scales = 'free_x')



resp <- log_lasso$Respiration_Rate_mg_DO_per_L_per_H

pred <- data.matrix(log_lasso[, c('Fe_mg_per_L', 'Percent_Fine_Sand', 'Percent_Coarse_Sand', 'Percent_Tot_Sand', 'Percent_Silt', 'Percent_Clay', 'mean_ssa', 'SpC', 'Temp', 'pH', 'Initial_Gravimetric_Water', 'Final_Gravimetric_Water')])

cv_model <- cv.glmnet(pred, resp, alpha = 1)

best_lambda <- cv_model$lambda.min
best_lambda

plot(cv_model)

best_model <- glmnet(pred, resp, alpha = 1, lambda = best_lambda)
coef(best_model)

resp_predict <- predict(best_model, s = best_lambda, newx = pred)

sst <- sum((resp - mean(resp))^2)
sse <- sum((resp_predict - resp)^2)

rsq = 1 - sse/sst

rsq

## LASSO ####

## Log LASSO

# lasso <- all_data %>% 
#   drop_na(Fe_mg_per_L) %>% 
#   drop_na(average_ssa) %>% 
#   drop_na(Initial_Gravimetric_Water) %>% 
#   mutate(Respiration_Rate_mg_DO_per_L_per_H = abs(Respiration_Rate_mg_DO_per_L_per_H))

log_lasso <- log_corr_samples %>% 
  rename(log_fe = `Log Fe (mg/L)`) %>% 
  rename(log_atp = `Log ATP (pmol/g)`) %>% 
  rename(log_fine_sand = `Log % Fine Sand`) %>% 
  rename(log_med_sand = `Log % Med. Sand`) %>% 
  rename(log_coarse_sand = `Log % Coarse Sand`) %>% 
  rename(log_clay = `Log % Clay`) %>% 
  rename(log_silt = `Log % Silt`) %>% 
  rename(log_ssa = `Log SSA`) %>% 
  rename(log_spc = `Log SpC`) %>% 
  rename(log_temp = `Log Temp`) %>% 
  rename(log_ph = `Log pH`) %>% 
  rename(log_in_grav_moi = `Log In.Grav.Moi.`) %>% 
  rename(log_fin_grav_moi = `Log Fin.Grav.Moi.`)

log_resp <- log_lasso$`Log Respiration (mg/L)`

#log_resp <- scale(log_resp)

log_pred <- data.matrix(log_lasso[, c("log_fe", "log_atp",  "log_fine_sand", "log_med_sand", "log_coarse_sand", "log_clay","log_silt", "log_ssa", "log_spc", "log_temp", "log_ph","log_in_grav_moi", "log_fin_grav_moi")])

#log_pred <- scale(log_pred)

log_cv_model <- cv.glmnet(log_pred, log_resp, alpha = 1)

log_best_lambda <- log_cv_model$lambda.min
log_best_lambda

plot(log_cv_model)

log_best_model <- glmnet(log_pred, log_resp, alpha = 1, lambda = log_best_lambda)
coef(log_best_model)

log_resp_predict <- predict(log_best_model, s = log_best_lambda, newx = log_pred)

log_sst <- sum((log_resp - mean(log_resp))^2)
log_sse <- sum((log_resp_predict - log_resp)^2)

log_rsq = 1 - log_sse/log_sst

log_rsq

## Log FS, SSA, Fe, Resp, Moi

log_pred_co <- data.matrix(log_lasso[, c("log_fe", "log_atp",  "log_fine_sand", "log_ssa", "log_spc", "log_temp", "log_ph", "log_fin_grav_moi")])

#log_pred <- scale(log_pred)

log_cv_model_co <- cv.glmnet(log_pred_co, log_resp, alpha = 1)

log_best_lambda_co <- log_cv_model_co$lambda.min
log_best_lambda_co

plot(log_cv_model_co)

log_best_model_co <- glmnet(log_pred_co, log_resp, alpha = 1, lambda = log_best_lambda_co)
coef(log_best_model_co)

log_resp_predict_co <- predict(log_best_model_co, s = log_best_lambda_co, newx = log_pred_co)

log_sst_co <- sum((log_resp - mean(log_resp))^2)
log_sse_co <- sum((log_resp_predict_co - log_resp)^2)

log_rsq_co = 1 - log_sse_co/log_sst_co

log_rsq_co

## Log LASSO Dry Samples

log_lasso_dry <- log_lasso %>% 
  rownames_to_column("Sample_Name") %>% 
  filter(grepl("D", Sample_Name)) %>%
  dplyr::select(c(`Log Respiration (mg/L)`, log_fe, log_atp, log_fine_sand, log_ssa, log_spc, log_temp, log_ph, log_fin_grav_moi))

log_resp_dry <- log_lasso_dry$`Log Respiration (mg/L)`

log_pred_dry <- data.matrix(log_lasso_dry[, c('log_fe', 'log_atp', 'log_fine_sand', 'log_ssa', 'log_spc', 'log_temp', 'log_ph', 'log_fin_grav_moi')])

#log_pred <- scale(log_pred)

log_cv_model_dry <- cv.glmnet(log_pred_dry, log_resp_dry, alpha = 1)

log_best_lambda_dry <- log_cv_model_dry$lambda.min
log_best_lambda_dry

plot(log_cv_model_dry)

log_best_model_dry <- glmnet(log_pred_dry, log_resp_dry, alpha = 1, lambda = log_best_lambda_dry)
coef(log_best_model_dry)

log_resp_predict_dry <- predict(log_best_model_dry, s = log_best_lambda_dry, newx = log_pred_dry)

log_sst_dry <- sum((log_resp_dry - mean(log_resp_dry))^2)
log_sse_dry <- sum((log_resp_predict_dry - log_resp_dry)^2)

log_rsq_dry = 1 - log_sse_dry/log_sst_dry

log_rsq_dry

## Log Wet LASSO
log_lasso_wet <- log_lasso %>% 
  rownames_to_column("Sample_Name") %>% 
  filter(grepl("W", Sample_Name)) %>%
  dplyr::select(c(`Log Respiration (mg/L)`, log_fe, log_atp, log_fine_sand, log_ssa, log_spc, log_temp, log_ph, log_fin_grav_moi))

log_resp_wet <- log_lasso_wet$`Log Respiration (mg/L)`

log_pred_wet <- data.matrix(log_lasso_wet[, c('log_fe', 'log_atp', 'log_fine_sand', 'log_ssa', 'log_spc', 'log_temp', 'log_ph', 'log_fin_grav_moi')])

#log_pred <- scale(log_pred)

log_cv_model_wet <- cv.glmnet(log_pred_wet, log_resp_wet, alpha = 1)

log_best_lambda_wet <- log_cv_model_wet$lambda.min
log_best_lambda_wet

plot(log_cv_model_wet)

log_best_model_wet <- glmnet(log_pred_wet, log_resp_wet, alpha = 1, lambda = log_best_lambda_wet)
coef(log_best_model_wet)

log_resp_predict_wet <- predict(log_best_model_wet, s = log_best_lambda_wet, newx = log_pred_wet)

log_sst_wet <- sum((log_resp_wet - mean(log_resp_wet))^2)
log_sse_wet <- sum((log_resp_predict_wet - log_resp_wet)^2)

log_rsq_wet = 1 - log_sse_wet/log_sst_wet

log_rsq_wet


## Effect LASSO 
effect_lasso <- log_effect

eff <- effect_lasso$`log_Effect Size`

#eff <- scale(eff)

pred <- data.matrix(effect_lasso[, c("log_Fe_mg_per_L_diff", "log_ATP_picomol_per_g_mean", "log_SpC Diff.", "log_Temp Diff.", "log_pH Diff.", "log_Fin.Grav.Moi. Diff.", "log_Percent_Fine_Sand_mean", "log_Percent_Med_Sand_mean", "log_Percent_Coarse_Sand_mean", "log_Percent_Clay_mean", "log_% Silt", "log_average_ssa_mean")]) 

#pred <- scale(pred)

cv_model <- cv.glmnet(pred, eff, alpha = 1)

best_lambda <- cv_model$lambda.min
best_lambda

plot(cv_model)

best_model <- glmnet(pred, eff, alpha = 1, lambda = best_lambda)
coef(best_model)

eff_predict <- predict(best_model, s = best_lambda, newx = pred)

sst <- sum((eff - mean(eff))^2)
sse <- sum((eff_predict - eff)^2)

rsq = 1 - sse/sst

rsq

eff <- effect_lasso$`log_Effect Size`

#eff <- scale(eff)

pred <- data.matrix(effect_lasso[, c("log_Fe_mg_per_L_diff", "log_ATP_picomol_per_g_mean", "log_SpC Diff.", "log_Temp Diff.", "log_pH Diff.", "log_Fin.Grav.Moi. Diff.", "log_Percent_Fine_Sand_mean", "log_average_ssa_mean")]) 

#pred <- scale(pred)

cv_model <- cv.glmnet(pred, eff, alpha = 1)

best_lambda <- cv_model$lambda.min
best_lambda

plot(cv_model)

best_model <- glmnet(pred, eff, alpha = 1, lambda = best_lambda)
coef(best_model)

eff_predict <- predict(best_model, s = best_lambda, newx = pred)

sst <- sum((eff - mean(eff))^2)
sse <- sum((eff_predict - eff)^2)

rsq = 1 - sse/sst

rsq


corr_samples_treat_real <- corr_samples %>%
  rownames_to_column("Sample_Name") %>% 
  mutate(Treat = case_when(grepl("W",Sample_Name)~"Wet",
                           grepl("D", Sample_Name) ~"Dry")) %>% 
  filter(`Respiration (mg/L)` < 90)

ggplot(corr_samples_treat_real, aes(x = `Respiration (mg/L)`, y = ATP_picomol_per_g, color = Treat)) + 
  geom_point()

ggplot(corr_samples_treat_real, aes(x = `Respiration (mg/L)`, y = Fe_mg_per_L, color = Treat)) + 
  geom_point()

ggplot(corr_samples_treat, aes(x = `Respiration (mg/L)`, y = Fe_mg_per_L, color = Treat)) + 
  geom_point()

dry_samples = corr_samples_treat %>% 
  filter()

ggplot(effect, aes(x = `Effect Size`, y = ATP_picomol_per_g_mean)) + 
  geom_point() + 
  ylab("ATP Diff.")

ggplot(effect, aes(x = `Effect Size`, y = Percent_Fine_Sand_mean)) + 
  geom_point() + 
  ylab("% Fine Sand")

ggplot(effect, aes(x = `Effect Size`, y = Fe_mg_per_L_diff)) + 
  geom_point() + 
  ylab("Fe Diff.")

ggplot(effect, aes(x = `Effect Size`, y = `Fin.Grav.Moi. Diff.`)) + 
  geom_point() + 
  ylab("Final Grav Moi. Diff.")
