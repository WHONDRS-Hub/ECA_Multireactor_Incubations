#### Sensitivity Analysis For ECA removals ####

# This script makes figures for ECA physical manuscript and performs PCA Analysis and LASSO regression

library(tidyverse)
library(corrplot)
library(ggpubr)
library(ggpmisc)
library(factoextra)
library(stringr)
library(glmnet)
library(magick)

rm(list=ls());graphics.off()

## Set image export

print.images = T

# Functions ---------------------------------------------------------------

# Transformation for normalization is cube root - have to cube root then add sign back to value to make it positive or negative
cube_root <- function(x) sign(x) * (abs(x))^(1/3)

# Read in/Merge Data ------------------------------------------------------------

## Individual Rate data for histograms ####
all_data = read.csv("C:/GitHub/ECA_Multireactor_Incubations/Data/EC_Sediment_SpC_pH_Temp_Respiration.csv", skip = 2) %>% 
  filter(grepl("EC", Sample_Name)) %>% 
  filter(!grepl("EC_011|EC_012|EC_023|EC_052|EC_053|EC_057", Sample_Name)) %>%  
  select(-c(Field_Name, IGSN, Material))

cube_respiration = all_data %>% 
  select(c(Sample_Name, Respiration_Rate_mg_DO_per_kg_per_H)) %>% 
  filter(Respiration_Rate_mg_DO_per_kg_per_H != -9999) %>% 
  mutate(cube_Respiration_mg_kg = cube_root(abs(as.numeric(Respiration_Rate_mg_DO_per_kg_per_H)))) %>% 
  mutate(Treat = if_else(grepl("D", Sample_Name), "Dry", "Wet"))

## Median Data ####

## Read in Median Data to get Dry moisture values

median = read.csv("C:/GitHub/ECA_Multireactor_Incubations/Data/EC_Sediment_Sample_Data_Summary.csv", skip = 2) %>% 
  filter(grepl("EC", Sample_Name)) %>% 
  separate(Sample_Name, c("Sample_Name", "Rep"), sep = "-") %>% 
  mutate(Sample_Name = paste0(Sample_Name, "_all")) %>% 
  filter(!grepl("EC_011|EC_012|EC_023|EC_052|EC_053|EC_057", Sample_Name)) %>%  
  select(-c(Field_Name, IGSN, Material, Median_Missing_Reps, Median_Respiration_Rate_mg_DO_per_kg_per_H)) %>% 
  select(-matches("per_L")) %>% 
  rename_with(~ str_remove_all(., "_[0-9]+")) %>% 
  rename_with(~ str_replace(., "^(([^_]*_){2}[^_]*).*", "\\1")) %>%
  rename_with(~ str_replace_all(., "Median", "median")) %>% 
  rename(median_SpC = median_SpC_microsiemens) %>% 
  rename(median_Temp = median_Temperature_degC)

## Effect Size Data ####

effect = read.csv("C:/GitHub/ECA_Multireactor_Incubations/Data/EC_Sediment_Effect_Size.csv", skip = 2) %>% 
  filter(grepl("EC", Sample_Name)) %>% 
  filter(!grepl("EC_011|EC_012|EC_023|EC_052|EC_053|EC_057", Sample_Name)) %>%
  mutate(Effect_Size = abs(as.numeric(Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H))) %>% 
  select(c(Sample_Name, Effect_Size)) 

## Read in grain size/ssa variables ####

grain = read.csv("C:/GitHub/ECA_Multireactor_Incubations/Data/v3_CM_SSS_Sediment_Sample_Data_Summary.csv", skip = 2) %>% 
  filter(grepl("CM", Sample_Name)) %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "CM", "EC")) %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "Sediment", "all")) %>% 
  select(c(Sample_Name, Percent_Tot_Sand, Percent_Coarse_Sand, Percent_Med_Sand, Percent_Fine_Sand, Percent_Silt, Percent_Clay, Mean_Specific_Surface_Area_m2_per_g))

## Join all data

effect_data = left_join(effect, grain, by = "Sample_Name") %>% 
  left_join(median, by = "Sample_Name") %>% 
  relocate(Rep, .after = "Sample_Name") %>% 
  mutate_at(vars(Percent_Tot_Sand:median_N_percent), as.numeric)  # make data numeric 

# Transform Data ----------------------------------------------------------

## Cube PCA for LASSO####
# Fe outlier not in analysis

cube_effect_data = effect_data %>% 
  mutate(across(where(is.numeric), cube_root)) %>% # cube root transform data
  rename_with(where(is.numeric), .fn = ~ paste0("cube_", .x)) #%>% 
  #column_to_rownames("Sample_Name") %>% 
  #filter(cube_Effect_Fe_mg_kg > -1) # remove Fe outlier for analysis

dry_cube_effect = effect_data %>% 
  filter(Rep == "D") %>% 
  column_to_rownames("Sample_Name") %>% 
  select(-c(Rep))

wet_cube_effect = effect_data %>% 
  filter(Rep == "W")%>% 
  column_to_rownames("Sample_Name") %>% 
  select(-c(Rep))

## Pearson Correlation Matrix ####

## Dry Medians 
# scale data before it goes into correlation matrix

scale_cube_dry_effect = as.data.frame(scale(dry_cube_effect))

scale_cube_dry_effect_pearson <- cor(scale_cube_dry_effect, method = "pearson")

if (print.images == T){
  
  png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Cube_Dry_Median_Effect_Pearson_Correlation_Matrix.png"), width = 12, height = 12, units = "in", res = 300)
  
  corrplot(scale_cube_dry_effect_pearson,type = "upper", method = "number", tl.col = "black", tl.cex = 1, cl.cex = 1,  title = "Effect Samples Pearson Correlation")
  
}

dev.off()

## Wet Medians
# scale data before it goes into correlation matrix

scale_cube_wet_effect = as.data.frame(scale(wet_cube_effect))

scale_cube_wet_effect_pearson <- cor(scale_cube_wet_effect, method = "pearson")

if (print.images == T){
  
  png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Cube_Wet_Median_Effect_Pearson_Correlation_Matrix.png"), width = 12, height = 12, units = "in", res = 300)
  
  corrplot(scale_cube_wet_effect_pearson,type = "upper", method = "number", tl.col = "black", tl.cex = 1, cl.cex = 1,  title = "Effect Samples Pearson Correlation")
  
}

dev.off()


## Downselected LASSO - Loop through coefficients to choose for LASSO ####

## Dry Medians 

# 1) Pivot data frame and sort highest to lowest

dry_pearson_df <- as.data.frame(scale_cube_dry_effect_pearson)

row_names_dry_pearson <- rownames(dry_pearson_df)

dry_pearson_df$Variable <- row_names_dry_pearson

# Melt the dataframe for plotting
dry_pearson_melted <- reshape2::melt(dry_pearson_df, id.vars = "Variable") %>% 
  filter(value != 1) %>% 
  mutate(value = abs(value)) %>% 
  filter(!grepl("Respiration", Variable))

dry_effect_melted <- dry_pearson_melted %>% 
  filter(grepl("Effect_Size", variable)) %>%
  filter(!grepl("Silt", Variable)) # remove silt from variables, lots of 0 values so not using

dry_choose_melted <- dry_pearson_melted %>% 
  filter(!grepl("Respiration", variable)) %>%
  filter(!grepl("Silt", variable)) %>%
  filter(!grepl("Silt", Variable)) %>% #try removing silt (0 values)
  #distinct(value, .keep_all = TRUE) %>% 
  left_join(dry_effect_melted, by = "Variable") %>% 
  rename(Variable_1 = Variable) %>% 
  rename(Variable_2 = variable.x) %>% 
  rename(Correlation = value.x) %>% 
  rename(Variable_1_Effect_Correlation = value.y) %>% 
  select(-c(variable.y)) %>% 
  left_join(dry_effect_melted, by = c("Variable_2" = "Variable")) %>% 
  rename(Variable_2_Effect_Correlation = value) %>% 
  select(-c(variable))

dry_loop_melt = dry_choose_melted %>% 
  arrange(desc(Correlation))

# Pearson correlation coefficient to remove above
correlation = 0.7

## Start loop to remove highly correlated (> 0.5)
dry_effect_filter = function(dry_loop_melt) {
  
  rows_to_keep = rep(TRUE, nrow(dry_loop_melt))
  
  for (i in seq_len(nrow(dry_loop_melt))) {
    
    if (!rows_to_keep[i]) next
    
    row = dry_loop_melt[i, ]
    
    if (row$Correlation < correlation) next
    
    if(row$Variable_1_Effect_Correlation >= row$Variable_2_Effect_Correlation) {
      
      var_to_keep = row$Variable_1
      var_to_remove = row$Variable_2
      
    } else {
      
      var_to_keep = row$Variable_2
      var_to_remove = row$Variable_1
      
    }
    
    dry_loop_melt$Variable_to_Keep[i] = var_to_keep
    dry_loop_melt$Variable_to_Remove[i] = var_to_remove
    
    for (j in seq(i + 1, nrow(dry_loop_melt))) {
      
      if(dry_loop_melt$Variable_1[j] == var_to_remove || dry_loop_melt$Variable_2[j] == var_to_remove) {
        
        rows_to_keep[j] = FALSE
        
      }
      
    }
    
    
  }
  
  return(dry_loop_melt[rows_to_keep, ])
  
}

dry_filtered_data = dry_effect_filter(dry_loop_melt) 

# pull out variables to remove
dry_removed_variables = dry_filtered_data %>% 
  distinct(Variable_to_Remove)

# pull out all variables 
dry_all_variables = dry_effect_melted %>% 
  select(c(Variable))

# remove variables from all variables to get variables to keep for LASSO 
dry_kept_variables = dry_effect_melted[!(dry_effect_melted$Variable %in% dry_removed_variables$Variable_to_Remove), ] #keeps SpC, Temp, pH, ATP, NPOC, TN (ext), TOC, TN (solid), med sand, silt

# if silt is removed, keeps SpC, Temp, pH, ATP, NPOC, TOC, TN (solid), med sand, fine sand, lost grav. moisture

## LASSO VARIABLES ####

# Keep variables selected from down-selected correlation matrix and add Cube_Effect_Size
dry_col_to_keep = unique(dry_kept_variables$Variable)
dry_col_to_keep = c(dry_col_to_keep, "Effect_Size")

dry_cube_variables = dry_cube_effect[, dry_col_to_keep, drop = FALSE]

## LASSO with Correlation Matrix Selected Variables ####
set.seed(42)
## Set response variable (Cube_Effect_Size) and scale
yvar <- data.matrix(scale(dry_cube_variables$Effect_Size, center = TRUE, scale = TRUE))
mean(yvar)
sd(yvar)

## Set predictor variables and scale
exclude_col = "Effect_Size"

dry_x_cube_variables = as.data.frame(scale(dry_cube_variables[, !(names(dry_cube_variables) %in% exclude_col)], center = T, scale = T))
#mean(x_cube_variables$Cube_SpC_Diff)
#sd(x_cube_variables$Cube_SpC_Diff)

xvars <- data.matrix(dry_x_cube_variables)

lasso = cv.glmnet(xvars, yvar, alpha = 1, nfolds = 5,
                  ,standardize = FALSE, standardize.response = FALSE, intercept = FALSE
                  #,standardize = TRUE, standardize.response = TRUE, intercept = FALSE
                  # , standardize = TRUE, standardize.response = FALSE, intercept = FALSE
)

best_lambda <- lasso$lambda.min
best_lambda

plot(lasso)

best_lasso_model <- glmnet(xvars, yvar, alpha = 1, lambda = best_lambda, family = "gaussian",
                           , standardize = FALSE, standardize.response = FALSE, intercept = FALSE
                           #  , standardize = TRUE, standardize.response = TRUE, intercept = FALSE
                           #, standardize = TRUE, standardize.response = FALSE, intercept = FALSE
)

#lasso_coefs = coef(best_lasso_model)
ds_lasso_coefs = coef(best_lasso_model)
#Fine Sand, ATP, SpC, N% Tot Sand
yvar_predict <- predict(best_lasso_model, s = best_lambda, newx = xvars)

sst <- sum((yvar - mean(yvar))^2)
sse <- sum((yvar_predict - yvar)^2)

rsq = 1 - sse/sst

rsq #0.455

## Wet Medians 

# 1) Pivot data frame and sort highest to lowest

wet_pearson_df <- as.data.frame(scale_cube_wet_effect_pearson)

row_names_wet_pearson <- rownames(wet_pearson_df)

wet_pearson_df$Variable <- row_names_wet_pearson

# Melt the dataframe for plotting
wet_pearson_melted <- reshape2::melt(wet_pearson_df, id.vars = "Variable") %>% 
  filter(value != 1) %>% 
  mutate(value = abs(value)) %>% 
  filter(!grepl("Respiration", Variable))

wet_effect_melted <- wet_pearson_melted %>% 
  filter(grepl("Effect_Size", variable)) %>%
  filter(!grepl("Silt", Variable)) # remove silt from variables, lots of 0 values so not using

wet_choose_melted <- wet_pearson_melted %>% 
  filter(!grepl("Respiration", variable)) %>%
  filter(!grepl("Silt", variable)) %>%
  filter(!grepl("Silt", Variable)) %>% #try removing silt (0 values)
  #distinct(value, .keep_all = TRUE) %>% 
  left_join(wet_effect_melted, by = "Variable") %>% 
  rename(Variable_1 = Variable) %>% 
  rename(Variable_2 = variable.x) %>% 
  rename(Correlation = value.x) %>% 
  rename(Variable_1_Effect_Correlation = value.y) %>% 
  select(-c(variable.y)) %>% 
  left_join(wet_effect_melted, by = c("Variable_2" = "Variable")) %>% 
  rename(Variable_2_Effect_Correlation = value) %>% 
  select(-c(variable))

wet_loop_melt = wet_choose_melted %>% 
  arrange(desc(Correlation))

# Pearson correlation coefficient to remove above
correlation = 0.7

## Start loop to remove highly correlated (> 0.5)
wet_effect_filter = function(wet_loop_melt) {
  
  rows_to_keep = rep(TRUE, nrow(wet_loop_melt))
  
  for (i in seq_len(nrow(wet_loop_melt))) {
    
    if (!rows_to_keep[i]) next
    
    row = wet_loop_melt[i, ]
    
    if (row$Correlation < correlation) next
    
    if(row$Variable_1_Effect_Correlation >= row$Variable_2_Effect_Correlation) {
      
      var_to_keep = row$Variable_1
      var_to_remove = row$Variable_2
      
    } else {
      
      var_to_keep = row$Variable_2
      var_to_remove = row$Variable_1
      
    }
    
    wet_loop_melt$Variable_to_Keep[i] = var_to_keep
    wet_loop_melt$Variable_to_Remove[i] = var_to_remove
    
    for (j in seq(i + 1, nrow(wet_loop_melt))) {
      
      if(wet_loop_melt$Variable_1[j] == var_to_remove || wet_loop_melt$Variable_2[j] == var_to_remove) {
        
        rows_to_keep[j] = FALSE
        
      }
      
    }
    
    
  }
  
  return(wet_loop_melt[rows_to_keep, ])
  
}

wet_filtered_data = wet_effect_filter(wet_loop_melt) 

# pull out variables to remove
wet_removed_variables = wet_filtered_data %>% 
  distinct(Variable_to_Remove)

# pull out all variables 
wet_all_variables = wet_effect_melted %>% 
  select(c(Variable))

# remove variables from all variables to get variables to keep for LASSO 
wet_kept_variables = wet_effect_melted[!(wet_effect_melted$Variable %in% wet_removed_variables$Variable_to_Remove), ] #keeps SpC, Temp, pH, ATP, NPOC, TN (ext), TOC, TN (solid), med sand, silt

# if silt is removed, keeps SpC, Temp, pH, ATP, NPOC, TOC, TN (solid), med sand, fine sand, lost grav. moisture

## LASSO VARIABLES ####

# Keep variables selected from down-selected correlation matrix and add Cube_Effect_Size
wet_col_to_keep = unique(wet_kept_variables$Variable)
wet_col_to_keep = c(wet_col_to_keep, "Effect_Size")

wet_cube_variables = wet_cube_effect[, wet_col_to_keep, drop = FALSE]

## LASSO with Correlation Matrix Selected Variables ####
set.seed(42)
## Set response variable (Cube_Effect_Size) and scale
yvar <- data.matrix(scale(wet_cube_variables$Effect_Size, center = TRUE, scale = TRUE))
mean(yvar)
sd(yvar)

## Set predictor variables and scale
exclude_col = "Effect_Size"

wet_x_cube_variables = as.data.frame(scale(wet_cube_variables[, !(names(wet_cube_variables) %in% exclude_col)], center = T, scale = T))
#mean(x_cube_variables$Cube_SpC_Diff)
#sd(x_cube_variables$Cube_SpC_Diff)

xvars <- data.matrix(wet_x_cube_variables)

lasso = cv.glmnet(xvars, yvar, alpha = 1, nfolds = 5,
                  ,standardize = FALSE, standardize.response = FALSE, intercept = FALSE
                  #,standardize = TRUE, standardize.response = TRUE, intercept = FALSE
                  # , standardize = TRUE, standardize.response = FALSE, intercept = FALSE
)

best_lambda <- lasso$lambda.min
best_lambda

plot(lasso)

best_lasso_model <- glmnet(xvars, yvar, alpha = 1, lambda = best_lambda, family = "gaussian",
                           , standardize = FALSE, standardize.response = FALSE, intercept = FALSE
                           #  , standardize = TRUE, standardize.response = TRUE, intercept = FALSE
                           #, standardize = TRUE, standardize.response = FALSE, intercept = FALSE
)

#lasso_coefs = coef(best_lasso_model)
ds_lasso_coefs = coef(best_lasso_model)
#Fine Sand, ATP, SpC, N% Tot Sand
yvar_predict <- predict(best_lasso_model, s = best_lambda, newx = xvars)

sst <- sum((yvar - mean(yvar))^2)
sse <- sum((yvar_predict - yvar)^2)

rsq = 1 - sse/sst

rsq #0.42

