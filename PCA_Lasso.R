#### Sensitivity Analysis For ECA removals ####
library(tidyverse)
library(corrplot)
library(ggpubr)
library(ggpmisc)
library(factoextra)
library(stringr)
library(glmnet)

rm(list=ls());graphics.off()

#### Read in Effect Size Data

#change this to published data when ready

effect_data <- read.csv("C:/Github/ECA_Multireactor_Incubations/Data/Cleaned Data/2024-05-29_Effect_Median_ECA_Data.csv",header = TRUE) %>% 
  dplyr::select(-c(X)) 

## Functions ####

cube_root <- function(x) sign(x) * (abs(x))^(1/3)

## Cube PCA for LASSO####
# Fe outlier not in analysis
# Effect Size not in PCA
cube_effect_data = effect_data %>% 
  mutate(across(where(is.numeric), cube_root)) %>% 
  rename_with(where(is.numeric), .fn = ~ paste0("cube_", .x))

cube_effect_data_corr = cube_effect_data %>% 
  column_to_rownames("Sample_ID") %>% 
  select(-c(Rep, cube_diff_median_Respiration_Rate_mg_DO_per_L_per_H, cube_diff_median_Fe_mg_per_L, cube_diff_median_ATP_nanomol_per_L, cube_diff_median_TN_mg_N_per_L, cube_diff_median_NPOC_mg_C_per_L)) %>% 
  rename(Cube_SpC_Diff = cube_diff_median_SpC) %>% 
  rename(Cube_pH_Diff = cube_diff_median_pH) %>%
  rename(Cube_Temp_Diff = cube_diff_median_Temp) %>%
  rename(Cube_Effect_Size = cube_diff_median_Respiration_Rate_mg_DO_per_kg_per_H) %>%
  rename(Cube_Fe_mg_kg_Diff = cube_diff_median_Fe_mg_per_kg) %>%
  rename(Cube_Dry_InGravMoi = cube_median_Dry_Initial_Gravimetric_Moisture) %>%
  rename(Cube_Dry_FinGravMoi = cube_median_Dry_Final_Gravimetric_Moisture) %>%
  rename(Cube_Dry_LostGravMoi = cube_median_Dry_Lost_Gravimetric_Moisture) %>%
  rename(Cube_ATP_pmol_g_Diff = cube_diff_median_ATP_picomol_per_g) %>%
  rename(Cube_NPOC_mg_C_per_kg_Diff = cube_diff_median_NPOC_mg_C_per_kg) %>% 
  rename(Cube_TN_mg_N_per_kg_Diff = cube_diff_median_TN_mg_N_per_kg) %>% 
  rename(Cube_Fine_Sand = cube_median_Percent_Fine_Sand) %>%
  rename(Cube_Med_Sand = cube_median_Percent_Med_Sand) %>%
  rename(Cube_Coarse_Sand = cube_median_Percent_Coarse_Sand) %>%
  rename(Cube_Tot_Sand = cube_median_Percent_Tot_Sand) %>%
  rename(Cube_Silt = cube_median_Percent_Silt) %>%
  rename(Cube_Clay = cube_median_Percent_Clay) %>%
  rename(Cube_SSA = cube_median_mean_ssa) %>% 
  rename(Cube_TN_Percent_Diff = cube_diff_median_tn_percent) %>% 
  rename(Cube_TOC_Percent_Diff = cube_diff_median_toc_percent) %>% 
  relocate(Cube_Effect_Size, .before = Cube_SpC_Diff)  %>% 
  filter(Cube_Fe_mg_kg_Diff > -1)

## Run PCA without effect size ####
cube_data_pca = cube_effect_data_corr %>% 
  select(-c(Cube_Effect_Size))

cube_effect_pca <- prcomp(cube_data_pca, scale = TRUE,
                          center = TRUE, retx = T)

summary(cube_effect_pca)

# Plot cumulative variation

cumulative_var = cumsum(cube_effect_pca$sdev^2)/sum(cube_effect_pca$sdev^2)

pca_plot = data.frame(PC = 1:length(cube_effect_pca$sdev), Cumulative_Variance = cumulative_var)

ggplot(pca_plot, aes(x = PC, y = Cumulative_Variance)) +
  geom_col()


# Look at PCA loadings

loadings = cube_effect_pca$rotation

print(loadings)

loadings_df <- as.data.frame(loadings) %>% 
  select(c(PC1, PC2, PC3, PC4, PC5))

# Generate heatmap

row_names <- rownames(loadings_df)
loadings_df$Variable <- row_names

# Melt the dataframe for plotting
loadings_melted <- reshape2::melt(loadings_df, id.vars = "Variable")

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Cubed_Root_PCA_Heatmap.png"), width = 12, height = 8, units = "in", res = 300)

# Plotting the heatmap
ggplot(loadings_melted, aes(x = variable, y = Variable, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(x = "Principal Component", y = "Variable", fill = "Loadings") +
  ggtitle("PCA Loadings Heatmap")+
  geom_text(aes(label = round(value, 2)), color = "black", size = 3)

dev.off()

## LASSO with PCA Loadings ####
set.seed(42)
## Set response variable (Cube_Effect_Size)
yvar <- data.matrix(scale(cube_effect_data_corr$Cube_Effect_Size, center = TRUE, scale = TRUE))
#yvar <- data.matrix(cube_effect_data_corr$Cube_Effect_Size)
mean(yvar)
sd(yvar)

## Set predictor variables
pca_scores = as.data.frame(scale(cube_effect_pca$x, center = TRUE, scale = TRUE))
#pca_scores = as.data.frame(cube_effect_pca$x)
mean(pca_scores$PC1)
sd(pca_scores$PC1)

xvars <- data.matrix(pca_scores[, c("PC1",  "PC2", "PC3", "PC4", "PC5")])

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

lasso_coefs = coef(best_lasso_model)

yvar_predict <- predict(best_lasso_model, s = best_lambda, newx = xvars)

sst <- sum((yvar - mean(yvar))^2)
sse <- sum((yvar_predict - yvar)^2)

rsq = 1 - sse/sst

rsq

lasso_df = as.data.frame(as.matrix(lasso_coefs))

colnames(lasso_df) = c("Coefficients")

lasso_df = lasso_df %>% 
  rownames_to_column(var = "variable") %>% 
  slice(-1)


weighted_pc = merge(lasso_df, loadings_melted, by = "variable") %>% 
  mutate(pc_weight = value * Coefficients)

weighted_total = weighted_pc %>% 
  group_by(Variable) %>% 
  mutate(total_weight = sum(pc_weight))

ggplot(weighted_total, aes(x = fct_reorder(Variable, total_weight))) +
  geom_col(aes(y = pc_weight, fill = variable)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggplot(weighted_total, aes(x = fct_reorder(Variable, total_weight), y = total_weight)) + 
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

abs_total = weighted_total %>% 
  mutate(across(where(is.numeric), abs))

ggplot(abs_total, aes(x = fct_reorder(Variable, total_weight), y = total_weight)) + 
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


cube_limits = c(-12,
                12)

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Cube_Median_Effect_PCA.png"), width = 8, height = 8, units = "in", res = 300)

cube_pca = fviz_pca_biplot(cube_effect_pca, col.var = "black",geom = "point"
) +
  geom_point(aes(color = cube_effect_data_corr$Cube_Effect_Size), size = 3.5) +
  scale_color_gradient2(limits = cube_limits, low = "firebrick2", mid = "goldenrod2",
                        high = "dodgerblue2", midpoint = (max(cube_limits)+min(cube_limits))/2) +
  labs(color = paste0("Cubed Root Wet - Dry Rate"))+ theme(legend.position = c(0.2, 0.85), 
                                                           legend.key.size = unit(0.25, "in"), 
                                                           legend.title = element_text(size = 8),
                                                           axis.title.x = element_text(size = 10))


cube_pca

dev.off()

## Pearson Correlation Matrix ####
scale_cube_effect_corr = as.data.frame(scale(cube_effect_data_corr))

cube_effect_samples_corr_pearson <- cor(scale_cube_effect_corr, method = "pearson")

# png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Cube_Median_Effect_Pearson_Correlation_Matrix.png"), width = 12, height = 12, units = "in", res = 300)

corrplot(cube_effect_samples_corr_pearson,type = "upper", method = "number", tl.col = "black", tl.cex = 1.6, cl.cex = 1.25,  title = "Effect Samples Pearson Correlation")
 
#dev.off()

# make one line correlation matrix with just effect size

corr_effect = matrix(cube_effect_samples_corr_pearson[1, ], nrow = 1)


colnames(corr_effect) = colnames(cube_effect_samples_corr_pearson)

rownames(corr_effect) = rownames(cube_effect_samples_corr_pearson)[1]

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Cube_Median_Effect_Pearson_Correlation_Matrix_One_Line.png"), width = 12, height = 12, units = "in", res = 300)
# 

corrplot(corr_effect, type = "upper", method = "number", tl.col = "black", tl.cex = 1.6, cl.cex = 1,  title = "Effect Samples Pearson Correlation", diag = FALSE, is.corr = FALSE, cl.pos = 'n')
# 
dev.off()

## Loop through coefficients to choose for LASSO ####

# 1) Pivot data frame and sort highest to lowest

pearson_df <- as.data.frame(cube_effect_samples_corr_pearson)

row_names_pearson <- rownames(pearson_df)

pearson_df$Variable <- row_names_pearson

# Melt the dataframe for plotting
pearson_melted <- reshape2::melt(pearson_df, id.vars = "Variable") %>% 
  filter(value != 1) %>% 
  mutate(value = abs(value)) %>% 
  filter(!grepl("Effect", Variable))

effect_melted <- pearson_melted %>% 
  filter(grepl("Effect", variable)) %>% 
  filter(!grepl("Silt", Variable))

choose_melted <- pearson_melted %>% 
  filter(!grepl("Effect", variable)) %>%
  filter(!grepl("Silt", variable)) %>% 
  filter(!grepl("Silt", Variable)) %>% #try removing silt (0 values)
  #distinct(value, .keep_all = TRUE) %>% 
  left_join(effect_melted, by = "Variable") %>% 
  rename(Variable_1 = Variable) %>% 
  rename(Variable_2 = variable.x) %>% 
  rename(Correlation = value.x) %>% 
  rename(Variable_1_Effect_Correlation = value.y) %>% 
  select(-c(variable.y)) %>% 
  left_join(effect_melted, by = c("Variable_2" = "Variable")) %>% 
  rename(Variable_2_Effect_Correlation = value) %>% 
  select(-c(variable))

loop_melt = choose_melted %>% 
  arrange(desc(Correlation))

## Start loop to remove highly correlated (> 0.5)
effect_filter = function(loop_melt) {
  
  rows_to_keep = rep(TRUE, nrow(loop_melt))
  
  for (i in seq_len(nrow(loop_melt))) {
    
    if (!rows_to_keep[i]) next
    
    row = loop_melt[i, ]
    
    if (row$Correlation < 0.7) next
    
    if(row$Variable_1_Effect_Correlation >= row$Variable_2_Effect_Correlation) {
      
      var_to_keep = row$Variable_1
      var_to_remove = row$Variable_2
      
    } else {
      
      var_to_keep = row$Variable_2
      var_to_remove = row$Variable_1
      
    }
    
    loop_melt$Variable_to_Keep[i] = var_to_keep
    loop_melt$Variable_to_Remove[i] = var_to_remove
   
    for (j in seq(i + 1, nrow(loop_melt))) {
      
      if(loop_melt$Variable_1[j] == var_to_remove || loop_melt$Variable_2[j] == var_to_remove) {
        
        rows_to_keep[j] = FALSE
        
      }
      
    }
     
    
  }
  
  return(loop_melt[rows_to_keep, ])
  
}
  
filtered_data = effect_filter(loop_melt) 

# pull out variables to remove
removed_variables = filtered_data %>% 
  distinct(Variable_to_Remove)

# pull out all variables 
all_variables = effect_melted %>% 
  select(c(Variable))

# remove variables from all variables to get variables to keep for LASSO 
kept_variables = effect_melted[!(effect_melted$Variable %in% removed_variables$Variable_to_Remove), ] #keeps SpC, Temp, pH, ATP, NPOC, TN (ext), TOC, TN (solid), med sand, silt

# if silt is removed, keeps SpC, Temp, pH, ATP, NPOC, TOC, TN (solid), med sand, fine sand, lost grav. moisture

## LASSO VARIABLES ####

# Keep variables selected from down-selected correlation matrix and add Cube_Effect_Size
col_to_keep = unique(kept_variables$Variable)
col_to_keep = c(col_to_keep, "Cube_Effect_Size")

cube_variables = cube_effect_data_corr[, col_to_keep, drop = FALSE]

## LASSO with Correlation Matrix Selected Variables ####
set.seed(42)
## Set response variable (Cube_Effect_Size) and scale
yvar <- data.matrix(scale(cube_variables$Cube_Effect_Size, center = TRUE, scale = TRUE))
mean(yvar)
sd(yvar)

## Set predictor variables and scale
exclude_col = "Cube_Effect_Size"

x_cube_variables = as.data.frame(scale(cube_variables[, !(names(cube_variables) %in% exclude_col)], center = T, scale = T))
#mean(x_cube_variables$Cube_SpC_Diff)
#sd(x_cube_variables$Cube_SpC_Diff)

xvars <- data.matrix(x_cube_variables)

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

lasso_coefs = coef(best_lasso_model)

yvar_predict <- predict(best_lasso_model, s = best_lambda, newx = xvars)

sst <- sum((yvar - mean(yvar))^2)
sse <- sum((yvar_predict - yvar)^2)

rsq = 1 - sse/sst

rsq

## LASSO with all variables to check for collinearity effects ####

## Set response variable (Cube_Effect_Size) and scale
yvar <- data.matrix(scale(cube_effect_data_corr$Cube_Effect_Size, center = TRUE, scale = TRUE))
mean(yvar)
sd(yvar)

## Set predictor variables and scale
x_cube_variables =  as.data.frame(scale(cube_effect_data_corr[, !(names(cube_effect_data_corr) %in% exclude_col)], center = T, scale = T))
#mean(x_cube_variables$Cube_SpC_Diff)
#sd(x_cube_variables$Cube_SpC_Diff)

xvars <- data.matrix(x_cube_variables)

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

lasso_coefs = coef(best_lasso_model)

yvar_predict <- predict(best_lasso_model, s = best_lambda, newx = xvars)

sst <- sum((yvar - mean(yvar))^2)
sse <- sum((yvar_predict - yvar)^2)

rsq = 1 - sse/sst

rsq
