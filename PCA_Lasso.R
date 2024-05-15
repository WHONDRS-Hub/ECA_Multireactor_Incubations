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

effect_data <- read.csv("C:/Github/ECA_Multireactor_Incubations/Data/Cleaned Data/Effect_Median_ECA_Data.csv",header = TRUE) %>% 
  select(-c(X)) 

## Functions ####

cube_root <- function(x) sign(x) * (abs(x))^(1/3)

## Cube PCA for LASSO####
# Fe outlier not in analysis
# Effect Size not in PCA
cube_effect_data = effect_data %>% 
  mutate(across(where(is.numeric), cube_root)) %>% 
  rename_with(where(is.numeric), .fn = ~ paste0("cube_", .x))

cube_effect_data_corr = cube_effect_data %>% 
  column_to_rownames("Sample_ID")%>% 
  select(-c(Rep, cube_diff_median_Respiration_Rate_mg_DO_per_L_per_H, cube_diff_median_Fe_mg_per_L, cube_diff_median_ATP_nanomol_per_L)) %>% 
  rename(Cube_SpC_Diff = cube_diff_median_SpC) %>% 
  rename(Cube_pH_Diff = cube_diff_median_pH) %>%
  rename(Cube_Temp_Diff = cube_diff_median_Temp) %>%
  rename(Cube_Effect_Size = cube_diff_median_Respiration_Rate_mg_DO_per_kg_per_H) %>%
  rename(Cube_Fe_mg_kg_Diff = cube_diff_median_Fe_mg_per_kg) %>%
  rename(Cube_InGravMoi_Diff = cube_diff_median_Initial_Gravimetric_Moisture) %>%
  rename(Cube_FinGravMoi_Diff = cube_diff_median_Final_Gravimetric_Moisture) %>%
  rename(Cube_LostGravMoi_Diff = cube_diff_median_Lost_Gravimetric_Water) %>%
  rename(Cube_ATP_pmol_g_Diff = cube_diff_median_ATP_picomol_per_g) %>%
  rename(Cube_Fine_Sand = cube_median_Percent_Fine_Sand) %>%
  rename(Cube_Med_Sand = cube_median_Percent_Med_Sand) %>%
  rename(Cube_Coarse_Sand = cube_median_Percent_Coarse_Sand) %>%
  rename(Cube_Tot_Sand = cube_median_Percent_Tot_Sand) %>%
  rename(Cube_Silt = cube_median_Percent_Silt) %>%
  rename(Cube_Clay = cube_median_Percent_Clay) %>%
  rename(Cube_SSA = cube_median_mean_ssa) %>% 
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

# Plotting the heatmap
ggplot(loadings_melted, aes(x = variable, y = Variable, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(x = "Principal Component", y = "Variable", fill = "Loadings") +
  ggtitle("PCA Loadings Heatmap")+
  geom_text(aes(label = round(value, 2)), color = "black", size = 3)


## LASSO with PCA Loadings ####
set.seed(42)
## Set response variable (Cube_Effect_Size)
#yvar <- data.matrix(scale(cube_effect_data_clean$Cube_Effect_Size, center = TRUE, scale = TRUE))
yvar <- data.matrix(cube_effect_data_corr$Cube_Effect_Size)
mean(yvar)
sd(yvar)

## Set predictor variables
#pca_scores = as.data.frame(scale(cube_effect_pca$x, center = TRUE, scale = TRUE))
pca_scores = as.data.frame(cube_effect_pca$x)
#mean(pca_scores$PC1)
#sd(pca_scores$PC1)

xvars <- data.matrix(pca_scores[, c("PC1",  "PC2", "PC3", "PC4", "PC5")])

lasso = cv.glmnet(xvars, yvar, alpha = 1, nfolds = 5,
                  #,standardize = FALSE, standardize.response = FALSE, intercept = FALSE
                  #,standardize = TRUE, standardize.response = TRUE, intercept = FALSE
                  , standardize = TRUE, standardize.response = FALSE, intercept = FALSE
                  )

best_lambda <- lasso$lambda.min
best_lambda

plot(lasso)

best_lasso_model <- glmnet(xvars, yvar, alpha = 1, lambda = best_lambda, family = "gaussian",
#   , standardize = FALSE, standardize.response = FALSE, intercept = FALSE
#  , standardize = TRUE, standardize.response = TRUE, intercept = FALSE
  , standardize = TRUE, standardize.response = FALSE, intercept = FALSE
)

coef(best_lasso_model)

yvar_predict <- predict(best_lasso_model, s = best_lambda, newx = xvars)

sst <- sum((yvar - mean(yvar))^2)
sse <- sum((yvar_predict - yvar)^2)

rsq = 1 - sse/sst

rsq


cube_limits = c(-12,
                12)

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Cube_Median_Effect_PCA.png"), width = 8, height = 8, units = "in", res = 300)

cube_pca = fviz_pca_biplot(cube_effect_pca, col.var = "black",geom = "point"
)+
  geom_point(aes(color = cube_effect_data_clean$cube_Effect_Size), size = 3.5)+
  scale_color_gradient2(limits = cube_limits, low = "firebrick2", mid = "goldenrod2",
                        high = "dodgerblue2", midpoint = (max(cube_limits)+min(cube_limits))/2) +
  labs(color = paste0("cube_Wet - Dry Rate"))

cube_pca

dev.off()
