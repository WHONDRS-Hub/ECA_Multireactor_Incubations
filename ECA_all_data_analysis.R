#### Sensitivity Analysis For ECA removals ####
library(tidyverse)

rm(list=ls());graphics.off()

#### Read in Data

#Individual samples 
all_data <- read.csv("C:/Github/ECA_Multireactor_Incubations/Data/Cleaned Data/All_ECA_Data_05-08-2024.csv",header = TRUE) 

#Summary Data 

sum_data <- read.csv("C:/Github/ECA_Multireactor_Incubations/Data/Cleaned Data/Medians_ECA_Data.csv",header = TRUE) %>% 
  select(-c(X))

#Effect Size Data
effect_data <- read.csv("C:/Github/ECA_Multireactor_Incubations/Data/Cleaned Data/Effect_Median_ECA_Data.csv",header = TRUE) %>% 
  select(-c(X))

## Functions ####

cube_root <- function(x) sign(x) * (abs(x))^(1/3)

panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r = (cor(x, y, method = c("spearman")))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) {cex.cor <- 0.8/strwidth(txt)} else {cex = cex.cor}
  text(0.5, 0.5, txt, cex = cex.cor * (1 + r)/1)
  
  # if(missing(cex.cor)) {cex <- 1.2/strwidth(txt)} else {cex = cex.cor}
  # text(0.5, 0.5, txt, cex = cex * sin(sqrt(abs(r))))
  
  test <- cor.test(x,y)
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " "))
  #text(0.5, 0.5, txt, cex = cex * r)
  text(.5, .8, Signif, cex=cex, col=2)
  
}


#### Rates histograms ####

clean_all_data = all_data %>%
  filter(Respiration_Rate_mg_DO_per_kg_per_H != 9999) %>% 
  mutate(Treat = case_when(grepl("W",INC)~"Wet",
                           grepl("D", INC) ~"Dry")) 

cube_all_data = clean_all_data %>% 
  mutate(across(where(is.numeric), cube_root)) %>% 
  rename_with(where(is.numeric), .fn = ~ paste0("cube_", .x))

## Respiration Histograms ####

color_pallete <- colorRampPalette(colors = c("#D55E00", "#0072B2"))

num_colors <- nlevels(clean_all_data$Treat)

samples_colors <- color_pallete(num_colors)

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Wet_Treatment_Histogram.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(subset(clean_all_data, Treat %in% "Wet"), aes(x = Respiration_Rate_mg_DO_per_kg_per_H)) +
  geom_histogram(fill = "#0072B2")+
  ggtitle("Wet Rates")+
  xlab(expression("Respiration Rate (mg O"[2]*" kg"^- 1*" H"^-1*")"))+
  theme(strip.text = element_text(
    size = 4))+
  ylim(0, 250)+
  theme_bw()

dev.off()

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Cube_Wet_Treatment_Histogram.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(subset(cube_all_data, Treat %in% "Wet"), aes(x = cube_Respiration_Rate_mg_DO_per_kg_per_H)) +
  geom_histogram(fill = "#0072B2")+
  ggtitle("Wet Rates")+
  xlab(expression("Respiration Rate (mg O"[2]*" kg"^- 1*" H"^-1*")"))+
  theme(strip.text = element_text(
    size = 4))+
  ylim(0, 50)+
  theme_bw()

dev.off()


png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Dry_Treatment_Histogram.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(subset(clean_all_data, Treat %in% "Dry"), aes(x = Respiration_Rate_mg_DO_per_kg_per_H)) +
  geom_histogram(fill = "#D55E00")+
  ggtitle("Dry Rates")+
  xlab(expression("Respiration Rate (mg O"[2]*" kg"^- 1*" H"^-1*")"))+
  theme(strip.text = element_text(
    size = 4))+
  ylim(0, 250) + 
  theme_bw()

dev.off()

####

## Effect Size Histogram ####

effect_limits = c(-1500, 1500)

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Median_Effect_Histogram.png"), width = 10, height = 10, units = "in", res = 300)

ggplot(effect_data, aes(x = diff_median_Respiration_Rate_mg_DO_per_kg_per_H))+
  # geom_histogram(binwidth = 0.15, fill = "#009E73")+
  geom_histogram(binwidth = 30, aes(fill = after_stat(x))) +
  scale_fill_gradient2(name = "Effect Size", limits = effect_limits, low = "firebrick2", mid = "goldenrod2",
                       high = "dodgerblue2", midpoint = (max(effect_limits)+min(effect_limits))/2) +
  theme_bw()+
  theme(axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size =18))+
  xlim(c(-1500,1500))+
  ylab("Count\n")+
  xlab(expression("\n Effect Size (Median Wet - Median Dry Rate; mg O"[2]*" kg"^-1*" H"^-1*")"))


dev.off()

cube_effect_data = effect_data %>% 
  mutate(across(where(is.numeric), cube_root)) %>% 
  rename_with(where(is.numeric), .fn = ~ paste0("cube_", .x))

cube_effect_limits <- c(-12, 12)

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Cube_Median_Effect_Histogram.png"), width = 10, height = 10, units = "in", res = 300)

ggplot(cube_effect_data, aes(x = cube_diff_median_Respiration_Rate_mg_DO_per_kg_per_H))+
  # geom_histogram(binwidth = 0.15, fill = "#009E73")+
  geom_histogram(binwidth = 0.5, aes(fill = after_stat(x))) +
  scale_fill_gradient2(name = "Cube Effect Size", limits = cube_effect_limits, low = "firebrick2", mid = "goldenrod2",
                       high = "dodgerblue2", midpoint = (max(cube_effect_limits)+min(cube_effect_limits))/2) +
  theme_bw()+
  theme(axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size =18))+
  xlim(c(-12, 12))+
  ylab("Count\n")+
  xlab(expression("\n Cube Effect Size (Median Wet - Median Dry Rate; mg O"[2]*" kg"^-1*" H"^-1*")"))

dev.off()

####

# James Correlation Matrix ####

cube_effect_data_corr = cube_effect_data %>% 
  column_to_rownames("Sample_ID")%>% 
  select(-c(Rep, cube_diff_median_Respiration_Rate_mg_DO_per_L_per_H, cube_diff_median_Fe_mg_per_L, cube_diff_median_ATP_nanomol_per_L))

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Cube_Median_Effect_Correlation_Matrix.png"), width = 12, height = 12, units = "in", res = 300)

pairs(cube_effect_data_corr, lower.panel = panel.smooth,upper.panel = panel.cor, gap = 0, cex.labels = 1, cex = 1)

dev.off()


## Correlation Ind ####

all_samples_clean_corr <- all_samples_clean %>% 
  dplyr::select(-c(kit, Treat, EC, Log_Mean_Rep_Fe_mg_L, Temp, pH, D50, `Geometric Mean`, `RUSLE Geometric Mean`))

all_samples_corr <- cor(all_samples_clean_corr, method = "spearman")

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/", as.character(Sys.Date()),"_All_Samples_Correlation_Matrix.png"), width = 8, height = 8, units = "in", res = 300)

corrplot(all_samples_corr,type = "upper", tl.col = "black", tl.cex = 1.6, cl.cex = 1.25,  title = "All Samples Correlation")

dev.off()

## Dry Corr Ind ####

all_samples_dry <- all_samples_clean %>%
  na.omit()  %>% 
  filter(!grepl("Wet", Treat)) %>% 
  dplyr::select(-c(kit, Treat, Log_Mean_Rep_Fe_mg_L, Temp, pH, D50, `Geometric Mean`, `RUSLE Geometric Mean`, EC))

all_samples_dry_corr <- cor(all_samples_dry,method = "spearman")

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/", as.character(Sys.Date()),"_All_Dry_Samples_Correlation_Matrix.png"), width = 8, height = 8, units = "in", res = 300)

corrplot(all_samples_dry_corr, type = "upper", tl.col = "black", tl.cex = 1.6, cl.cex = 1.25, title = "All Dry Samples Correlation")

dev.off()

## Wet Correlation Ind ####
all_samples_wet <- all_samples_clean %>%
  na.omit()  %>% 
  filter(!grepl("Dry", Treat)) %>% 
  dplyr::select(-c(kit, Treat, Log_Mean_Rep_Fe_mg_L, Temp, pH, D50, `Geometric Mean`, `RUSLE Geometric Mean`, EC))

all_samples_wet_corr <- cor(all_samples_wet, method = "spearman")

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/", as.character(Sys.Date()),"_All_Wet_Samples_Correlation_Matrix.png"), width = 8, height = 8, units = "in", res = 300)

corrplot(all_samples_wet_corr,type = "upper", tl.col = "black", tl.cex = 1.6, cl.cex = 1.25, title = "All Wet Samples Correlation")

dev.off()

#### Wet/Dry Correlation Matrices ####
wd_list <- list(rem_resp_avg, mean_fe_treat, mean_chem, average_grav)

#merge all data frames in list
wet_dry <- wd_list %>% 
  reduce(merge, by = c("kit", "Treat"), remove = FALSE)

mean_wet_dry <- merge(wet_dry, grn_ssa, by = "kit")

mean_wet_dry_clean <- mean_wet_dry %>% 
  rename(`Average Rate (mg/L)` = Average_Rate) %>% 
  #rename(`Effect Size` = effect) %>% 
    rename(`Mean Initial Gravimetric Water` = average_grav_intial) %>% 
  rename(`Mean Final Gravimetric Water` = average_grav_final) %>% 
  rename(`Mean Lost Gravimetric Water` = average_grav_lost) %>% 
 # rename(`Mean Fe (II) (mg/kg)`  = Mean_Treat_Fe_mg_kg) %>% 
  rename(`Mean Sp. Conducitivity` = Mean_Sp_Conductivity) %>% 
  rename(`Mean pH` = Mean_pH) %>% 
  rename(`Mean Temp.` = Mean_Temperature) %>% 
  rename(`% Fine Sand` = Percent_Fine_Sand) %>% 
  rename(`% Med. Sand` = Percent_Med_Sand) %>% 
  rename(`% Coarse Sand` = Percent_Coarse_Sand) %>% 
  rename(`% Mud` = Percent_Mud) %>% 
  rename(`Geometric Mean` = geom) %>% 
  rename(`RUSLE Geometric Mean` = geom_rusle) %>% 
  rename(`D50` = d50) %>% 
  dplyr::select(-c(#Mean_Treat_Fe_mg_L,
    #Log_Mean_Treat_Fe_mg_kg,
    `RUSLE Geometric Mean`, `Geometric Mean`,  Percent_Clay, `% Fine Sand`,`% Med. Sand`, `% Coarse Sand`,`% Mud`, Percent_Silt, Percent_Tot_Sand)) %>% 
  na.omit()

mean_wet_dry_clean_corr <- mean_wet_dry_clean %>%
  unite(kit_treat, c("kit", "Treat"), sep = "_", remove = TRUE) %>% 
  dplyr::select(-c(EC)) %>% 
  na.omit() %>% 
  remove_rownames %>% 
  column_to_rownames(var = c("kit_treat")) 

mean_samples_corr <- cor(mean_wet_dry_clean_corr, method = "spearman")

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/", as.character(Sys.Date()),"_Mean_Wet_Dry_Samples_Correlation_Matrix.png"), width = 8, height = 8, units = "in", res = 300)

corrplot(mean_samples_corr,type = "upper", tl.col = "black", tl.cex = 1.6, cl.cex = 1.25,  title = "Mean Samples Correlation")

dev.off()

## Mean Wet Corr ####
mean_wet <- mean_wet_dry_clean %>%
  filter(!grepl("Dry", Treat)) %>% 
  remove_rownames %>% 
  column_to_rownames(var = "kit") %>% 
  dplyr::select(-c(Treat, EC))
 
mean_wet_corr <- cor(mean_wet, method = "spearman")

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/", as.character(Sys.Date()),"_Mean_Wet_Treatment_Correlation_Matrix.png"), width = 8, height = 8, units = "in", res = 300)

corrplot(mean_wet_corr, title = "Wet Treatment Correlation", type = "upper", tl.col = "black", tl.cex = 1.6, cl.cex = 1.25)

dev.off()

## Mean Dry Corr ####
mean_dry <-  mean_wet_dry_clean %>%
  filter(!grepl("Wet", Treat)) %>% 
  remove_rownames %>% 
  column_to_rownames(var = "kit") %>% 
  dplyr::select(-c(Treat, EC))

mean_dry_corr <- cor(mean_dry, method = "spearman")

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/", as.character(Sys.Date()),"_Mean_Dry_Treatment_Correlation_Matrix.png"), width = 8, height = 8, units = "in", res = 300)

corrplot(mean_dry_corr, title = "Dry Treatment Correlation", type = "upper", tl.col = "black", tl.cex = 1.6, cl.cex = 1.25)

dev.off()

#### Effect Differences Correlation Matrix ####

effect_list <- list(effect_all, mean_fe_diff, mean_chem_diff, average_grav_lost,grn, ssa_clean)

#merge all data frames in list
effect <- effect_list %>% 
  reduce(merge, by = "kit") %>% 
  mutate(log_ssa = log10(mean_ssa)) %>% 
  mutate(log_mud = log10(Percent_Mud)) %>% 
  rename(`Effect Size` = Effect_Size) %>% 
  rename(`Log Effect Size` = Log_Effect_Size) %>% 
  #rename(`Fe (II) Diff.`  = Fe_Difference_mg_kg) %>% 
  rename(`SpC Diff.` = Sp_Conductivity_Difference) %>% 
  rename(`Temp. Diff.` = Temp_Difference) %>% 
  rename(`pH Diff.` = pH_Difference) %>% 
  rename(`Moisture Diff.` = Final_Gravimetric_Moisture_Difference) %>% 
  rename(`% Fine Sand` = Percent_Fine_Sand) %>% 
  rename(`% Med. Sand` = Percent_Med_Sand) %>% 
  rename(`% Coarse Sand` = Percent_Coarse_Sand) %>% 
  rename(`% Mud` = Percent_Mud) %>% 
  rename(D50 = d50) %>% 
  rename(`Fe Difference (mg/L)` = Fe_Difference_mg_L) %>% 
  na.omit()  %>% 
  remove_rownames %>% 
  column_to_rownames(var = "kit") %>% 
  dplyr::select(-c("Percent_Tot_Sand", "Percent_Silt", "Percent_Clay",D50, EC))%>% 
  dplyr::select(-c(#"% Fine Sand", "% Med. Sand", "% Coarse Sand", "% Mud",
    "geom_rusle", "geom"
    #, "D50"
    ))

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/", as.character(Sys.Date()),"_Log_Effect_Log_SSA.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(effect, aes(y = `Log Effect Size`, x = `log_ssa`)) +
  geom_point() +
  #stat_poly_eq()+
  #stat_poly_line()+
  #geom_smooth(method = "lm")+
  ylab("Log Effect Size") +
  xlab(expression("Log Specific Surface Area (m"^2*" g"^-1*")")) + 
  #xlab("Log % Mud")+
  theme_bw()

dev.off()

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/", as.character(Sys.Date()),"_Effect_SSA.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(effect, aes(y = `Effect Size`, x = `mean_ssa`)) +
  geom_point(size = 4, aes(color = `Effect Size`)) +
  scale_color_gradient2(name = "Effect Size", limits = effect_limits, low = "firebrick2", mid = "goldenrod2",
                       high = "dodgerblue2", midpoint = (max(effect_limits)+min(effect_limits))/2)+
  #stat_poly_eq()+
  #stat_poly_line()+
  #geom_smooth(method = "lm")+
  ylab("Effect Size") +
  xlab(expression("Specific Surface Area (m"^2*" g"^-1*")")) + 
  #xlab("Log % Mud")+
  theme_bw()+
  theme(axis.title = element_text(size=20))

dev.off()

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/", as.character(Sys.Date()),"_Effect_Mud.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(effect, aes(y = `Effect Size`, x = `% Mud`)) +
  geom_point() +
  #stat_poly_eq()+
  #stat_poly_line()+
  #geom_smooth(method = "lm")+
  ylab("Effect Size") +
  xlab(expression("% Mud")) + 
  #xlab("Log % Mud")+
  theme_bw()

dev.off()

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/", as.character(Sys.Date()),"_Log_Effect_Log_Mud.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(effect, aes(y = `Log Effect Size`, x = `log_mud`)) +
  geom_point() +
  #stat_poly_eq()+
  #stat_poly_line()+
  #geom_smooth(method = "lm")+
  ylab("Log Effect Size") +
  xlab(expression("Log % Mud")) + 
  #xlab("Log % Mud")+
  theme_bw()

dev.off()
  
effect_corr <- cor(effect, method = "spearman")

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/", as.character(Sys.Date()),"_Effect_Difference_Correlation_Matrix_Coulson.png"), width = 8, height = 8, units = "in", res = 300)

corrplot(effect_corr, type = 'upper', tl.col = "black", tl.cex = 1.6, cl.cex = 1.25)

dev.off()

rcorr(cbind(effect_corr))




## Effect Size PCA ####

effect_pca <- prcomp(effect, scale = TRUE,
                 center = TRUE, retx = T)

# Summary
summary(effect_pca)

# See the principal components
dim(effect_pca$x)
effect_pca$x

limits = c(-4,
           4)

ind <- get_pca_ind(effect_pca)
ind

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/", as.character(Sys.Date()),"_Effect_PCA.png"), width = 10, height = 10, units = "in", res = 500)

fviz_pca_biplot(effect_pca, col.var = "black",geom = "point"
)+
  geom_point(aes(color = effect$`Effect Size`), size = 3.5)+
  scale_color_gradient2(limits = limits, low = "firebrick2", mid = "goldenrod2",
                        high = "dodgerblue2", midpoint = (max(limits)+min(limits))/2) +
  labs(color = paste0("Wet - Dry Rate"))

dev.off()
