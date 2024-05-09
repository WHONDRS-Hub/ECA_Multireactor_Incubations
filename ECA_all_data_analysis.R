#### Sensitivity Analysis For ECA removals ####
library(tidyverse)
library(corrplot)
library(ggpubr)
library(ggpmisc)
library(factoextra)
library(stringr)

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
                           grepl("D", INC) ~"Dry")) %>% 
  mutate(Fe_mg_per_L = if_else(grepl("SFE_Below", Fe_mg_per_L), "0.001", Fe_mg_per_L)) %>% 
  mutate(Fe_mg_per_kg = if_else(grepl("SFE_Below", Fe_mg_per_kg), "0.003", Fe_mg_per_kg)) %>%
  mutate(Fe_mg_per_L = if_else(grepl("SFE_Above", Fe_mg_per_L), str_extract(Fe_mg_per_L, "(?<=\\|[^|]{1,100}\\|)\\d+\\.\\d+"), Fe_mg_per_L)) %>% 
  mutate(Fe_mg_per_kg = if_else(grepl("SFE_Above", Fe_mg_per_kg), str_extract(Fe_mg_per_kg, "(?<=\\|[^|]{1,100}\\|)\\d+\\.\\d+"), Fe_mg_per_kg)) %>% 
  mutate(Fe_mg_per_L = as.numeric(Fe_mg_per_L)) %>% 
  mutate(Fe_mg_per_kg = as.numeric(Fe_mg_per_kg))

cube_all_data = clean_all_data %>% 
  mutate(across(where(is.numeric), cube_root)) %>% 
  rename_with(where(is.numeric), .fn = ~ paste0("cube_", .x))

## Respiration Histograms ####

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
  xlab(expression("Cube Respiration Rate (mg O"[2]*" kg"^- 1*" H"^-1*")"))+
  theme(strip.text = element_text(
    size = 4))+
  ylim(0, 87.5)+
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

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Cube_Wet_Treatment_Histogram.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(subset(cube_all_data, Treat %in% "Dry"), aes(x = cube_Respiration_Rate_mg_DO_per_kg_per_H)) +
  geom_histogram(fill = "#D55E00")+
  ggtitle("Wet Rates")+
  xlab(expression("Cube Respiration Rate (mg O"[2]*" kg"^- 1*" H"^-1*")"))+
  theme(strip.text = element_text(
    size = 4))+
  ylim(0, 87.5)+
  theme_bw()

dev.off()

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_All_Rates_Histogram.png"), width = 8, height = 8, units = "in", res = 300)

all_cube_hist = ggplot(cube_all_data, aes(x = cube_Respiration_Rate_mg_DO_per_kg_per_H)) +
  geom_histogram(position = "identity", alpha = 0.8, aes(fill = Treat))+
  scale_fill_manual(values = c("#D55E00","#0072B2"))  +
  #ggtitle("Wet Rates")+
  xlab(expression("Cube Respiration Rate (mg O"[2]*" kg"^- 1*" H"^-1*")"))+
  theme(strip.text = element_text(
    size = 4))+
  ylim(0, 87.5)+
  theme_bw()

all_cube_hist

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

cube_effect_hist = ggplot(cube_effect_data, aes(x = cube_diff_median_Respiration_Rate_mg_DO_per_kg_per_H))+
  # geom_histogram(binwidth = 0.15, fill = "#009E73")+
  geom_histogram(binwidth = 0.5, aes(fill = after_stat(x))) +
  scale_fill_gradient2(name = "Cube Effect Size", limits = cube_effect_limits, low = "firebrick2", mid = "goldenrod2",
                       high = "dodgerblue2", midpoint = (max(cube_effect_limits)+min(cube_effect_limits))/2) +
  theme_bw()+
  #theme(axis.title.x = element_text(size = 4),
      #  axis.title.y = element_text(size = 4),
      #  axis.text.x = element_text(size = 4),
      #  axis.text.y = element_text(size =4))+
  xlim(c(-12, 12))+
  ylab("count\n")+
  xlab(expression("\n Cube Effect Size (Median Wet - Median Dry Rate; mg O"[2]*" kg"^-1*" H"^-1*")"))
  

cube_effect_hist

dev.off()

####


## Effect Size + Rate Combined Figure

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Combined_Histogram.png"), width = 8, height = 4, units = "in", res = 300)

cube_effect_hist_new = cube_effect_hist +
  theme(legend.position = c(0.85, 0.8),
        legend.key.size = unit(0.15, "in"), 
        legend.title = element_text(size = 8),
        axis.title.x = element_text(size = 10)) + 
  xlab(expression(atop("\n Cube Effect Size", "(Median Wet - Median Dry Rate; mg O"[2]*" kg"^-1*" H"^-1*")"))) 
  

all_cube_hist_new = all_cube_hist + 
  theme(legend.position = c(0.85, 0.8), 
        legend.key.size = unit(0.15, "in"), 
        legend.title = element_text(size = 8),
        axis.title.x = element_text(size = 10)) +
  guides(fill = guide_legend(title="Treatment")) + 
  xlab(expression(atop("\n Cube Respiration Rate", "(mg O"[2]*" kg"^-1*" H"^-1*")")))

combine_hist = ggarrange(all_cube_hist_new, cube_effect_hist_new, labels = c("A", "B"), label.x = c(0.15, 0.175), label.y = c(0.97, 0.97))

combine_hist

dev.off()

# James Correlation Matrix ####
all_samples_cube_corr <- cube_all_data %>% 
  dplyr::select(-c(Sample_ID, INC, Treat, Methods_Deviation, cube_Incubation_Water_Mass_g, cube_Dry_Sediment_Mass_g, cube_Initial_Water_mass_g, cube_Final_Water_mass_g, cube_Respiration_Rate_mg_DO_per_L_per_H, cube_Fe_mg_per_L, cube_ATP_nanomol_per_L)) %>% 
  relocate(cube_Respiration_Rate_mg_DO_per_kg_per_H, .before = cube_SpC) %>% 
  column_to_rownames("Sample_Name") %>% 
  filter(!is.na(cube_mean_ssa))


cube_all_samples_corr <- cor(all_samples_cube_corr, method = "spearman")

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Cube_All_Samples_Correlation_Matrix.png"), width = 12, height = 12, units = "in", res = 300)

pairs(all_samples_cube_corr, lower.panel = panel.smooth,upper.panel = panel.cor, gap = 0, cex.labels = 0.5, cex = 1)

#corrplot(cube_all_samples_corr,type = "upper", tl.col = "black", tl.cex = 1.6, cl.cex = 1.25,  title = "All Samples Correlation")

dev.off()

## NEED TO BE UPDATED - WET/DRY Individual Matrices #### ####

## Dry Corr Ind ###

all_samples_dry <- all_samples_clean %>%
  na.omit()  %>% 
  filter(!grepl("Wet", Treat)) %>% 
  dplyr::select(-c(kit, Treat, Log_Mean_Rep_Fe_mg_L, Temp, pH, D50, `Geometric Mean`, `RUSLE Geometric Mean`, EC))

all_samples_dry_corr <- cor(all_samples_dry,method = "spearman")

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/", as.character(Sys.Date()),"_All_Dry_Samples_Correlation_Matrix.png"), width = 8, height = 8, units = "in", res = 300)

#corrplot(all_samples_dry_corr, type = "upper", tl.col = "black", tl.cex = 1.6, cl.cex = 1.25, title = "All Dry Samples Correlation")

dev.off()

## Wet Correlation Ind ###
all_samples_wet <- all_samples_clean %>%
  na.omit()  %>% 
  filter(!grepl("Dry", Treat)) %>% 
  dplyr::select(-c(kit, Treat, Log_Mean_Rep_Fe_mg_L, Temp, pH, D50, `Geometric Mean`, `RUSLE Geometric Mean`, EC))

all_samples_wet_corr <- cor(all_samples_wet, method = "spearman")

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/", as.character(Sys.Date()),"_All_Wet_Samples_Correlation_Matrix.png"), width = 8, height = 8, units = "in", res = 300)

corrplot(all_samples_wet_corr,type = "upper", tl.col = "black", tl.cex = 1.6, cl.cex = 1.25, title = "All Wet Samples Correlation")

dev.off()



## NEED TO BE UPDATED - WET/DRY Median Matrices ####
## Effect Size ####
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
  relocate(Cube_Effect_Size, .before = Cube_SpC_Diff) %>% 
  filter(!is.na(Cube_SSA))
  
png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Cube_Median_Effect_Correlation_Matrix.png"), width = 12, height = 12, units = "in", res = 300)

pairs(cube_effect_data_corr, lower.panel = panel.smooth,upper.panel = panel.cor, gap = 0, cex.labels = 0.5, cex = 1)

dev.off()

## Scatter Plots ####

##Fe vs. Effect ##

fe_cube_out = cube_effect_data_corr %>% 
  filter(Cube_Fe_mg_kg_Diff > -1)

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Cube_Median_Effect_vs_Fe_Scatter.png"), width = 6, height = 6, units = "in", res = 300)

fe_cube = ggplot(cube_effect_data_corr, aes(x = Cube_Fe_mg_kg_Diff, y = Cube_Effect_Size)) +
  geom_point() +
  theme_bw() +
  stat_cor(data = fe_cube_out, size = 5, digits = 2, aes(label = paste(..rr.label.., ..p.label.., sep = "~`;`~")))+
  stat_poly_line(data = fe_out, se = FALSE)+
  xlab("Cube Root Fe (II) (mg/kg) Difference (Wet - Dry)") +
  ylab("Cube Root Effect Size (mg/kg) (Wet - Dry)")+ 
  theme(text = element_text(size = 13)) 

dev.off()

fe_out = effect_data %>% 
  filter(diff_median_Fe_mg_per_kg > -1)

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Median_Effect_vs_Fe_Scatter.png"), width = 6, height = 6, units = "in", res = 300)

ggplot(effect_data, aes(x = diff_median_Fe_mg_per_kg, y = diff_median_Respiration_Rate_mg_DO_per_kg_per_H)) +
  geom_point() +
  theme_bw() +
  stat_cor(data = fe_out, size = 5, digits = 2, aes(label = paste(..rr.label.., ..p.label.., sep = "~`;`~")))+
  stat_poly_line(data = fe_out, se = FALSE)+
  xlab("Fe (II) (mg/kg) Difference (Wet - Dry)") +
  ylab("Effect Size (mg/kg) (Wet - Dry)")+ 
  theme(text = element_text(size = 13)) 

dev.off()

##ATP vs. Effect ##

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Cube_Median_Effect_vs_ATP_Scatter.png"), width = 6, height = 6, units = "in", res = 300)

atp_cube = ggplot(cube_effect_data_corr, aes(x = Cube_ATP_pmol_g_Diff, y = Cube_Effect_Size)) +
  geom_point() +
  theme_bw() +
  stat_cor(data = cube_effect_data_corr, size = 5, digits = 2, aes(label = paste(..rr.label.., ..p.label.., sep = "~`;`~")))+
  stat_poly_line(data = cube_effect_data_corr, se = FALSE)+
  xlab("Cube Root ATP (picomol/g) Difference (Wet - Dry)") +
  ylab("Cube Root Effect Size (mg/kg) (Wet - Dry)")+ 
  theme(text = element_text(size = 13)) 

dev.off()

##Fine Sand vs. Effect ##

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Cube_Median_Effect_vs_ATP_Scatter.png"), width = 6, height = 6, units = "in", res = 300)

fs_cube = ggplot(cube_effect_data_corr, aes(x = Cube_Fine_Sand, y = Cube_Effect_Size)) +
  geom_point() +
  theme_bw() +
  stat_cor(data = cube_effect_data_corr, size = 5, digits = 2, aes(label = paste(..rr.label.., ..p.label.., sep = "~`;`~")))+
  stat_poly_line(data = cube_effect_data_corr, se = FALSE)+
  xlab("Cube Root Fine Sand (%)") +
  ylab("Cube Root Effect Size (mg/kg) (Wet - Dry)")+ 
  theme(text = element_text(size = 13)) 

dev.off()


## Effect Size PCA ####

effect_data_clean = effect_data %>% 
  filter(!is.na(median_mean_ssa)) %>% 
  select(-c(Rep, diff_median_ATP_nanomol_per_L, diff_median_Respiration_Rate_mg_DO_per_L_per_H, diff_median_Fe_mg_per_L)) %>% 
  column_to_rownames("Sample_ID") %>% 
  rename(SpC_Diff = diff_median_SpC) %>% 
  rename(pH_Diff = diff_median_pH) %>%
  rename(Temp_Diff = diff_median_Temp) %>%
  rename(Effect_Size = diff_median_Respiration_Rate_mg_DO_per_kg_per_H) %>%
  rename(Fe_mg_kg_Diff = diff_median_Fe_mg_per_kg) %>%
  rename(InGravMoi_Diff = diff_median_Initial_Gravimetric_Moisture) %>%
  rename(FinGravMoi_Diff = diff_median_Final_Gravimetric_Moisture) %>%
  rename(LostGravMoi_Diff = diff_median_Lost_Gravimetric_Water) %>%
  rename(ATP_pmol_g_Diff = diff_median_ATP_picomol_per_g) %>%
  rename(Fine_Sand = median_Percent_Fine_Sand) %>%
  rename(Med_Sand = median_Percent_Med_Sand) %>%
  rename(Coarse_Sand = median_Percent_Coarse_Sand) %>%
  rename(Tot_Sand = median_Percent_Tot_Sand) %>%
  rename(Silt = median_Percent_Silt) %>%
  rename(Clay = median_Percent_Clay) %>%
  rename(SSA = median_mean_ssa) %>% 
  relocate(Effect_Size, .before = SpC_Diff)

effect_pca <- prcomp(effect_data_clean, scale = TRUE,
                 center = TRUE, retx = T)

# Summary
summary(effect_pca)

# See the principal components
dim(effect_pca$x)
effect_pca$x

limits = c(-1450,
           1450)

ind <- get_pca_ind(effect_pca)
ind

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Median_Effect_PCA.png"), width = 8, height = 8, units = "in", res = 300)

pca = fviz_pca_biplot(effect_pca, col.var = "black",geom = "point"
)+
  geom_point(aes(color = effect_data_clean$Effect_Size), size = 3.5)+
  scale_color_gradient2(limits = limits, low = "firebrick2", mid = "goldenrod2",
                        high = "dodgerblue2", midpoint = (max(limits)+min(limits))/2) +
  labs(color = paste0("Wet - Dry Rate"))

pca

dev.off()

## PCA with Scatter Plots ####

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Combined_PCA.png"), width = 12, height = 8, units = "in", res = 300)

pca_new = pca + theme(legend.position = c(0.85, 0.2), 
      legend.key.size = unit(0.25, "in"), 
      legend.title = element_text(size = 8),
      axis.title.x = element_text(size = 10))

fe_cube_new = fe_cube + 
  ylab("") + 
  theme(axis.title.x = element_text(size = 10))


atp_cube_new = atp_cube + 
  ylab("") + 
  theme(axis.title.x = element_text(size = 10))

fs_cube_new = fs_cube + 
  ylab("") + 
  theme(axis.title.x = element_text(size = 10))

combine_scatter = ggarrange(fe_cube_new, atp_cube_new, fs_cube_new, common.legend = TRUE, nrow = 3, labels = c("B", "C", "D"), label.x = c(0.9, 0.9, 0.9), label.y = c(0.3, 0.3, 0.3), heights = c(1,1,1)) 

  annotate_scatter = annotate_figure(combine_scatter, left = text_grob("Cube Root Effect Size (mg/kg) (Wet - Dry)", x = 0.75, y = 0.5, size = 15, vjust = 1, rot = 90))


combine_pca = ggarrange(pca_new, annotate_scatter, labels = c("A"), ncol = 2, label.x = c(0.08), label.y = c(0.95), widths = c(2, 1), heights = c(2,2)) 

combine_pca

dev.off()
