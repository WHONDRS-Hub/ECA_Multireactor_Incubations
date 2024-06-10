#### Sensitivity Analysis For ECA removals ####
library(tidyverse)
library(corrplot)
library(ggpubr)
library(ggpmisc)
library(factoextra)
library(stringr)
library(glmnet)
library(pdftools)
library(magick)


rm(list=ls());graphics.off()
#set.seed(26)

#### Read in Data
##Switch to published DP when ready

#Individual samples 
all_data <- read.csv("C:/Github/ECA_Multireactor_Incubations/Data/Cleaned Data/All_ECA_Data_05-29-2024.csv",header = TRUE) 

#Summary Data 

sum_data <- read.csv("C:/Github/ECA_Multireactor_Incubations/Data/Cleaned Data/2024-05-29_Medians_ECA_Data.csv",header = TRUE) %>% 
  select(-c(X))

#Effect Size Data
effect_data <- read.csv("C:/Github/ECA_Multireactor_Incubations/Data/Cleaned Data/2024-05-29_Effect_Median_ECA_Data.csv",header = TRUE) %>% 
 dplyr::select(-c(X)) 

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

# png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Wet_Treatment_Histogram.png"), width = 8, height = 8, units = "in", res = 300)
# 
# ggplot(subset(clean_all_data, Treat %in% "Wet"), aes(x = Respiration_Rate_mg_DO_per_kg_per_H)) +
#   geom_histogram(fill = "#0072B2")+
#   ggtitle("Wet Rates")+
#   xlab(expression("Respiration Rate (mg O"[2]*" kg"^- 1*" H"^-1*")"))+
#   theme(strip.text = element_text(
#     size = 4))+
#   ylim(0, 250)+
#   theme_bw()
# 
# dev.off()

# png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Cube_Wet_Treatment_Histogram.png"), width = 8, height = 8, units = "in", res = 300)
# 
# ggplot(subset(cube_all_data, Treat %in% "Wet"), aes(x = cube_Respiration_Rate_mg_DO_per_kg_per_H)) +
#   geom_histogram(fill = "#0072B2")+
#   ggtitle("Wet Rates")+
#   xlab(expression("Cube Respiration Rate (mg O"[2]*" kg"^- 1*" H"^-1*")"))+
#   theme(strip.text = element_text(
#     size = 4))+
#   ylim(0, 87.5)+
#   theme_bw()
# 
# dev.off()


# png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Dry_Treatment_Histogram.png"), width = 8, height = 8, units = "in", res = 300)
# 
# ggplot(subset(clean_all_data, Treat %in% "Dry"), aes(x = Respiration_Rate_mg_DO_per_kg_per_H)) +
#   geom_histogram(fill = "#D55E00")+
#   ggtitle("Dry Rates")+
#   xlab(expression("Respiration Rate (mg O"[2]*" kg"^- 1*" H"^-1*")"))+
#   theme(strip.text = element_text(
#     size = 4))+
#   ylim(0, 250) + 
#   theme_bw()
# 
# dev.off()

# png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Cube_Dry_Treatment_Histogram.png"), width = 8, height = 8, units = "in", res = 300)
# 
# ggplot(subset(cube_all_data, Treat %in% "Dry"), aes(x = cube_Respiration_Rate_mg_DO_per_kg_per_H)) +
#   geom_histogram(fill = "#D55E00")+
#   ggtitle("Wet Rates")+
#   xlab(expression("Cube Respiration Rate (mg O"[2]*" kg"^- 1*" H"^-1*")"))+
#   theme(strip.text = element_text(
#     size = 4))+
#   ylim(0, 87.5)+
#   theme_bw()
# 
# dev.off()

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_All_Rates_Histogram.png"), width = 4, height = 4, units = "in", res = 300)

all_cube_hist = ggplot(cube_all_data, aes(x = cube_Respiration_Rate_mg_DO_per_kg_per_H)) +
  geom_histogram(position = "identity", alpha = 0.8, aes(fill = Treat))+
  scale_fill_manual(values = c("#D55E00","#0072B2"))  +
  #ggtitle("Wet Rates")+
  theme(strip.text = element_text(
    size = 4))+
  ylim(0, 87.5)+
  theme_bw() + 
  theme(legend.position = c(0.85, 0.8), 
        legend.key.size = unit(0.15, "in"), 
        legend.title = element_text(size = 8),
        axis.title.x = element_text(size = 10)) +
  guides(fill = guide_legend(title="Treatment")) + 
  xlab(expression(atop("\n Cubed Root Respiration Rate", "(mg O"[2]*" kg"^-1*" H"^-1*")"))) +
  ylab("Count")

all_cube_hist

dev.off()


####

## Effect Size Histogram ####

# effect_limits = c(-1500, 1500)
# 
# png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Median_Effect_Histogram.png"), width = 10, height = 10, units = "in", res = 300)
# 
# ggplot(effect_data, aes(x = diff_median_Respiration_Rate_mg_DO_per_kg_per_H))+
#   # geom_histogram(binwidth = 0.15, fill = "#009E73")+
#   geom_histogram(binwidth = 30, aes(fill = after_stat(x))) +
#   scale_fill_gradient2(name = "Effect Size", limits = effect_limits, low = "firebrick2", mid = "goldenrod2",
#                        high = "dodgerblue2", midpoint = (max(effect_limits)+min(effect_limits))/2) +
#   theme_bw()+
#   theme(axis.title.x = element_text(size = 24),
#         axis.title.y = element_text(size = 24),
#         axis.text.x = element_text(size = 18),
#         axis.text.y = element_text(size =18))+
#   xlim(c(-1500,1500))+
#   ylab("Count\n")+
#   xlab(expression("\n Effect Size (Median Wet - Median Dry Rate; mg O"[2]*" kg"^-1*" H"^-1*")"))
# 
# 
# dev.off()

## EFFECT CUBE ROOT
cube_effect_data = effect_data %>% 
  mutate(across(where(is.numeric), cube_root)) %>% 
  rename_with(where(is.numeric), .fn = ~ paste0("cube_", .x))

cube_effect_limits <- c(-12, 12)

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Cube_Median_Effect_Histogram.png"), width = 4, height = 4, units = "in", res = 300)

cube_effect_hist = ggplot(cube_effect_data, aes(x = cube_diff_median_Respiration_Rate_mg_DO_per_kg_per_H))+
  # geom_histogram(binwidth = 0.15, fill = "#009E73")+
  geom_histogram(binwidth = 0.5, aes(fill = after_stat(x))) +
  scale_fill_gradient2(name = "Cubed Root Effect Size", limits = cube_effect_limits, low = "firebrick2", mid = "goldenrod2",
                       high = "dodgerblue2", midpoint = (max(cube_effect_limits)+min(cube_effect_limits))/2) +
  theme_bw()+
  #theme(axis.title.x = element_text(size = 4),
      #  axis.title.y = element_text(size = 4),
      #  axis.text.x = element_text(size = 4),
      #  axis.text.y = element_text(size =4))+
  xlim(c(-12, 12))+
  ylab("Count\n")+
    theme(legend.position = c(0.8, 0.8),
        legend.key.size = unit(0.15, "in"), 
        legend.title = element_text(size = 8),
        axis.title.x = element_text(size = 10)) + 
  xlab(expression(atop("\n Cubed Root Effect Size", "(Median Wet - Median Dry Rate; mg O"[2]*" kg"^-1*" H"^-1*")"))) 
  

cube_effect_hist

dev.off()

####

## Effect Size 

# png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Combined_Histogram.png"), width = 8, height = 4, units = "in", res = 300)
# 
# cube_effect_hist_new = cube_effect_hist +
#   theme(legend.position = c(0.8, 0.8),
#         legend.key.size = unit(0.15, "in"), 
#         legend.title = element_text(size = 8),
#         axis.title.x = element_text(size = 10)) + 
#   xlab(expression(atop("\n Cubed Root Effect Size", "(Median Wet - Median Dry Rate; mg O"[2]*" kg"^-1*" H"^-1*")"))) 
#   
# 
# all_cube_hist_new = all_cube_hist + 
#   theme(legend.position = c(0.85, 0.8), 
#         legend.key.size = unit(0.15, "in"), 
#         legend.title = element_text(size = 8),
#         axis.title.x = element_text(size = 10)) +
#   guides(fill = guide_legend(title="Treatment")) + 
#   xlab(expression(atop("\n Cubed Root Respiration Rate", "(mg O"[2]*" kg"^-1*" H"^-1*")")))

#combine_hist = ggarrange(all_cube_hist_new, cube_effect_hist_new, labels = c("A", "B", "C"), label.x = c(0.15, 0.15, 0.15), label.y = c(0.97, 0.97, 0.97))

# combine_hist
# 
# dev.off()


## Effect Size + Rate + Map Combined Figure

# Read in Map .pdf

map_path = "C:/GitHub/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/Map/Ecoregion_EffectSize_Map.pdf"

map_image = image_read_pdf(map_path, density = 300)

scale_map_image = image_scale(map_image, "80%")

scale_map_label_image = image_annotate(scale_map_image, "C", size = 15, location = "+10+80", color = "black")

# Read in Rate histograms figure .png

rate_image = image_read("C:/GitHub/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/2024-06-10_All_Rates_Histogram.png")

rate_label_image = image_annotate(rate_image, "A", size = 65, location = "+30+20", color = "black")



# Read in Effect Size Figure .pdf

effect_image = image_read("C:/GitHub/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/2024-06-10_Cube_Median_Effect_Histogram.png")

effect_label_image = image_annotate(effect_image, "B", size = 65, location = "+30+20", color = "black")


com_image = image_append(c(rate_label_image, effect_label_image))

com_map_image = image_append(c(com_image, scale_map_label_image), stack = TRUE)

image_write(com_map_image, path = "C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/2024-06-10_Combined_map.png")


# James Correlation Matrix ####
all_samples_cube_corr <- cube_all_data %>% 
  dplyr::select(-c(Sample_ID, INC, Treat, Methods_Deviation, cube_Incubation_Water_Mass_g, cube_Dry_Sediment_Mass_g, cube_Initial_Water_mass_g, cube_Final_Water_mass_g, cube_Respiration_Rate_mg_DO_per_L_per_H, cube_Fe_mg_per_L, cube_ATP_nanomol_per_L)) %>% 
  relocate(cube_Respiration_Rate_mg_DO_per_kg_per_H, .before = cube_SpC) %>% 
  column_to_rownames("Sample_Name") 


cube_all_samples_corr <- cor(all_samples_cube_corr, method = "spearman")

# png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Cube_All_Samples_Correlation_Matrix.png"), width = 12, height = 12, units = "in", res = 300)
# 
#pairs(all_samples_cube_corr, lower.panel = panel.smooth,upper.panel = panel.cor, gap = 0, cex.labels = 0.5, cex = 1)
# 
#corrplot(cube_all_samples_corr,type = "upper", tl.col = "black", tl.cex = 1.6, cl.cex = 1.25,  title = "All Samples Correlation")
# 
# dev.off()


## Effect Size ####
cube_effect_data_corr = cube_effect_data %>% 
  column_to_rownames("Sample_ID") %>% 
  dplyr::select(-c(Rep, cube_diff_median_Respiration_Rate_mg_DO_per_L_per_H, cube_diff_median_Fe_mg_per_L, cube_diff_median_ATP_nanomol_per_L, cube_diff_median_NPOC_mg_C_per_L, cube_diff_median_TN_mg_N_per_L)) %>% 
  rename(Cube_SpC_Diff = cube_diff_median_SpC) %>% 
  rename(Cube_pH_Diff = cube_diff_median_pH) %>%
  rename(Cube_Temp_Diff = cube_diff_median_Temp) %>%
  rename(Cube_Effect_Size = cube_diff_median_Respiration_Rate_mg_DO_per_kg_per_H) %>%
  rename(Cube_Fe_mg_kg_Diff = cube_diff_median_Fe_mg_per_kg) %>%
  rename(Cube_Dry_InGravMoi= cube_median_Dry_Initial_Gravimetric_Moisture) %>%
  rename(Cube_Dry_FinGravMoi = cube_median_Dry_Final_Gravimetric_Moisture) %>%
  rename(Cube_Dry_LostGravMoi = cube_median_Dry_Lost_Gravimetric_Moisture) %>%
  rename(Cube_ATP_pmol_g_Diff = cube_diff_median_ATP_picomol_per_g) %>%
  rename(Cube_NPOC_mg_kg_Diff = cube_diff_median_NPOC_mg_C_per_kg) %>% 
  rename(Cube_TN_mg_kg_Diff = cube_diff_median_TN_mg_N_per_kg) %>% 
  rename(Cube_TN_solid_perc_Diff = cube_diff_median_tn_percent) %>% 
  rename(Cube_TOC_solid_perc_Diff = cube_diff_median_toc_percent) %>% 
  rename(Cube_Fine_Sand = cube_median_Percent_Fine_Sand) %>%
  rename(Cube_Med_Sand = cube_median_Percent_Med_Sand) %>%
  rename(Cube_Coarse_Sand = cube_median_Percent_Coarse_Sand) %>%
  rename(Cube_Tot_Sand = cube_median_Percent_Tot_Sand) %>%
  rename(Cube_Silt = cube_median_Percent_Silt) %>%
  rename(Cube_Clay = cube_median_Percent_Clay) %>%
  rename(Cube_SSA = cube_median_mean_ssa) %>% 
  relocate(Cube_Effect_Size, .before = Cube_SpC_Diff) 

cube_effect_samples_corr <- cor(cube_effect_data_corr, method = "spearman")
  
 png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Cube_Median_Effect_Correlation_Matrix.png"), width = 12, height = 12, units = "in", res = 300)
# 
# pairs(cube_effect_data_corr, lower.panel = panel.smooth,upper.panel = panel.cor, gap = 0, cex.labels = 0.5, cex = 1)

corrplot(cube_effect_samples_corr,type = "upper", method = "number", tl.col = "black", tl.cex = 1.6, cl.cex = 1.25,  title = "Effect Samples Correlation")
# 
 dev.off()

## Scatter Plots ####

##Fe vs. Effect ##

fe_cube_out = cube_effect_data_corr %>% 
  filter(Cube_Fe_mg_kg_Diff > -1)

#png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Cube_Median_Effect_vs_Fe_Scatter.png"), width = 6, height = 6, units = "in", res = 300)

fe_cube = ggplot(cube_effect_data_corr, aes(x = Cube_Fe_mg_kg_Diff, y = Cube_Effect_Size)) +
  geom_point() +
  theme_bw() +
  #stat_cor(data = fe_cube_out, size = 5, digits = 2, aes(label = paste(..rr.label.., ..p.label.., sep = "~`;`~")))+
  stat_cor(data = fe_cube_out, label.x = -2.5, label.y = 11, size = 4, digits = 2, aes(label = paste(..rr.label..)))+
  stat_cor(data = fe_cube_out, label.x = -2.5, label.y = 9.5, size = 4, digits = 2, aes(label = paste(..p.label..)))+
  stat_poly_line(data = fe_cube_out, se = FALSE)+
  xlab("Cubed Root Fe (II) (mg/kg) Difference (Wet - Dry)") +
  ylab("Cubed Root Effect Size (mg/kg) (Wet - Dry)")+ 
  theme(text = element_text(size = 12)) 

#dev.off()

# fe_out = effect_data %>% 
#   filter(diff_median_Fe_mg_per_kg > -1)
# 
# png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Median_Effect_vs_Fe_Scatter.png"), width = 6, height = 6, units = "in", res = 300)
# 
# ggplot(effect_data, aes(x = diff_median_Fe_mg_per_kg, y = diff_median_Respiration_Rate_mg_DO_per_kg_per_H)) +
#   geom_point() +
#   theme_bw() +
#   stat_cor(data = fe_out, size = 5, digits = 2, aes(label = paste(..rr.label.., ..p.label.., sep = "~`;`~")))+
#   stat_poly_line(data = fe_out, se = FALSE)+
#   xlab("Fe (II) (mg/kg) Difference (Wet - Dry)") +
#   ylab("Effect Size (mg/kg) (Wet - Dry)")+ 
#   theme(text = element_text(size = 13)) 
# 
# dev.off()

##ATP vs. Effect ##

# png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Cube_Median_Effect_vs_ATP_Scatter.png"), width = 6, height = 6, units = "in", res = 300)
# 
atp_cube = ggplot(cube_effect_data_corr, aes(x = Cube_ATP_pmol_g_Diff, y = Cube_Effect_Size)) +
  geom_point() +
  theme_bw() +
  stat_cor(data = cube_effect_data_corr, size = 5, digits = 2, aes(label = paste(..rr.label.., ..p.label.., sep = "~`;`~")))+
  stat_poly_line(data = cube_effect_data_corr, se = FALSE)+
  xlab("Cube Root ATP (picomol/g) Difference (Wet - Dry)") +
  ylab("Cube Root Effect Size (mg/kg) (Wet - Dry)")+
  theme(text = element_text(size = 13))

# dev.off()

##Fine Sand vs. Effect ##

#png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Cube_Median_Effect_vs_ATP_Scatter.png"), width = 6, height = 6, units = "in", res = 300)

fs_cube = ggplot(cube_effect_data_corr, aes(x = Cube_Fine_Sand, y = Cube_Effect_Size)) +
  geom_point() +
  theme_bw() +
  #stat_cor(data = cube_effect_data_corr, size = 5, digits = 2, aes(label = paste(..rr.label.., ..p.label.., sep = "~`\n`~")))+ #sep = "~`;`~"
  stat_cor(data = cube_effect_data_corr, label.x = 0.9, label.y = 11.5, size = 4, digits = 2, aes(label = paste(..rr.label..)))+
  stat_cor(data = cube_effect_data_corr, label.x = 0.9, label.y = 11, size = 4, digits = 2, aes(label = paste(..p.label..)))+
  stat_poly_line(data = cube_effect_data_corr, se = FALSE)+
  xlab("Cubed Root Fine Sand (%)") +
  ylab("Cubed Root Effect Size (mg/kg) (Wet - Dry)")+ 
  theme(text = element_text(size = 12)) 

# silt_cube = ggplot(cube_effect_data_corr, aes(x = Cube_Silt, y = Cube_Effect_Size)) +
#   geom_point() +
#   theme_bw() +
#   #stat_cor(data = cube_effect_data_corr, size = 5, digits = 2, aes(label = paste(..rr.label.., ..p.label.., sep = "~`\n`~")))+ #sep = "~`;`~"
#   stat_cor(data = cube_effect_data_corr, label.x = 0.9, label.y = 11.5, size = 4, digits = 2, aes(label = paste(..rr.label..)))+
#   stat_cor(data = cube_effect_data_corr, label.x = 0.9, label.y = 11, size = 4, digits = 2, aes(label = paste(..p.label..)))+
#   stat_poly_line(data = cube_effect_data_corr, se = FALSE)+
#   xlab("Cubed Root Silt (%)") +
#   ylab("Cubed Root Effect Size (mg/kg) (Wet - Dry)")+ 
#   theme(text = element_text(size = 12))

# tot_sand_cube = ggplot(cube_effect_data_corr, aes(x = Cube_Tot_Sand, y = Cube_Effect_Size)) +
#   geom_point() +
#   theme_bw() +
#   #stat_cor(data = cube_effect_data_corr, size = 5, digits = 2, aes(label = paste(..rr.label.., ..p.label.., sep = "~`\n`~")))+ #sep = "~`;`~"
#   stat_cor(data = cube_effect_data_corr, label.x = 3.2, label.y = 11.5, size = 4, digits = 2, aes(label = paste(..rr.label..)))+
#   stat_cor(data = cube_effect_data_corr, label.x = 3.2, label.y = 11, size = 4, digits = 2, aes(label = paste(..p.label..)))+
#   stat_poly_line(data = cube_effect_data_corr, se = FALSE)+
#   xlab("Cubed Root Total Sand (%)") +
#   ylab("Cubed Root Effect Size (mg/kg) (Wet - Dry)")+ 
#   theme(text = element_text(size = 12))

# clay_cube = ggplot(cube_effect_data_corr, aes(x = Cube_Clay, y = Cube_Effect_Size)) +
#   geom_point() +
#   theme_bw() +
#   #stat_cor(data = cube_effect_data_corr, size = 5, digits = 2, aes(label = paste(..rr.label.., ..p.label.., sep = "~`\n`~")))+ #sep = "~`;`~"
#   stat_cor(data = cube_effect_data_corr, label.x = 3.2, label.y = 11.5, size = 4, digits = 2, aes(label = paste(..rr.label..)))+
#   stat_cor(data = cube_effect_data_corr, label.x = 3.2, label.y = 11, size = 4, digits = 2, aes(label = paste(..p.label..)))+
#   stat_poly_line(data = cube_effect_data_corr, se = FALSE)+
#   xlab("Cubed Root Clay (%)") +
#   ylab("Cubed Root Effect Size (mg/kg) (Wet - Dry)")+ 
#   theme(text = element_text(size = 12))

# coarse_sand_cube = ggplot(cube_effect_data_corr, aes(x = Cube_Coarse_Sand, y = Cube_Effect_Size)) +
#   geom_point() +
#   theme_bw() +
#   #stat_cor(data = cube_effect_data_corr, size = 5, digits = 2, aes(label = paste(..rr.label.., ..p.label.., sep = "~`\n`~")))+ #sep = "~`;`~"
#   stat_cor(data = cube_effect_data_corr, label.x = 3.2, label.y = 11.5, size = 4, digits = 2, aes(label = paste(..rr.label..)))+
#   stat_cor(data = cube_effect_data_corr, label.x = 3.2, label.y = 11, size = 4, digits = 2, aes(label = paste(..p.label..)))+
#   stat_poly_line(data = cube_effect_data_corr, se = FALSE)+
#   xlab("Cubed Root Coarse Sand (%)") +
#   ylab("Cubed Root Effect Size (mg/kg) (Wet - Dry)")+ 
#   theme(text = element_text(size = 12))

#dev.off()

# cn_cube = ggplot(cube_effect_data_corr, aes(x = Cube_TOC_solid_perc_Diff, y = Cube_Effect_Size)) +
#   geom_point() +
#   theme_bw() +
#   #stat_cor(data = cube_effect_data_corr, size = 5, digits = 2, aes(label = paste(..rr.label.., ..p.label.., sep = "~`\n`~")))+ #sep = "~`;`~"
#   stat_cor(data = cube_effect_data_corr, label.x = -1.5, label.y = 11.5, size = 4, digits = 2, aes(label = paste(..rr.label..)))+
#   stat_cor(data = cube_effect_data_corr, label.x = -1.5, label.y = 11, size = 4, digits = 2, aes(label = paste(..p.label..)))+
#   stat_poly_line(data = cube_effect_data_corr, se = FALSE)+
#   xlab("Cubed Root TOC Difference Solid (%)") +
#   ylab("Cubed Root Effect Size (mg/kg) (Wet - Dry)")+ 
#   theme(text = element_text(size = 12))


png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Cube_Median_Effect_vs_Fin_Grav_Moi_Scatter.png"), width = 6, height = 6, units = "in", res = 300)

grav_cube = ggplot(cube_effect_data_corr, aes(x = Cube_Dry_LostGravMoi, y = Cube_Effect_Size)) +
  geom_point() +
  theme_bw() +
  #stat_cor(data = cube_effect_data_corr, size = 5, digits = 2, aes(label = paste(..rr.label.., ..p.label.., sep = "~`;`~")))+
  stat_cor(data = cube_effect_data_corr, label.x = 0.425, label.y = 11, size = 4, digits = 2, aes(label = paste(..rr.label..)))+
  stat_cor(data = cube_effect_data_corr, label.x = 0.425, label.y = 10.5, size = 4, digits = 2, aes(label = paste(..p.label..)))+
  stat_poly_line(data = cube_effect_data_corr, se = FALSE)+
  xlab("Cubed Root Lost Gravimetric Moisture Difference (Dry)") +
  ylab("Cubed Root Effect Size (mg/kg) (Wet - Dry)")+ 
  theme(text = element_text(size = 12)) 

dev.off()



## Effect Size PCA ####

effect_data_clean = effect_data %>% 
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
  relocate(Effect_Size, .before = SpC_Diff) %>% 
  filter(Fe_mg_kg_Diff > -1)

effect_pca <- prcomp(effect_data_clean, scale = TRUE,
                 center = TRUE, retx = T)

# Summary
summary(effect_pca)

# See the principal components
dim(effect_pca$x)
effect_pca$x

## PCA into LASSO ####
#Put prop. of variation in DF for LASSO

pca_lasso_data = effect_data_clean %>% 
  select(-c(Effect_Size))

effect_pca_lasso <- prcomp(pca_lasso_data, scale = TRUE,
                     center = TRUE, retx = T)

summary(effect_pca_lasso)

pca_scores = as.data.frame(effect_pca_lasso$x)

correlations = cor(pca_lasso_data, pca_scores)
print(correlations)

pca_scores_scale = as.data.frame(scale(pca_scores))

#mean(pca_scores$PC1)
#sd(pca_scores$PC1)



## Set response variable
resp_pca_lasso <- scale(effect_data_clean$Effect_Size)

## Set predictor variables
pred_pca_lasso <- data.matrix(pca_scores_scale[, c("PC1",  "PC2",  "PC3", "PC4", "PC5", "PC6", "PC7",  "PC8", "PC9", "PC10", "PC11", "PC12")])



lasso = cv.glmnet(pred_pca_lasso, resp_pca_lasso, alpha = 1)

best_lambda <- lasso$lambda.min
best_lambda

plot(lasso)

best_lasso_model <- glmnet(pred_pca_lasso, resp_pca_lasso, alpha = 1, lambda = best_lambda)

coef(best_lasso_model)

resp_predict <- predict(best_lasso_model, s = best_lambda, newx = pred_pca_lasso)

sst <- sum((resp_pca_lasso - mean(resp_pca_lasso))^2)
sse <- sum((resp_predict - resp_pca_lasso)^2)

rsq = 1 - sse/sst

rsq


####

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

# fe_cube_new = fe_cube + 
#   ylab("") + 
#   theme(axis.title.x = element_text(size = 8))

atp_cube_new = atp_cube +
  ylab("") +
  theme(axis.title.x = element_text(size = 12))



moi_cube_new = grav_cube + 
  ylab("") + 
  theme(axis.title.x = element_text(size = 12))

fs_cube_new = fs_cube + 
  ylab("") + 
  theme(axis.title.x = element_text(size = 12))

#combine_scatter = ggarrange(fs_cube_new, atp_cube_new, moi_cube_new, common.legend = TRUE, nrow = 3, labels = c("B", "C", "D"), label.x = c(0.9, 0.9, 0.9), label.y = c(0.3, 0.3, 0.3), heights = c(1,1,1)) 

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Combined_Scatter_SFS.png"), width = 16, height = 8, units = "in", res = 300)

combine_scatter = ggarrange(fs_cube, grav_cube, common.legend = TRUE, nrow = 1, ncol = 2, heights = c(1,1)) 

 # annotate_scatter = annotate_figure(combine_scatter, left = text_grob("Cube Root Effect Size (mg/kg) (Wet - Dry)", x = 0.75, y = 0.5, size = 15, vjust = 1, rot = 90))
  
  combine_scatter

  dev.off()

combine_pca = ggarrange(pca_new, annotate_scatter, labels = c("A"), ncol = 2, label.x = c(0.08), label.y = c(0.95), widths = c(2, 1), heights = c(2,2)) 

combine_pca

dev.off()

## Cube PCA for LASSO####

# Fe outlier not in analysis
# Effect Size not in PCA

cube_effect_data_clean = cube_effect_data_corr %>% 
  #relocate(cube_Effect_Size, .before = cube_SpC_Diff) %>% 
  filter(Cube_Fe_mg_kg_Diff > -1)

#limits = c(0, 1450)
#ggplot(effect_data_clean, aes(x = FinGravMoi_Diff, y = Fine_Sand, color = Effect_Size)) +
#  geom_point() +
#  scale_color_gradient2(limits = limits, low = "firebrick2", mid = "goldenrod2",
 #                       high = "dodgerblue2", midpoint = (max(limits)+min(limits))/2)

##Run PCA without effect size
cube_data_pca = cube_effect_data_clean %>% 
  select(-c(Cube_Effect_Size))

cube_effect_pca <- prcomp(cube_data_pca, scale = TRUE,
                     center = TRUE, retx = T)

# Summary
summary(cube_effect_pca)

#cumulative variation

cumulative_var = cumsum(cube_effect_pca$sdev^2)/sum(cube_effect_pca$sdev^2)

pca_plot = data.frame(PC = 1:length(cube_effect_pca$sdev), Cumulative_Variance = cumulative_var)

ggplot(pca_plot, aes(x = PC, y = Cumulative_Variance)) +
  geom_col()


# Look at PCA loadings

loadings = cube_effect_pca$rotation

print(loadings)

loadings_df <- as.data.frame(loadings) %>% 
  select(c(PC1, PC2, PC3, PC4, PC5)) 

# Plotting the heatmap

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




set.seed(42)
## Set response variable
yvar <- data.matrix(scale(cube_effect_data_clean$Cube_Effect_Size, center = TRUE, scale = TRUE))
#yvar <- data.matrix(cube_effect_data_clean$Cube_Effect_Size)
mean(yvar)
sd(yvar)

## Set predictor variables
pca_scores = as.data.frame(scale(cube_effect_pca$x, center = TRUE, scale = TRUE))
#pca_scores = as.data.frame(cube_effect_pca$x)
mean(pca_scores$PC1)
sd(pca_scores$PC1)

xvars <- data.matrix(pca_scores[, c("PC1",  "PC2", "PC3", "PC4", "PC5")])

lasso = cv.glmnet(xvars, yvar, alpha = 1, 
                  standardize = FALSE, standardize.response = FALSE, intercept = FALSE, nfolds = 5)

best_lambda <- lasso$lambda.min
best_lambda

plot(lasso)

best_lasso_model <- glmnet(xvars, yvar, alpha = 1, lambda = best_lambda, standardize = FALSE, standardize.response = FALSE, intercept = FALSE)

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
  geom_point(aes(color = cube_effect_data_clean$Cube_Effect_Size), size = 3.5)+
  scale_color_gradient2(limits = cube_limits, low = "firebrick2", mid = "goldenrod2",
                        high = "dodgerblue2", midpoint = (max(cube_limits)+min(cube_limits))/2) +
  labs(color = paste0("Cubed Root Effect Size\n(Wet - Dry Rate)"))

cube_pca

dev.off()
