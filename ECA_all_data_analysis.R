#### Sensitivity Analysis For ECA removals ####

rm(list=ls());graphics.off()

#### Read in Data

#Individual samples 
all_data <- read.csv("C:/Github/ECA_Multireactor_Incubations/Data/Cleaned Data/All_ECA_Data.csv",header = TRUE) %>% 
  select(-c(X))

#Summary Data 

sum_data <- read.csv("C:/Github/ECA_Multireactor_Incubations/Data/Cleaned Data/Summary_ECA_Data.csv",header = TRUE) %>% 
  select(-c(X))

#Effect Size Data
effect_data <- read.csv("C:/Github/ECA_Multireactor_Incubations/Data/Cleaned Data/Effect_ECA_Data.csv",header = TRUE) %>% 
  select(-c(X))

#### Rates histograms ####

respiration_sep <- all_respiration %>%
  filter(Respiration_Rate_mg_DO_per_L_per_H != -9999) %>% 
  mutate(Respiration_Rate_mg_DO_per_L_per_H = abs(Respiration_Rate_mg_DO_per_L_per_H)) %>% 
  separate(Sample_Name, into = c("Sample", "rep"), sep = "-", remove = FALSE) %>% 
  mutate(Treat = case_when(grepl("W",rep)~"Wet",
                           grepl("D", rep) ~"Dry")) %>% 
  mutate(log_rate = (log10(Respiration_Rate_mg_DO_per_L_per_H + 1)))

respiration_sep$Treat <- as.factor(respiration_sep$Treat)

## Respiration Histograms ####

color_pallete <- colorRampPalette(colors = c("#D55E00", "#0072B2"))

num_colors <- nlevels(a$Treat)

samples_colors <- color_pallete(num_colors)

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/", as.character(Sys.Date()),"_Wet_Treatment_Histogram.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(subset(respiration_sep, Treat %in% "Wet"), aes(x = Respiration_Rate_mg_DO_per_L_per_H)) +
  geom_histogram(fill = "#0072B2")+
  ggtitle("Wet Rates")+
  xlab(expression("Respiration Rate (mg O"[2]*" L"^- 1*" H"^-1*")"))+
  theme(strip.text = element_text(
    size = 4))+
  ylim(0, 215)+
  theme_bw()

dev.off()


png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/", as.character(Sys.Date()),"_Dry_Treatment_Histogram.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(subset(respiration_sep, Treat %in% "Dry"), aes(x = Respiration_Rate_mg_DO_per_L_per_H)) +
  geom_histogram(fill = "#D55E00")+
  ggtitle("Dry Rates")+
  xlab(expression("Respiration Rate (mg O"[2]*" L"^- 1*" H"^-1*")"))+
  theme(strip.text = element_text(
    size = 4))+
  ylim(0, 215) + 
  theme_bw()

dev.off()

####

## Effect Size Histogram ####

effect_limits = c(-300, 300)

png(file = paste0("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/", as.character(Sys.Date()),"_effect_histogram.png"), width = 10, height = 10, units = "in", res = 300)

ggplot(effect_all, aes(x = Effect_Size))+
  # geom_histogram(binwidth = 0.15, fill = "#009E73")+
  geom_histogram(binwidth = 15, aes(fill = after_stat(x))) +
  scale_fill_gradient2(name = "Effect Size", limits = effect_limits, low = "firebrick2", mid = "goldenrod2",
                       high = "dodgerblue2", midpoint = (max(effect_limits)+min(effect_limits))/2) +
  theme_bw()+
  theme(axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size =18))+
  xlim(c(-300,300))+
  ylab("Count\n")+
  xlab(expression("\n Effect Size (Wet - Dry Rate; mg O"^2*" L"^-1*" H"^-1*")"))


dev.off()

log_effect_limits <- c(-2.48, 2.48)

png(file = paste0("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/", as.character(Sys.Date()),"_log_effect_histogram.png"), width = 10, height = 10, units = "in", res = 300)

ggplot(effect_all, aes(x = Log_Effect_Size))+
  # geom_histogram(binwidth = 0.15, fill = "#009E73")+
  geom_histogram(binwidth = 0.15, aes(fill = after_stat(x))) +
  scale_fill_gradient2(name = "Log Effect Size", limits = log_effect_limits, low = "firebrick2", mid = "goldenrod2",
                       high = "dodgerblue2", midpoint = (max(log_effect_limits)+min(log_effect_limits))/2) +
  theme_bw()+
  theme(axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size =18))+
  xlim(c(-2.48, 2.48))+
  ylab("Count\n")+
  xlab(expression("\n Log Effect Size (Wet - Dry Rate; mg O"^2*" L"^-1*" H"^-1*")"))

dev.off()

####


##D50/RUSLE Calculations - Need to be updated ####
grn_all <- grn %>% 
  dplyr::select(-c(geom_rusle, geom, Percent_Tot_Sand, Percent_Mud)) %>% 
  mutate(Percent_Coarse_Sand_Finer = Percent_Coarse_Sand + Percent_Med_Sand + Percent_Fine_Sand + Percent_Silt + Percent_Clay) %>% 
  mutate(Percent_Med_Sand_Finer = Percent_Med_Sand + Percent_Fine_Sand + Percent_Silt + Percent_Clay) %>% 
  mutate(Percent_Fine_Sand_Finer = Percent_Fine_Sand + Percent_Silt + Percent_Clay) %>% 
  mutate(Percent_Silt_Finer = Percent_Silt + Percent_Clay) %>% 
  mutate(Percent_Clay_Finer = Percent_Clay) %>% 
  dplyr::select(-c(Percent_Coarse_Sand, Percent_Med_Sand, Percent_Fine_Sand, Percent_Silt, Percent_Clay)) %>% 
  pivot_longer(cols = c("Percent_Coarse_Sand_Finer", "Percent_Med_Sand_Finer", "Percent_Fine_Sand_Finer", "Percent_Silt_Finer", "Percent_Clay_Finer"), names_to = "Category", values_to = "fraction") %>% 
  mutate(size = case_when(
    grepl("Coarse", Category) ~ 1.25,
    grepl("Med", Category) ~ 0.375,
    grepl("Fine_Sand", Category) ~ 0.1,
    grepl("Silt", Category) ~ 0.026, 
    grepl("Clay", Category) ~ 0.01)) %>% 
  mutate(y_int = 50) %>% 
  group_by(kit) %>% 
  mutate(slope = if_else(fraction[Category == "Percent_Med_Sand_Finer"] < 50, ((fraction[Category == "Percent_Coarse_Sand_Finer"] - fraction[Category == "Percent_Med_Sand_Finer"]) / (size[Category == "Percent_Coarse_Sand_Finer"] - size[Category == "Percent_Med_Sand_Finer"])), 
                if_else(fraction[Category == "Percent_Fine_Sand_Finer"] < 50, ((fraction[Category == "Percent_Med_Sand_Finer"] - fraction[Category == "Percent_Fine_Sand_Finer"]) / (size[Category == "Percent_Med_Sand_Finer"] - size[Category == "Percent_Fine_Sand_Finer"])), 
                if_else(fraction[Category == "Percent_Silt_Finer"] < 50, ((fraction[Category == "Percent_Fine_Sand_Finer"] - fraction[Category == "Percent_Silt_Finer"]) / (size[Category == "Percent_Fine_Sand_Finer"] - size[Category == "Percent_Silt_Finer"])), ((fraction[Category == "Percent_Silt_Finer"] - fraction[Category == "Percent_Clay_Finer"]) / (size[Category == "Percent_Silt_Finer"] - size[Category == "Percent_Clay_Finer"]))
    )))) %>% 
  mutate(x_int = if_else(fraction[Category == "Percent_Med_Sand_Finer"] < 50, (((50 - fraction[Category == "Percent_Med_Sand_Finer"])/slope[Category == "Percent_Coarse_Sand_Finer"])+size[Category == "Percent_Med_Sand_Finer"]),
                  if_else(fraction[Category == "Percent_Fine_Sand_Finer"] < 50, (((50 - fraction[Category == "Percent_Fine_Sand_Finer"])/slope[Category == "Percent_Med_Sand_Finer"])+size[Category == "Percent_Fine_Sand_Finer"]),
                  if_else(fraction[Category == "Percent_Silt_Finer"] < 50,(((50 - fraction[Category == "Percent_Silt_Finer"])/slope[Category == "Percent_Fine_Sand_Finer"])+size[Category == "Percent_Silt_Finer"]), (((50 - fraction[Category == "Percent_Clay_Finer"])/slope[Category == "Percent_Silt_Finer"])+size[Category == "Percent_Clay_Finer"])))))
  

d50 <- grn_all %>% 
  distinct(kit, .keep_all = TRUE) %>% 
  dplyr::select(c(kit, x_int)) %>% 
  rename(d50 = x_int)

grn <- merge(d50, grn)

grn <- grn %>% 
  relocate(d50, .after = geom)

# geom <- grn %>% 
#   select(c(kit, geom_rusle,geom)) %>% 
#   rename(log_geom_mean = geom_rusle) %>% 
#   rename(geom_mean = geom)
# 
# all_geom <- left_join(d50, geom, by = c("kit"))
# 
# geom_vs_d50 <- ggplot(all_geom, aes(x = geom_mean, y = d50))+
#   geom_point()+
#   geom_smooth()
# 
# geom_vs_rusle <- ggplot(all_geom, aes(x = geom_mean, y = log_geom_mean))+
#   geom_point()+
#   geom_smooth()
# 
# d50_vs_rusle <- ggplot(all_geom, aes(x = d50, y = log_geom_mean))+
#   geom_point()+
#   geom_smooth()


#Coarse sand = 0.5 - 2 mm (1.25 mm), med sand = 0.25 - 0.499 (0.375) , fine sand = 0.05 - 0.25 (0.1), silt = , clay =

## Cleaned Data - Needs to be updated ####
all_samples_clean <- all_samples_grn %>% 
  na.omit()  %>% 
  #filter(!is.na(Sample_ID)) %>% 
  remove_rownames %>% 
  column_to_rownames(var = "Sample_Name") %>% 
  rename(`Rate (mg/L/H)` = Respiration_Rate_mg_DO_per_L_per_H) %>% 
  rename(`Initial Gravimetric Water` = grav_initial) %>% 
  rename(`Final Gravimetric Water` = grav_final) %>% 
  rename(`Lost Gravimetric Water` = lost_grav_perc) %>% 
 # rename(`Fe (II) (mg/kg)`  = Mean_Rep_Fe_mg_kg) %>% 
  rename(`Fe (II) (mg/L)` = Mean_Rep_Fe_mg_per_L) %>% 
  rename(`Sp. Conducitivity` = SpC) %>% 
  rename(`pH` = pH) %>% 
  rename(`% Fine Sand` = Percent_Fine_Sand) %>% 
  rename(`% Med. Sand` = Percent_Med_Sand) %>% 
  rename(`% Coarse Sand` = Percent_Coarse_Sand) %>% 
  rename(`% Mud` = Percent_Mud) %>% 
  rename(`Geometric Mean` = geom) %>% 
  rename(`RUSLE Geometric Mean` = geom_rusle) %>% 
  rename(`D50` = d50)


## Respiration vs. Mud ####
png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/ESS-PI_EGU/", as.character(Sys.Date()),"_Respiration_vs_Mud_lowess.png"), width = 8, height = 8, units = "in", res = 300)

par(mar = c(5 ,6 , 6, 6), xpd = FALSE)

plot(all_samples_clean$`% Mud`, all_samples_clean$`Rate (mg/L)`, cex= 1.8, cex.lab = 2.2, cex.axis = 1.8 ,xlab = "% Mud" , ylab = expression("Respiration Rate (mg O"[2]*" L"^- 1*" min"^-1*")"), lwd = 2, col = samples_colors[all_samples_clean$Treat], pch = c(16,17)[as.numeric(all_samples_clean$Treat)])

#lines(lowess(all_samples_clean$`% Mud`[all_samples_clean$Treat=="Dry"], all_samples_clean$`Rate (mg/L)`[all_samples_clean$Treat=="Dry"]), col = "#D55E00", lwd = 3)


#lines(lowess(all_samples_clean$`% Mud`[all_samples_clean$Treat=="Wet"], all_samples_clean$`Rate (mg/L)`[all_samples_clean$Treat=="Wet"]), col = "#0072B2", lwd = 3)

par(xpd = TRUE)

legend(x = "topright", inset = c(-0.15,0.45), legend = paste(levels(all_samples_clean$Treat)), col = samples_colors, pch = c(16,17), cex = 1)

dev.off()

## Fe (mg/kg) vs. Respiration ####

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/ESS-PI_EGU/", as.character(Sys.Date()),"_Respiration_vs_Fe_mg_kg_lowess.png"), width = 8, height = 8, units = "in", res = 300)

par(mar = c(5 ,6 , 6, 6), xpd = FALSE)

plot(all_samples_clean$`Fe (II) (mg/L)`, all_samples_clean$`Rate (mg/L)`, cex= 1.8, cex.lab = 2.2, cex.axis = 1.8 ,xlab = expression("Fe (II) (mg kg"^-1*")"), ylab = expression("Respiration Rate (mg O"[2]*" L"^- 1*" min"^-1*")"), lwd = 2, col = samples_colors[all_samples_clean$Treat], pch = c(16,17)[as.numeric(all_samples_clean$Treat)])

#lines(lowess(all_samples_clean$`Mean_Fe_mg_per_L`, all_samples_clean$`Rate (mg/L)`), col = "blue", lwd = 3)

par(xpd = TRUE)

legend(x = "topright", inset = c(-0.15,0.45), legend = paste(levels(all_samples_clean$Treat)), col = samples_colors,  pch = c(16,17), cex = 1)

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



## Fine Sand, Mud, Coarse Sand, Fe Difference vs. Effect Size ####

fe <- ggplot(log_effect, aes(x = `Fe Difference (mg/L)`, y = `Effect Size`)) +
  geom_point()+
  geom_smooth()+
  xlab("Log Fe Difference (mg/L)")+
  ylab("Log Effect Size")

fine <- ggplot(log_effect, aes(x = `% Fine Sand`, y = `Effect Size`)) +
  geom_point()+
  geom_smooth()+
  xlab("Log % Fine Sand")+
  ylab("Log Effect Size")+
  ylim(0, 0.8)

coarse <- ggplot(log_effect, aes(x = `% Coarse Sand`, y = `Effect Size`)) +
  geom_point()+
  geom_smooth()+
  xlab("Log % Coarse Sand")+
  ylab("Log Effect Size")

mud <- ggplot(log_effect, aes(x = `% Mud`, y = `Effect Size`)) +
  geom_point()+
  geom_smooth()+
  xlab("Log % Mud")+
  ylab("Log Effect Size")+
  ylim(0, 0.75)

all <- plot_grid(fe, fine, coarse, mud)

ggsave("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/all_regression.png", all, width = 8, height = 8)

## Fine Sand vs. Effect Size ####

fine_sand_hist <- ggplot(effect, aes(x = `% Fine Sand`))+
 geom_histogram()

effect_hist <- ggplot(effect, aes(x = `Effect Size`))+
  geom_histogram()

mud_hist <- ggplot(effect, aes(x = `% Mud`))+
  geom_histogram()

all_hist <- plot_grid(fine_sand_hist, mud_hist, effect_hist)

ggsave("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/2023-07-07_texture_hist.png", all_hist, width = 8, height = 8)

##this constant (+1) not enough to correct SpC/Fe differences

log_effect <- log10(effect +1)

log_fine_sand_hist <- ggplot(log_effect, aes(x = `% Fine Sand`))+
  geom_histogram()+
  xlab("Log % Fine Sand")

log_effect_hist <- ggplot(log_effect, aes(x = `Effect Size`))+
  geom_histogram()+
  xlab("Log Effect Size")

log_mud_hist <- ggplot(log_effect, aes(x = `% Mud`))+
  geom_histogram() + 
  xlab("Log % Mud")

all_log_hist <- plot_grid( log_fine_sand_hist, log_effect_hist, log_mud_hist)

ggsave("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/2023-07-07_log_texture_hist.png", all_log_hist, width = 8, height = 8)

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/", as.character(Sys.Date()),"_Comparison_Effect_vs_Fine_Sand_Scatter.png"), width = 16, height = 8, units = "in", res = 300)

par(mfrow = c(1, 2))
par(mar = c(5 ,6 , 4, 1))

log_comp <- plot(log_effect$`% Fine Sand`, log_effect$`Effect Size`, cex= 1.8, cex.lab = 1.5, cex.axis = 1.8 ,xlab = "Log % Fine Sand" , ylab = expression(paste("Log Effect Size (Wet - Dry Rate) (mg O"[2]*" L"^-1*" min"^-1*")")), lwd = 2)

log_comp <- log_comp + lines(lowess(log_effect$`% Fine Sand`, log_effect$`Effect Size`), col = "blue", lwd = 3)

comp <- plot(effect$`% Fine Sand`, effect$`Effect Size`, cex= 1.8, cex.lab = 1.5, cex.axis = 1.8 ,xlab = "% Fine Sand" , ylab = expression(paste("Effect Size (Wet - Dry Rate) (mg O"[2]*" L"^-1*" min"^-1*")")), lwd = 2)

comp <- comp + lines(lowess(effect$`% Fine Sand`, effect$`Effect Size`), col = "blue", lwd = 3)

dev.off()


ggsave("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/2023-07-07_log_texture_hist.png", all_log_hist, width = 8, height = 8)

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/", as.character(Sys.Date()),"_Effect_vs_Fine_Sand_Scatter.png"), width = 8, height = 8, units = "in", res = 300)

par(mar = c(5 ,6 , 4, 1))

plot(effect$`% Fine Sand`, effect$`Effect Size`, cex= 1.8, cex.lab = 1.8, cex.axis = 1.8 ,xlab = "% Fine Sand" , ylab = expression(paste("Effect Size (Wet - Dry Rate) (mg O"[2]*" L"^-1*" min"^-1*")")), lwd = 2)

#lines(lowess(effect$`% Fine Sand`, effect$`Effect Size`), col = "blue", lwd = 3)

dev.off()


png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/", as.character(Sys.Date()),"_Effect_vs_Fine_Sand_Scatter.png"), width = 8, height = 8, units = "in", res = 300)

# ggplot(effect, aes(x = `% Fine Sand`, y = `Effect Size`))+
#   geom_point()+
#   #geom_smooth(method = lm)+
#   theme_bw()

par(mar = c(5 ,6 , 4, 1))

plot(effect$`% Fine Sand`, effect$`Effect Size`, cex= 1.8, cex.lab = 1.8, cex.axis = 1.8 ,xlab = "Log % Fine Sand" , ylab = expression(paste("Log Effect Size (Wet - Dry Rate) (mg O"[2]*" L"^-1*" min"^-1*")")), lwd = 2)

lines(lowess(effect$`% Fine Sand`, effect$`Effect Size`), col = "blue", lwd = 3)

dev.off()

## Mud vs. Effect Size ####

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/", as.character(Sys.Date()),"_Effect_vs_Mud.png"), width = 8, height = 8, units = "in", res = 300)

par(mar = c(5 ,6 , 4, 1))

plot(effect$`% Mud`, effect$`Effect Size`, cex= 1.8, cex.lab = 1.8, cex.axis = 1.8 ,xlab = "% Mud" , ylab = expression(paste("Effect Size (Wet - Dry Rate) (mg O"[2]*" L"^-1*" min"^-1*")")), lwd = 2)

#lines(lowess(effect$`% Mud`, effect$`Effect Size`), col = "blue", lwd = 3)

dev.off()


png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/", as.character(Sys.Date()),"_Log_Effect_vs_Log_Mud.png"), width = 8, height = 8, units = "in", res = 300)

par(mar = c(5 ,6 , 4, 1))

plot(log_effect$`% Mud`, log_effect$`Effect Size`, cex= 1.8, cex.lab = 1.8, cex.axis = 1.8 ,xlab = "Log % Mud" , ylab = expression(paste("Log Effect Size (Wet - Dry Rate) (mg O"[2]*" L"^-1*" min"^-1*")")), lwd = 2)

#lines(lowess(log_effect$`% Mud`, log_effect$`Effect Size`), col = "blue", lwd = 3)

dev.off()

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/", as.character(Sys.Date()),"_Comparison_Effect_vs_Mud_Scatter.png"), width = 16, height = 8, units = "in", res = 300)

par(mfrow = c(1, 2))
par(mar = c(5 ,6 , 4, 1))

mud <- plot(effect$`% Mud`, effect$`Effect Size`, cex= 1.8, cex.lab = 1.8, cex.axis = 1.8 ,xlab = "% Mud" , ylab = expression(paste("Effect Size (Wet - Dry Rate) (mg O"[2]*" L"^-1*" min"^-1*")")), lwd = 2)

mud <- mud + lines(lowess(effect$`% Mud`, effect$`Effect Size`), col = "blue", lwd = 3)

par(mar = c(5 ,6 , 4, 1))

log_mud <- plot(log_effect$`% Mud`, log_effect$`Effect Size`, cex= 1.8, cex.lab = 1.8, cex.axis = 1.8 ,xlab = "Log % Mud" , ylab = expression(paste("Log Effect Size (Wet - Dry Rate) (mg O"[2]*" L"^-1*" min"^-1*")")), lwd = 2)

log_mud <- log_mud + lines(lowess(log_effect$`% Mud`, log_effect$`Effect Size`), col = "blue", lwd = 3)

dev.off()

# James Correlation Matrix ####

panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r = (cor(x, y))
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

png(file = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/", as.character(Sys.Date()),"_James_Effect_Correlation_Matrix.png"), width = 12, height = 12, units = "in", res = 300)

pairs(log_effect, lower.panel = panel.smooth,upper.panel = panel.cor, gap = 0, cex.labels = 1, cex = 1)

dev.off()

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
