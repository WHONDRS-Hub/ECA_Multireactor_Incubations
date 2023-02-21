###### Load Library ######

library(dplyr); library(ggplot2);library(ggsignif)
library(ggpubr);library(reshape2);library(ggpmisc)
library(segmented);library(broom)
library(ggpmisc);library(segmented);library(lubridate); library(readxl);
library(tidyverse);library(patchwork)

##### Load data ######
rm(list=ls());graphics.off()

# Set working directory to data file
#Example:
pnnl.user = 'laan208'

input.path = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/")

setwd(input.path)

#path for raw data
path <- ("rates/")

#path for reading in 100% saturation values for each kit based on pressure/temperature during disk calibration

fast.sat <- paste0("C:/Users/", pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/rates/fast_rate_calculations.xlsx")

##what threshold for bad first point?

fast.rates <- read_excel(fast.sat)  
  
#kits that go to 0 quickly - 23W, 23D, 27W, 40W, 32D2, 32D4, 12W, 11W,  34W, 35W, 35D1245, 69W

fast.rates.kits <- fast.rates %>% 
  rename("DO_mg_L" = "DO_sat_mg_L") %>% 
  filter(grepl("EC_11_INC-W|EC_12_INC-W|EC_23_INC|EC_27_INC-W|EC_32_INC-W|EC_32_INC-D2|EC_32_INC-D4|EC_34_INC-W|EC_35_INC-W|EC_35_INC-D1|EC_35_INC-D2|EC_35_INC-D4|EC_35_INC-D5|EC_40_INC-W|EC_51_INC-W", Sample_ID))

#read in respiration data and clean
import_data = function(path){
  
  # pull a list of file names in the target folder with the target pattern
  # then read all files and combine
  
  filePaths <- list.files(path = path, recursive = T, pattern = "\\.csv$", full.names = TRUE)
  
  filePaths <- filePaths[!grepl("results_linear_fit|elapsed.time.ratios|ECA_Sediment", filePaths)]
  
  #filePaths <- files[!grepl("elapsed.time.ratios", files)]
  
  # dat <- 
  do.call(rbind, lapply(filePaths, function(path){
    # then add a new column `source` to denote the file name
    df <- read.csv(path, skip = 4)
    df[["source_file"]] <- rep(path, nrow(df)) # add a column for the source file

    df %>%
      na.omit () %>% 
      as_tibble(row.names = 1:nrow(df)) %>% 
        # since all rows are in 2-min increments, just multiply the row number by 2
       tibble::rownames_to_column() %>%
       mutate(elapsed_min = as.numeric(rowname)*2) %>%
       dplyr::select(-rowname) %>% 
      # make longer, so all the data are in a single column
      pivot_longer(-c(elapsed_min,source_file), names_to = "disc_number", values_to = "DO_mg_L", values_transform = as.numeric) %>% 
      # remove unnecessary strings from the source_file name
      mutate(source_file = str_remove_all(source_file, paste0(path, "/")))
  }
  ))
}
data = import_data(path)


data_long = 
  data %>% 
  mutate(disc_number = str_remove_all(disc_number, "X")) %>%
  mutate(source_file = str_remove_all(source_file, "optode data/"),source_file = str_remove_all(source_file, ".csv")) %>%  
  filter(elapsed_min > 0) %>% 
  separate(col = source_file, into = c("fol1", "fol2", "fol3","fol4", "source_file"), sep = "/") %>% 
  dplyr::select(-fol1, -fol2, -fol3,-fol4)

vials <- data_long %>% 
  filter(source_file != "results _03a") %>% 
  filter(source_file != "results_03b") %>% 
  filter(source_file != "results_01") %>% 
  filter(source_file != "results_02") %>% 
  filter(source_file != "results_04_08") %>% 
  filter(source_file != "results_10_15") %>% 
  filter(source_file != "results_06_07") %>% 
  na.omit(vials)

#import mapping files
import_data = function(input.path){ 
  
  #pull a list of files in target folder with correct pattern
  #read all files and combine
  
  map.file <-  list.files(input.path, recursive = T, pattern = "\\.xlsx$",full.names = T)
  
  map.file <- map.file[!grepl("red|Red|EC_01_|EC_02_|EC_03_|EC_06_07|EC_10_15|EC_04_08|fast|combined|salinity", map.file)]
  
  mapping <- lapply(map.file, read_xlsx)
  
  for (i in 1:length(mapping)){mapping[[i]] <- cbind(mapping[[i]], map.file[i])}
  
  all.map <- 
    do.call(rbind,mapping)
}

all.map = import_data(input.path)

all.map = all.map %>% 
  rename("source_file" = "map.file[i]") %>% 
  rename("disc_number" = "Disk_ID") %>% 
  mutate(source_file = str_remove_all(source_file, paste0(input.path, "/"))) %>% 
  separate(source_file, sep = "/", c("source_file", "file")) %>% 
  mutate(source_file = str_replace(source_file, "EC", "results")) %>% 
  dplyr::select(-`Disk Calibration date`, -Location, -`Time on`, -`Time off`, -SpC, -pH, -Temp, -Notes, -file)

all.samples <- merge(vials, all.map)

bind <- merge(all.samples, fast.rates.kits, all = TRUE) %>% 
  dplyr::select(-`calibration date`) 

#fill in samples, remove last two columns for cleaner data frame

#generate another dataset (w/ time, DO) with everything that has been removed: high, low (median), at the end, then can plot everything on top of each other with different colors to see what has been removed

#still need to have it populate into respiration how many points that it removed

min.points = 2
threshold = 0.99
res.threshold = 0.01
slope.thresh = 0
do.thresh = 0.5
ymax = max(na.omit(bind$DO_mg_L))
ymin = min(na.omit(bind$DO_mg_L))

respiration <- as.data.frame(matrix(NA, ncol = 13, nrow =1))

colnames(respiration) = c("Sample_ID", "slope_of_the_regression", "rate_mg_per_L_per_min", "rate_mg_per_L_per_h", "R_squared", "R_squared_adj","residuals", "p_value", "total_incubation_time_min", "number_of_points", "removed_points_high", "removed_points_beg", "removed_points_end")

location = c("-W1", "-W2", "-W3", "-W4", "-W5", "-D1", "-D2", "-D3", "-D4", "-D5")

rate = as.data.frame(matrix(NA, ncol = 9, nrow = length(unique(bind$Sample_ID))))

colnames(rate) = c("Sample_ID", "slope_of_the_regression", "rate_mg_per_L_per_min", "rate_mg_per_L_per_h", "R_squared", "R_squared_adj","residuals", "p_value", "total_incubation_time_min")


#bad kits, fix later
bind <- bind %>% 
  filter(Sample_ID != "EC_27_INC-W5") %>% 
  filter(Sample_ID != "EC_27_INC-D1") %>% 
  filter(Sample_ID != "EC_14_INC-W5") %>% 
  filter(Sample_ID != "EC_13_INC-W5") %>% 
  filter(Sample_ID != "EC_13_INC-D4") %>% 
  filter(Sample_ID != "EC_27_INC-D2") %>% 
  filter(!(elapsed_min < 8 & Sample_ID == "EC_14_INC-W1")) %>% 
  filter(!(elapsed_min <= 2 & Sample_ID == "EC_05_INC-D4")) %>% 
  filter(!(elapsed_min <= 2 & Sample_ID == "EC_05_INC-D5")) %>% 
  filter(!(elapsed_min <= 2 & Sample_ID == "EC_32_INC-D3")) %>% 
  filter(!(elapsed_min <= 2 & Sample_ID == "EC_56_INC-W3")) %>% 
  filter(Sample_ID != "EC_72_INC-W5") %>% 
  filter(Sample_ID != "EC_72_INC-D5")
  

for (i in 1:length(location)){
  
  data_location_subset = bind[grep(location[i],bind$Sample_ID),]
  
  unique.incubations = unique(data_location_subset$Sample_ID)
  
  rate = as.data.frame(matrix(NA, ncol = 9, nrow = length(unique(bind$Sample_ID))))
  
  colnames(rate) = c("Sample_ID", "slope_of_the_regression", "rate_mg_per_L_per_min", "rate_mg_per_L_per_h", "R_squared", "R_squared_adj","residuals", "p_value", "total_incubation_time_min")
  
  for (j in 1:length(unique.incubations)){
    
    data_site_subset = subset(data_location_subset, data_location_subset$Sample_ID == unique.incubations[j])
    data_site_subset = data_site_subset[order(data_site_subset$elapsed_min, decreasing = FALSE),]
    data_site_subset$elapsed_min = as.numeric(data_site_subset$elapsed_min)
    
    data_site_subset_low = data_site_subset %>% 
      filter(DO_mg_L < 14) 
    
    data_site_subset_beg = data_site_subset_low
    
    for (k in 1:3) {
      
      if(data_site_subset_beg$DO_mg_L[k] < median(data_site_subset_beg$DO_mg_L)) {
        
        data_site_subset_beg = data_site_subset_beg[-k,]
        
      }
      else if (data_site_subset_beg$DO_mg_L[k] >= median(data_site_subset_beg$DO_mg_L)){
        #break()
        data_site_subset_beg = data_site_subset_beg
        #ask James if we want to assume that if the first point is fine, no other points need to be checked
        
      }
      
    }
    
    ##remove samples if at >4 minutes, they are below the DO threshold. This is trying to remove low values from the back end of curves
    
    data_site_subset_thresh = data_site_subset_beg  %>% 
      filter(!(elapsed_min > 4 & DO_mg_L < do.thresh))
    
    data_site_subset_fin = data_site_subset_thresh
    
    fit = lm(data_site_subset_fin$DO_mg_L~data_site_subset_fin$elapsed_min)
    u = fit$coefficients
    b = u[[1]]
    c = u[[2]]
    r = summary(fit)$r.squared
    residuals = deviance(fit)
    r.adj = summary(fit)$adj.r.squared
    p = summary(fit)$coefficients[4]
    rog = r
    resog = residuals
    
    temp = data_site_subset_fin[-nrow(data_site_subset_fin),]
    fit.temp = lm(temp$DO_mg_L~temp$elapsed_min)
    rtemp = summary(fit.temp)$r.squared
    restemp = deviance(fit.temp)
    
    #points being removed without r2 increasing with removal
    # for(l in 1:60){
    # 
    #   if (r < threshold & nrow(data_site_subset_fin)>=min.points){
    # 
    #     data_site_subset_fin = data_site_subset_fin[-nrow(data_site_subset_fin),]
    # 
    #     fit = lm(data_site_subset_fin$DO_mg_L~data_site_subset_fin$elapsed_min)
    #     u = fit$coefficients
    #     b = u[[1]] #Intercept
    #     c = u[[2]] #rate mg/L min
    #     r = summary(fit)$r.squared
    #     r.adj = summary(fit)$adj.r.squared
    #     p = summary(fit)$coefficients[4]
    #     r2 = r
    #   }
    # }
    
    

    #points being removed with increase in r2 after removal
     # if (((r < threshold & rtemp >= rog) | (residuals > res.threshold & restemp <= residuals)) & nrow(data_site_subset_fin) >= min.points){
     # 
     #  data_site_subset_fin = data_site_subset_fin[-nrow(data_site_subset_fin),]
     #  fit = lm(data_site_subset_fin$DO_mg_L~data_site_subset_fin$elapsed_min)
     #  u = fit$coefficients
     #  b = u[[1]] #Intercept
     #  c = u[[2]] #rate mg/L min
     #  r = summary(fit)$r.squared
     #  residuals = deviance(fit)
     #  r.adj = summary(fit)$adj.r.squared
     #  p = summary(fit)$coefficients[4]
     #  r2 = r
     #  res2 = residuals
     # 
     # 
     #  if (((r < threshold & r >= rog) | (residuals > res.threshold & residuals <= resog)) & nrow(data_site_subset_fin)>=min.points){
     #    data_site_subset_fin = data_site_subset_fin[-nrow(data_site_subset_fin),]
     #    fit = lm(data_site_subset_fin$DO_mg_L~data_site_subset_fin$elapsed_min)
     #    u = fit$coefficients
     #    b = u[[1]] #Intercept
     #    c = u[[2]] #rate mg/L min
     #    r = summary(fit)$r.squared
     #    residuals = deviance(fit)
     #    r.adj = summary(fit)$adj.r.squared
     #    p = summary(fit)$coefficients[4]
     #    r3 = r
     #    res3 = residuals
     #  }
     #  if (((r < threshold&r >= r2) | (residuals > res.threshold & residuals >=res2)) & nrow(data_site_subset_fin)>=min.points){
     #    data_site_subset_fin = data_site_subset_fin[-nrow(data_site_subset_fin),]
     #    fit = lm(data_site_subset_fin$DO_mg_L~data_site_subset_fin$elapsed_min)
     #    u = fit$coefficients
     #    b = u[[1]] #Intercept
     #    c = u[[2]] #rate mg/L min
     #    r = summary(fit)$r.squared
     #    r.adj = summary(fit)$adj.r.squared
     #    p = summary(fit)$coefficients[4]
     #  }
     # }
    
 
    
    # my.formula <- y ~ x
    # final <- ggplot(data_site_subset_fin, aes(x = elapsed_min, y = DO_mg_L)) + coord_cartesian(ylim = c(0,15))+ geom_point(size = 2) + expand_limits(x = 0, y = 0) +
    #   geom_smooth(method = "lm", se=F, formula = my.formula) +
    #   stat_poly_eq(formula = my.formula,label.y = "top",label.x = "right", aes(label = paste( ..rr.label.., sep = "~~~"),size=1), parse = TRUE)+stat_fit_glance(data=data_site_subset_fin, method = 'lm', method.args = list(formula = my.formula),geom = 'text',aes(label =paste("          p = ",signif(..p.value.., digits = 1), sep = ""),size=1),label.y = c(14.25),label.x = "left") +
    #   theme_bw()+theme(legend.title = element_blank(),legend.background = element_rect(fill = 'NA'), legend.text = element_text(size = 12,face="bold"))+
    #   labs(y = expression(Dissolved_Oxygen_mg_per_L), x = expression(Time_Elapsed_min))+ theme(axis.text.x=element_text(size = 12,face="bold"))+
    #   ggtitle("Final " ,data_site_subset_fin$Sample_ID[1]) +
    #   theme(plot.title = element_text(lineheight=.8, face="bold"))+
    #   theme(axis.text.x=element_text(colour = c("black","black")))+
    #   theme(aspect.ratio=1)+
    #   theme(axis.text.y=element_text(size = 12,face="bold"))+
    #   theme(axis.title.x =element_text(size = 12,face="bold"))+
    #   theme(axis.title =element_text(size = 12,face="bold"))+
    #   theme(axis.title.y =element_text(size = 12,face="bold"))
    # 
    # high_rem <- ggplot(data_site_subset_low, aes(x = elapsed_min, y = DO_mg_L)) + coord_cartesian(ylim = c(0,15))+ geom_point(size = 2) + expand_limits(x = 0, y = 0) +
    #   geom_smooth(method = "lm", se=F, formula = my.formula) +
    #   stat_poly_eq(formula = my.formula,label.y = "top",label.x = "right", aes(label = paste( ..rr.label.., sep = "~~~"),size=1), parse = TRUE)+stat_fit_glance(data=data_site_subset_low, method = 'lm', method.args = list(formula = my.formula),geom = 'text',aes(label =paste("         p = ",signif(..p.value.., digits = 1), sep = ""),size=1),label.y = c(14.25),label.x = "left") +
    #   theme_bw()+theme(legend.title = element_blank(),legend.background = element_rect(fill = 'NA'), legend.text = element_text(size = 12,face="bold"))+
    #   labs(y = expression(Dissolved_Oxygen_mg_per_L), x = expression(Time_Elapsed_min))+ theme(axis.text.x=element_text(size = 12,face="bold"))+
    #   ggtitle("High Points Removed " ,data_site_subset_low$Sample_ID[1]) +
    #   theme(plot.title = element_text(lineheight=.8, face="bold"))+
    #   theme(axis.text.x=element_text(colour = c("black","black")))+
    #   theme(aspect.ratio=1)+
    #   theme(axis.text.y=element_text(size = 12,face="bold"))+
    #   theme(axis.title.x =element_text(size = 12,face="bold"))+
    #   theme(axis.title =element_text(size = 12,face="bold"))+
    #   theme(axis.title.y =element_text(size = 12,face="bold"))
    # 
    # beg_rem <- ggplot(data_site_subset_beg, aes(x = elapsed_min, y = DO_mg_L)) + coord_cartesian(ylim = c(0,15))+ geom_point(size = 2) + expand_limits(x = 0, y = 0) +
    #   geom_smooth(method = "lm", se=F, formula = my.formula) +
    #   stat_poly_eq(formula = my.formula,label.y = "top",label.x = "right", aes(label = paste( ..rr.label.., sep = "~~~"),size=1), parse = TRUE)+stat_fit_glance(data=data_site_subset_beg, method = 'lm', method.args = list(formula = my.formula),geom = 'text',aes(label =paste("         p = ",signif(..p.value.., digits = 1), sep = ""),size=1),label.y = c(14.25),label.x = "left") +
    #   theme_bw()+theme(legend.title = element_blank(),legend.background = element_rect(fill = 'NA'), legend.text = element_text(size = 12,face="bold"))+
    #   labs(y = expression(Dissolved_Oxygen_mg_per_L), x = expression(Time_Elapsed_min))+ theme(axis.text.x=element_text(size = 12,face="bold"))+
    #   ggtitle("Beginning Points Removed ", data_site_subset_beg$Sample_ID[1]) +
    #   theme(plot.title = element_text(lineheight=.8, face="bold"))+
    #   theme(axis.text.x=element_text(colour = c("black","black")))+
    #   theme(aspect.ratio=1)+
    #   theme(axis.text.y=element_text(size = 12,face="bold"))+
    #   theme(axis.title.x =element_text(size = 12,face="bold"))+
    #   theme(axis.title =element_text(size = 12,face="bold"))+
    #   theme(axis.title.y =element_text(size = 12,face="bold"))
    # 
    # all <- ggplot(data_site_subset, aes(x = elapsed_min, y = DO_mg_L)) + coord_cartesian(ylim = c(0,15))+ geom_point(size = 2) + expand_limits(x = 0, y = 0) +
    #   geom_smooth(method = "lm", se=F, formula = my.formula) +
    #   stat_poly_eq(formula = my.formula,label.y = "top",label.x = "right", aes(label = paste( ..rr.label.., sep = "~~~"),size=1), parse = TRUE)+stat_fit_glance(data=data_site_subset, method = 'lm', method.args = list(formula = my.formula),geom = 'text',aes(label =paste("         p = ",signif(..p.value.., digits = 1), sep = ""),size=1),label.y = c(14.25),label.x = "left") +
    #   theme_bw()+theme(legend.title = element_blank(),legend.background = element_rect(fill = 'NA'), legend.text = element_text(size = 12,face="bold"))+
    #   labs(y = expression(Dissolved_Oxygen_mg_per_L), x = expression(Time_Elapsed_min))+ theme(axis.text.x=element_text(size = 12,face="bold"))+
    #   ggtitle("No Points Removed " ,data_site_subset$Sample_ID[1]) +
    #   theme(plot.title = element_text(lineheight=.8, face="bold"))+
    #   theme(axis.text.x=element_text(colour = c("black","black")))+
    #   theme(aspect.ratio=1)+
    #   theme(axis.text.y=element_text(size = 12,face="bold"))+
    #   theme(axis.title.x =element_text(size = 12,face="bold"))+
    #   theme(axis.title =element_text(size = 12,face="bold"))+
    #   theme(axis.title.y =element_text(size = 12,face="bold"))
    # 
    # multi <- (final + high_rem + beg_rem + all) +
    #  plot_layout(widths = c(2,2))
    #ggsave(file=paste0(path,"Plots/combined_figures_pts_removed/20230217_Plots/inc_R2/DO_vs_Incubation_Time_",data_site_subset$Sample_ID[1],".pdf"))
    
    rate$Sample_ID[j] = as.character(data_site_subset_fin$Sample_ID[1])
    rate$slope_of_the_regression[j] = round(as.numeric((c)),3) #in mg O2/L min
    rate$rate_mg_per_L_per_min[j] = round(abs(as.numeric((c))),3) #in mg O2/L min
    rate$rate_mg_per_L_per_h[j] = round(abs(as.numeric((c))*60),3) #in mg O2/L h 
    rate$R_squared[j] = round(as.numeric(abs(r)),3)
    rate$R_squared_adj[j] = round(as.numeric(abs(r.adj)),3)
    rate$p_value[j] = p
    rate$residuals[j] = round(as.numeric(residuals),3)
    #rate$total_incubation_time_min[j] = as.numeric(difftime(data_site_subset$elapsed_min[nrow(data_site_subset_fin)],data_site_subset_fin$elapsed_min[1],units="mins"))#in minutes
    rate$number_of_points[j] = nrow(data_site_subset_fin)
    rate$removed_points_high[j] = (as.numeric(nrow(data_site_subset)) - as.numeric(nrow(data_site_subset_low)))
    
    rate$removed_points_beg[j] = (as.numeric(nrow(data_site_subset_low)) - as.numeric(nrow(data_site_subset_beg)))
    
    rate$removed_points_end[j] = (as.numeric(nrow(data_site_subset_beg)) - as.numeric(nrow(data_site_subset_fin)))
  }
  # Removing rows where the Sample_ID was an NA  
  rate = rate[!is.na(rate$Sample_ID),]   
  # Combining the regression metrics across locations and unique incubations
  respiration = rbind(respiration,rate)
}

respiration = respiration[-1,]

respiration = respiration %>% 
  mutate(slope_of_the_regression = if_else(slope_of_the_regression>0,0,slope_of_the_regression)) %>% 
  mutate(rate_mg = if_else(slope_of_the_regression>0,0,slope_of_the_regression)) 

write.csv(respiration,paste0(path,"Plots/ECA_Sediment_Incubations_Respiration_Rates_merged_by_",pnnl.user,"_on_",Sys.Date(),"_pts_rem.csv"), row.names = F)
