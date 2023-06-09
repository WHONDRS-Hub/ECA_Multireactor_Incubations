###### Load Library ######

library(dplyr); library(ggplot2);library(ggsignif)
library(ggpubr);library(reshape2);library(ggpmisc)
library(segmented);library(broom);library(lmtest)
library(ggpmisc);library(segmented);library(lubridate); library(readxl);
library(tidyverse);library(patchwork)
library(readr)

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

fast.rates <- read_excel(fast.sat)  
  
#kits that go to 0 quickly - 23W, 23D, 27W, 40W, 32D2, 32D4, 12W, 11W,  34W, 35W, 35D1245, 69W

fast.rates.kits <- fast.rates %>% 
  rename("DO_mg_L" = "DO_sat_mg_L") 

#read in respiration data and clean
import_data = function(path){
  
  # pull a list of file names in the target folder with the target pattern
  # then read all files and combine
  
  filePaths <- list.files(path = path, recursive = T, pattern = "\\.csv$", full.names = TRUE)
  
  filePaths <- filePaths[!grepl("results_linear_fit|elapsed.time.ratios|archive|ECA_Sediment|removed", filePaths)]
  
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


##### Clean Data ####

data_long = 
  data %>% 
  mutate(disc_number = str_remove_all(disc_number, "X")) %>%
  mutate(source_file = str_remove_all(source_file, "optode data/"),source_file = str_remove_all(source_file, ".csv")) %>%  
  filter(elapsed_min > 0) %>% 
  separate(col = source_file, into = c("fol1", "fol2", "fol3","fol4", "source_file"), sep = "/") %>% 
  dplyr::select(-fol1, -fol2, -fol3,-fol4)


##take out samples incubated in jars
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
  
  map.file <- map.file[!grepl("red|Red|EC_01_|EC_02_|EC_03_|EC_06_07|EC_10_15|EC_04_08|fast|combined|SPC|IC|QA", map.file)]
  
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

#fill in samples, remove last two columns for cleaner data frame

bind <- merge(all.samples, fast.rates.kits, all = TRUE) %>% 
  dplyr::select(-`calibration date`) 

#13 - overexposed samples
#27 - overexposed
#14 kit - W1 second low point not being removed in script currently
#32-D3 kit - low first point, might be partially fixed with heteroscedasticity

bind <- bind %>% 
  filter(Sample_ID != "EC_13_INC-W5") %>% 
  filter(Sample_ID != "EC_13_INC-D4") %>% 
  filter(!(elapsed_min < 8 & Sample_ID == "EC_14_INC-W1")) %>%
  filter(Sample_ID != "EC_14_INC-W5") %>% 
  filter(Sample_ID != "EC_27_INC-D1") %>%
  filter(Sample_ID != "EC_27_INC-D2") #%>% 
 # filter(!(elapsed_min == 2 & Sample_ID == "EC_32_INC-D3")) 


##### Read in times of pictures #####

import_data = function(input.path){ 
  
  #pull a list of files in target folder with correct pattern
  #read all files and combine
  
  time.map.file <-  list.files(input.path, recursive = T, pattern = "\\.xlsx$",full.names = T)
  
  time.map.file <- time.map.file[!grepl("red|Red|EC_01_|EC_02_|EC_03_|EC_06_07|EC_10_15|EC_04_08|fast|combined|SPC|IC|QAQC", time.map.file)]
  
  time.mapping <- lapply(time.map.file, read_xlsx)
  
  for (i in 1:length(time.mapping)){time.mapping[[i]] <- cbind(time.mapping[[i]], time.map.file[i])}
  
  time.map <- 
    do.call(rbind,time.mapping)
}

time.map = import_data(input.path)

time.map$`Time on` <- as.POSIXct(time.map$`Time on`, format = "%Y/%m/%d %H:%M:%%S")

time.map$`Time on` <- format(time.map$`Time on`, format = "%H:%M")

time.map$`Time off` <- as.POSIXct(time.map$`Time off`, format = "%Y/%m/%d %H:%M:%%S")

time.map$`Time off` <- format(time.map$`Time off`, format = "%H:%M")

time.map = time.map %>% 
  rename("source_file" = "time.map.file[i]") %>% 
  rename("disc_number" = "Disk_ID") %>% 
  mutate(source_file = str_remove_all(source_file, paste0(input.path, "/"))) %>% 
  separate(source_file, sep = "/", c("source_file", "file")) %>% 
  mutate(source_file = str_replace(source_file, "EC", "results")) %>% 
  dplyr::select(-`Disk Calibration date`, -Location, -`Time off`, -SpC, -pH, -Temp, -Notes, -file)

import_data = function(input.path){ 
  
  filePaths <- list.files(path = input.path, recursive = T, pattern = "\\.txt$", full.names = TRUE)
  
  filePaths <- filePaths[!grepl("cal", filePaths)]
  filePaths <- filePaths[!grepl("images", filePaths)]
  
  mapping <- lapply(filePaths, read.delim)
  
  for (i in 1:length(mapping)){mapping[[i]] <- cbind(mapping[[i]], filePaths[i])}
  
  all.txt <- 
    do.call(rbind,mapping)
}

all.txt = import_data(input.path)  

img_all <- all.txt %>% 
  mutate(`filePaths[i]` = str_remove_all(`filePaths[i]`, paste0(input.path, "/"))) %>% 
  separate(col = `filePaths[i]`, into = c("source_file", "photo"), sep = "/") %>% 
  mutate(source_file = str_replace(source_file,"EC", "results")) 

img_time <- img_all[!grepl("Custom|type|.tif|----", img_all$TIFF.image.set.saved.with.Look.RGB.v0.1),]

img_time <- rename(img_time, Time = TIFF.image.set.saved.with.Look.RGB.v0.1)

img_time$Time <- gsub('[AMP]','', img_time$Time)

img_time$Time <- as.POSIXct(img_time$Time, format = "%m/%d/%Y %H:%M:%S")

img_time$Day <- format(img_time$Time, format = "%m/%d/%Y")

img_time$Day <- as.Date(img_time$Day, format = "%m/%d/%Y")

img_time$Time_Corr <- ifelse(img_time$Day < "2022-11-06" | img_time$Day > "2023-03-12", img_time$Time + 3600, img_time$Time)

class(img_time$Time_Corr) <- c("POSIXct", "POSIXt")

img_time$Time_HMS <- format(as.POSIXct(img_time$Time_Corr), format = "%H:%M:%S")

img_time$Time_HM <- format(as.POSIXct(img_time$Time_Corr), format = "%H:%M")

img_time$Time_S <- format(as.POSIXct(img_time$Time_Corr), format = "%S")

all.samples <- merge(img_time, time.map)

all.samples$min_bef <- format(strptime(all.samples$Time_HM, format = "%H:%M") - 60, "%H:%M")

all.samples$time_same <- NA

for (i in 1:nrow(all.samples)) {
  
  if (all.samples$Time_HM[i] == all.samples$`Time on`[i]) {
    
    all.samples$time_same[i] <- "yes"
    
  }
  
  else if (all.samples$Time_S[i] <= 20 & all.samples$`Time on`[i] == all.samples$min_bef[i]) {
  
    all.samples$time_same[i] <- "maybe"
    
  }
  
  else { 
    
    all.samples$time_same[i] <- "no"
    
    }
  
}

corr.time <- all.samples %>% 
  dplyr::select(c(source_file, Sample_ID, disc_number, Time_HMS, Time_HM, `Time on`, time_same)) %>% 
  filter(!grepl("Blank", Sample_ID)) %>% 
  group_by(Sample_ID) %>% 
  filter(Time_HM >= `Time on`) %>% 
  arrange(Time_HM) %>% 
  filter(row_number() == 1)

time_samples <- merge(corr.time, bind, by = "Sample_ID")

time_samples <- time_samples %>% 
  rename(source_file = source_file.x) %>% 
  rename(disc_number = disc_number.x) %>% 
  dplyr::select(-c(source_file.y, disc_number.y, Time_HMS, Time_HM, `Time on`))

# start loop to remove samples ####

#generate another dataset (w/ time, DO) with everything that has been removed: high, low (median), at the end, then can plot everything on top of each other with different colors to see what has been removed


min.points = 2
threshold = 0.99
res.threshold = 0.25
slope.thresh = -0.006
do.thresh = 1
time.thresh = 2
fast = 5
ymax = max(na.omit(bind$DO_mg_L))
ymin = min(na.omit(bind$DO_mg_L))
bpfit = 0.1
high.do = 14
p.value = 0.05
mean.conc = 7.5


respiration <- as.data.frame(matrix(NA, ncol = 20, nrow =1))

colnames(respiration) = c("Sample_ID","slope_beginning", "slope_of_the_regression", "rate_mg_per_L_per_min", "rate_mg_per_L_per_h","Initial_R_squared", "Final_R_squared", "R_squared_adj","residuals","initial_p_value", "final_p_value", "total_incubation_time_min", "number_of_points", "removed_points_high", "removed_points_beg", "removed_points_end", "breusch_p_value","flag_r2", "flag_pos_slope", "flag_heteroscedastic")

location = c("-W1", "-W2", "-W3", "-W4", "-W5", "-D1", "-D2", "-D3", "-D4", "-D5")

rate = as.data.frame(matrix(NA, ncol = 16, nrow = length(unique(bind$Sample_ID))))

colnames(rate) = c("Sample_ID","slope_beginning", "slope_of_the_regression", "rate_mg_per_L_per_min", "rate_mg_per_L_per_h","Initial_R_squared", "Final_R_squared", "R_squared_adj","residuals", "initial_p_value", "final_p_value", "total_incubation_time_min", "breusch_p_value","flag_r2", "flag_pos_slope", "flag_heteroscedastic")


for (i in 1:length(location)){
  
  data_location_subset = time_samples[grep(location[i],time_samples$Sample_ID),]
  
  unique.incubations = unique(data_location_subset$Sample_ID)
  
  rate = as.data.frame(matrix(NA, ncol = 16, nrow = length(unique(bind$Sample_ID))))
  
  colnames(rate) = c("Sample_ID","slope_beginning", "slope_of_the_regression", "rate_mg_per_L_per_min", "rate_mg_per_L_per_h","Initial_R_squared", "Final_R_squared", "R_squared_adj","residuals", "initial_p_value", "final_p_value", "total_incubation_time_min", "breusch_p_value","flag_r2", "flag_pos_slope", "flag_heteroscedastic")
  
  
  for (j in 1:length(unique.incubations)){
    
    data_site_subset = subset(data_location_subset, data_location_subset$Sample_ID == unique.incubations[j])
    data_site_subset = data_site_subset[order(data_site_subset$elapsed_min, decreasing = FALSE),]
    data_site_subset$elapsed_min = as.numeric(data_site_subset$elapsed_min)
    
    fit_all = lm(data_site_subset$DO_mg_L~data_site_subset$elapsed_min)
   
    #remove saturation point from kits that didn't go to 0 quickly  
    
    data_site_subset_fast = data_site_subset
    
    for(m in 1:2){
      
      if (data_site_subset_fast$DO_mg_L[2] <= fast & data_site_subset_fast$elapsed_min[2] == time.thresh){
        
        data_site_subset_fast = data_site_subset_fast
        
      } 
      
      else if (data_site_subset_fast$DO_mg_L[2] >= fast & data_site_subset_fast$elapsed_min[2] == time.thresh){
        
        data_site_subset_fast = data_site_subset_fast[-1,]
        
      }
    }
    
    #remove points greater than 14 mg/L
    
    data_site_subset_low = data_site_subset_fast %>% 
      filter(DO_mg_L < high.do) 
    
    fit_low = lm(data_site_subset_low$DO_mg_L~data_site_subset_low$elapsed_min)
    
     ##remove samples if at >4 minutes, they are below the DO threshold. This is trying to remove low values from the back end of curves
    
    data_site_subset_thresh = data_site_subset_low %>% 
      filter(!(elapsed_min > time.thresh & DO_mg_L < do.thresh)) 
    
    data_site_subset_rem = data_site_subset_thresh
    
    for(n in 1:2){
      if(data_site_subset_rem$elapsed_min[2] == 2 & data_site_subset_rem$DO_mg_L[2] <= 5){
        
        data_site_subset_rem = data_site_subset_rem %>% 
          filter(!elapsed_min > 2 )
        
      }
    }
    
    data_site_subset_beg = data_site_subset_rem
    
    for (k in 1:2) {
      
      if (nrow(data_site_subset_beg) <= 2) {
        
        data_site_subset_beg = data_site_subset_beg
      }
      
      else if ((data_site_subset_beg$time_same[k] == "yes" | data_site_subset_beg$time_same[k] == "maybe") & ((data_site_subset_beg$DO_mg_L[k] > fast & data_site_subset_beg$elapsed_min[k] == 2))) {
        
        data_site_subset_beg = data_site_subset_beg[-k,]
        
      }
      
  }
 
    fitog = lm(data_site_subset_beg$DO_mg_L~data_site_subset_beg$elapsed_min)
    
    data_site_subset_fin = data_site_subset_beg
    
     fit = lm(data_site_subset_fin$DO_mg_L~data_site_subset_fin$elapsed_min)
    u = fit$coefficients
    b = u[[1]]
    c = u[[2]]
    slope_beginning = c
    r = summary(fit)$r.squared
    residuals = deviance(fit)
    r.adj = summary(fit)$adj.r.squared
    p = summary(fit)$coefficients[4]
    pog = p
    rog = r
    resog = residuals
    bp = bptest(fit)[[4]]
    bpog = bptest(fit)[[4]]
    
    temp = data_site_subset_fin[-nrow(data_site_subset_fin),]
    fit.temp = lm(temp$DO_mg_L~temp$elapsed_min)
    rtemp = summary(fit.temp)$r.squared
    restemp = deviance(fit.temp)
    bptemp = bptest(fit.temp)[[4]]
    

    if (slope_beginning >= slope.thresh | all(data_site_subset_fin$DO_mg_L >= mean.conc)) {
      
      data_site_subset_fin = data_site_subset_fin
      
    }
    
  
   
    else {
     
      for(l in 1:60){
      
      if (nrow(data_site_subset_fin)<= min.points){
        
        data_site_subset_fin = data_site_subset_fin
        
      }
        
      
     
      else if (bp < bpfit & nrow(data_site_subset_fin)>=min.points){

        data_site_subset_fin = data_site_subset_fin[-nrow(data_site_subset_fin),]

        fit = lm(data_site_subset_fin$DO_mg_L~data_site_subset_fin$elapsed_min)
        u = fit$coefficients
        b = u[[1]] #Intercept
        c = u[[2]] #rate mg/L min
        r = summary(fit)$r.squared
        r.adj = summary(fit)$adj.r.squared
        residuals = deviance(fit)
        p = summary(fit)$coefficients[4]
        r2 = r
        res2 = residuals
        bp = bptest(fit)[[4]]
        
      }
      
      }
      
      if (slope_beginning < 0 & c > 0){
        
        data_site_subset_fin = data_site_subset_thresh
        
      }
      
    }
   
  my.format <- "Slope: %s\nR2: %s\np: %s"  
  my.formula <- y ~ x
     final <- ggplot(data_site_subset_fin, aes(x = elapsed_min, y = DO_mg_L)) + coord_cartesian(ylim = c(0,15))+ geom_point(size = 2) + expand_limits(x = 0, y = 0) +
      geom_smooth(method = "lm", se=F, formula = my.formula) +
      geom_label(aes(x = 0, y = 13), size = 2.5, hjust = 0,
      label = paste("R2 = ", signif(summary(fit)$r.squared, 3),
                 "\nP = ", signif(summary(fit)$coefficients[[4]], 3),
                 "\nSlope = ", signif(fit$coefficients[[2]], 3)))+
      #geom_text(size = 10)+
      #stat_poly_eq(formula = my.formula,label.y = "top",label.x = "right", aes(label = paste( ..rr.label.., sep = "~~~"),size=1), parse = TRUE)+stat_fit_glance(data=data_site_subset_fin, method = 'lm', method.args = list(formula = my.formula),geom = 'text',aes(label =paste("          p = ",signif(..p.value.., digits = 1), sep = ""),size=1),label.y = c(14.25),label.x = "left") +
      theme_bw()+theme(legend.title = element_blank(),legend.background = element_rect(fill = 'NA'), legend.text = element_text(size = 12,face="bold"))+
      labs(y = expression(Dissolved_Oxygen_mg_per_L), x = expression(Time_Elapsed_min))+ theme(axis.text.x=element_text(size = 12,face="bold"))+
      ggtitle("Final " ,data_site_subset_fin$Sample_ID[1]) +
      theme(plot.title = element_text(lineheight=.8, face="bold"))+
      theme(axis.text.x=element_text(colour = c("black","black")))+
      theme(aspect.ratio=1)+
      theme(axis.text.y=element_text(size = 12,face="bold"))+
      theme(axis.title.x =element_text(size = 12,face="bold"))+
      theme(axis.title =element_text(size = 12,face="bold"))+
      theme(axis.title.y =element_text(size = 12,face="bold"))

    high_rem <- ggplot(data_site_subset_low, aes(x = elapsed_min, y = DO_mg_L)) + coord_cartesian(ylim = c(0,15))+ geom_point(size = 2) + expand_limits(x = 0, y = 0) +
      geom_smooth(method = "lm", se=F, formula = my.formula) +
      geom_label(aes(x = 0, y = 13), size = 2.5, hjust = 0,
                 label = paste("R2 = ", signif(summary(fit_low)$r.squared, 3),
                               "\nP = ", signif(summary(fit_low)$coefficients[[4]], 3),
                               "\nSlope = ", signif(fit_low$coefficients[[2]], 3)))+
      # stat_poly_eq(formula = my.formula,label.y = "top",label.x = "right", aes(label = paste( ..rr.label.., sep = "~~~"),size=1), parse = TRUE)+stat_fit_glance(data=data_site_subset_low, method = 'lm', method.args = list(formula = my.formula),geom = 'text',aes(label =paste("         p = ",signif(..p.value.., digits = 1), sep = ""),size=1),label.y = c(14.25),label.x = "left") +
      theme_bw()+theme(legend.title = element_blank(),legend.background = element_rect(fill = 'NA'), legend.text = element_text(size = 12,face="bold"))+
      labs(y = expression(Dissolved_Oxygen_mg_per_L), x = expression(Time_Elapsed_min))+ theme(axis.text.x=element_text(size = 12,face="bold"))+
      ggtitle("High Points Removed " ,data_site_subset_low$Sample_ID[1]) +
      theme(plot.title = element_text(lineheight=.8, face="bold"))+
      theme(axis.text.x=element_text(colour = c("black","black")))+
      theme(aspect.ratio=1)+
      theme(axis.text.y=element_text(size = 12,face="bold"))+
      theme(axis.title.x =element_text(size = 12,face="bold"))+
      theme(axis.title =element_text(size = 12,face="bold"))+
      theme(axis.title.y =element_text(size = 12,face="bold"))

    beg_rem <- ggplot(data_site_subset_beg, aes(x = elapsed_min, y = DO_mg_L)) + coord_cartesian(ylim = c(0,15))+ geom_point(size = 2) + expand_limits(x = 0, y = 0) +
      geom_smooth(method = "lm", se=F, formula = my.formula) +
      geom_label(aes(x = 0, y = 13), size = 2.5, hjust = 0,
                 label = paste("R2 = ", signif(summary(fitog)$r.squared, 3),
                               "\nP = ", signif(summary(fitog)$coefficients[[4]], 3),
                               "\nSlope = ", signif(fitog$coefficients[[2]], 3)))+
      # stat_poly_eq(formula = my.formula,label.y = "top",label.x = "right", aes(label = paste( ..rr.label.., sep = "~~~"),size=1), parse = TRUE)+stat_fit_glance(data=data_site_subset_beg, method = 'lm', method.args = list(formula = my.formula),geom = 'text',aes(label =paste("         p = ",signif(..p.value.., digits = 1), sep = ""),size=1),label.y = c(14.25),label.x = "left") +
      theme_bw()+theme(legend.title = element_blank(),legend.background = element_rect(fill = 'NA'), legend.text = element_text(size = 12,face="bold"))+
      labs(y = expression(Dissolved_Oxygen_mg_per_L), x = expression(Time_Elapsed_min))+ theme(axis.text.x=element_text(size = 12,face="bold"))+
      ggtitle("First Point Removed ", data_site_subset_beg$Sample_ID[1]) +
      theme(plot.title = element_text(lineheight=.8, face="bold"))+
      theme(axis.text.x=element_text(colour = c("black","black")))+
      theme(aspect.ratio=1)+
      theme(axis.text.y=element_text(size = 12,face="bold"))+
      theme(axis.title.x =element_text(size = 12,face="bold"))+
      theme(axis.title =element_text(size = 12,face="bold"))+
      theme(axis.title.y =element_text(size = 12,face="bold"))

    all <- ggplot(data_site_subset, aes(x = elapsed_min, y = DO_mg_L)) + coord_cartesian(ylim = c(0,15))+ geom_point(size = 2) + expand_limits(x = 0, y = 0) +
      geom_smooth(method = "lm", se=F, formula = my.formula) +
      geom_label(aes(x = 0, y = 13), size = 2.5, hjust = 0,
                 label = paste("R2 = ", signif(summary(fit_all)$r.squared, 3),
                               "\nP = ", signif(summary(fit_all)$coefficients[[4]], 3),
                               "\nSlope = ", signif(fit_all$coefficients[[2]], 3)))+
      # stat_poly_eq(formula = my.formula,label.y = "top",label.x = "right", aes(label = paste( ..rr.label.., sep = "~~~"),size=1), parse = TRUE)+stat_fit_glance(data=data_site_subset, method = 'lm', method.args = list(formula = my.formula),geom = 'text',aes(label =paste("         p = ",signif(..p.value.., digits = 1), sep = ""),size=1),label.y = c(14.25),label.x = "left") +
      theme_bw()+theme(legend.title = element_blank(),legend.background = element_rect(fill = 'NA'), legend.text = element_text(size = 12,face="bold"))+
      labs(y = expression(Dissolved_Oxygen_mg_per_L), x = expression(Time_Elapsed_min))+ theme(axis.text.x=element_text(size = 12,face="bold"))+
      ggtitle("No Points Removed " ,data_site_subset$Sample_ID[1]) +
      theme(plot.title = element_text(lineheight=.8, face="bold"))+
      theme(axis.text.x=element_text(colour = c("black","black")))+
      theme(aspect.ratio=1)+
      theme(axis.text.y=element_text(size = 12,face="bold"))+
      theme(axis.title.x =element_text(size = 12,face="bold"))+
      theme(axis.title =element_text(size = 12,face="bold"))+
      theme(axis.title.y =element_text(size = 12,face="bold"))

    multi <- (final + beg_rem + high_rem + all) +
     plot_layout(widths = c(2,2))
    ggsave(file=paste0(path,"Plots/breusch_test_fits/4-26-2023/Same Time Removed/DO_vs_Incubation_Time_",data_site_subset$Sample_ID[1],".pdf"))

    rate$Sample_ID[j] = as.character(data_site_subset_fin$Sample_ID[1])
    rate$slope_of_the_regression[j] = round(as.numeric((c)),3) #in mg O2/L min
    rate$rate_mg_per_L_per_min[j] = round(abs(as.numeric((c))),3) #in mg O2/L min
    rate$rate_mg_per_L_per_h[j] = round(abs(as.numeric((c))*60),3) #in mg O2/L h 
    rate$Initial_R_squared[j] = round(abs(as.numeric((rog))),3)
    rate$Final_R_squared[j] = round(as.numeric(abs(r)),3)
    rate$R_squared_adj[j] = round(as.numeric(abs(r.adj)),3)
    rate$final_p_value[j] = p
    rate$initial_p_value[j] = pog
    rate$residuals[j] = round(as.numeric(residuals),3)
    rate$slope_beginning[j] = round(as.numeric(slope_beginning),3)
    #rate$total_incubation_time_min[j] = as.numeric(difftime(data_site_subset$elapsed_min[nrow(data_site_subset_fin)],data_site_subset_fin$elapsed_min[1],units="mins"))#in minutes
    rate$number_of_points[j] = nrow(data_site_subset_fin)
    rate$removed_points_high[j] = (as.numeric(nrow(data_site_subset)) - as.numeric(nrow(data_site_subset_low)))
    
    rate$removed_points_beg[j] = (as.numeric(nrow(data_site_subset_low)) - as.numeric(nrow(data_site_subset_beg)))
    
    rate$removed_points_end[j] = (as.numeric(nrow(data_site_subset_beg)) - as.numeric(nrow(data_site_subset_fin)))
    
    rate$breusch_p_value[j] = round(as.numeric(bp),3)
    
    rate$flag_r2[j] = case_when(rate$Initial_R_squared[j] > rate$Final_R_squared[j] ~ "Final R2 < Initial R2")
    
    rate$flag_pos_slope[j] = case_when(rate$slope_beginning[j] < 0 & rate$slope_of_the_regression[j] > 0~ "Final Slope positive - initial slope negative, reporting initial slope")
    
    rate$flag_heteroscedastic[j] = case_when( rate$breusch_p_value[j] <= 0.1 ~ "Heteroscedastic but not removing data")
    
    }
  
  # Removing rows where the Sample_ID was an NA  
  #rate$flag[j] = case_when(
    #rate$initial_p_value[j] < p.value & rate$final_p_value[j] > p.value ~ "Final P-Value significantly different from 0",
    
  
  rate = rate[!is.na(rate$Sample_ID),]   
  # Combining the regression metrics across locations and unique incubations
  respiration = rbind(respiration,rate)
}

respiration = respiration[-1,]


#moved this to effect size code

# respiration = respiration %>% 
#   mutate(slope_of_the_regression = if_else(slope_of_the_regression>0,0,slope_of_the_regression)) %>% 
#   mutate(rate_mg = if_else(slope_of_the_regression>0,0,slope_of_the_regression)) 

write.csv(respiration,paste0(path,"Plots/All_Respiration_Rates/ECA_Sediment_Incubations_Respiration_Rates_merged_by_",pnnl.user,"_on_",Sys.Date(),".csv"), row.names = F)
