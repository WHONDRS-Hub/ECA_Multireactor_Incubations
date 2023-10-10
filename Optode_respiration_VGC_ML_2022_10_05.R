###### Load Library ######

library(dplyr); library(ggplot2);library(ggsignif)
library(ggpubr);library(reshape2);library(ggpmisc)
library(segmented);library(broom);library(lmtest);library(car)
library(ggpmisc);library(lubridate); library(readxl);
library(tidyverse);library(patchwork)
library(readr)


##### Load data ######
rm(list=ls());graphics.off()

# Set working directory to data file
#Example:
pnnl.user = 'laan208'

input.path = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/rates")

map.path = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/")

#setwd(input.path)

#path for raw data
path <- ("rates/")

#path for reading in 100% saturation values for each kit based on pressure/temperature during disk calibration
fast.sat <- paste0("C:/Users/", pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/rates/fast_rate_calculations.xlsx")

fast.rates <- read_excel(fast.sat)  

fast.rates.kits <- fast.rates %>% 
  rename("DO_mg_L" = "DO_sat_mg_L") 

#read in respiration data and clean
import_data = function(input.path){
  
  # pull a list of file names in the target folder with the target pattern
  # then read all files and combine
  
  filePaths <- list.files(path = input.path, recursive = T, pattern = "\\.csv$", full.names = TRUE)
  
  filePaths <- filePaths[!grepl("results_linear_fit|elapsed.time.ratios|archive|ECA_Sediment|removed|All_Respiration|results _03a|results_03b|results_01|results_02|results_04_08|results_10_15|results_06_07", filePaths)]
  
  # dat <- 
  do.call(rbind, lapply(filePaths, function(input.path){
    # then add a new column `source` to denote the file name
    df <- read.csv(input.path, skip = 4)
    df[["source_file"]] <- rep(input.path, nrow(df)) # add a column for the source file

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
      mutate(source_file = str_remove_all(source_file, paste0(input.path, "/")))
  }
  ))
}
data = import_data(input.path)

##### Clean Data ####

data_long = 
  data %>% 
  mutate(disc_number = str_remove_all(disc_number, "X")) %>%
  mutate(source_file = str_remove_all(source_file, ".*optode_data/"),source_file = str_remove_all(source_file, ".csv")) %>%  
  filter(elapsed_min > 0)

#import mapping files
import_data = function(map.path){ 
  
  #pull a list of files in target folder with correct pattern
  #read all files and combine
  
  map.file <-  list.files(map.path, recursive = T, pattern = "\\.xlsx$",full.names = T)
  
  map.file <- map.file[!grepl("red|Red|EC_01_|EC_02_|EC_03_|EC_06_07|EC_10_15|EC_04_08|fast|combined|SPC|IC|QA|Theoretical|theoretical", map.file)]
  
  mapping <- lapply(map.file, read_xlsx)
  
  for (i in 1:length(mapping)){mapping[[i]] <- cbind(mapping[[i]], map.file[i])}
  
  all.map <- 
    do.call(rbind,mapping)
}

all.map = import_data(map.path)

all.map = all.map %>% 
  rename("source_file" = "map.file[i]") %>% 
  rename("disc_number" = "Disk_ID") %>% 
  mutate(source_file = str_remove_all(source_file, ".*//")) %>% 
  separate(source_file, sep = "/", c("source_file", "file")) %>% 
  mutate(source_file = str_replace(source_file, "EC", "results")) %>% 
  dplyr::select(-`Disk Calibration date`, -Location, -`Time on`, -`Time off`, -SpC, -pH, -Temp, -Notes, -file)

all.samples <- merge(data_long, all.map)

# Fill in samples, remove last two columns for cleaner data frame

bind <- merge(all.samples, fast.rates.kits, all = TRUE) 

#not sure what happened with blank 1 from 14-18. starts to go to 0 at minute 10
blanks <- bind %>% 
  filter(grepl("Blank", Sample_ID)) %>% 
  filter(!grepl("results_14_18",source_file))

min <- min(blanks$DO_mg_L)
med <- median(blanks$DO_mg_L)

##### Read in times of pictures #####

import_data = function(map.path){ 
  
  #pull a list of files in target folder with correct pattern
  #read all files and combine
  
  time.map.file <-  list.files(map.path, recursive = T, pattern = "\\.xlsx$",full.names = T)
  
  time.map.file <- time.map.file[!grepl("red|Red|EC_01_|EC_02_|EC_03_|EC_06_07|EC_10_15|EC_04_08|fast|combined|SPC|IC|QAQC|EV|Theoretical|theoretical", time.map.file)]
  
  time.mapping <- lapply(time.map.file, read_xlsx)
  
  for (i in 1:length(time.mapping)){time.mapping[[i]] <- cbind(time.mapping[[i]], time.map.file[i])}
  
  time.map <- 
    do.call(rbind,time.mapping)
}

time.map = import_data(map.path)

time.map$`Time on` <- as.POSIXct(time.map$`Time on`, format = "%Y/%m/%d %H:%M:%%S")

time.map$`Time on` <- format(time.map$`Time on`, format = "%H:%M")

time.map$`Time off` <- as.POSIXct(time.map$`Time off`, format = "%Y/%m/%d %H:%M:%%S")

time.map$`Time off` <- format(time.map$`Time off`, format = "%H:%M")

time.map = time.map %>% 
  rename("source_file" = "time.map.file[i]") %>% 
  rename("disc_number" = "Disk_ID") %>% 
  mutate(source_file = str_remove_all(source_file, ".*//")) %>% 
  separate(source_file, sep = "/", c("source_file", "file")) %>% 
  mutate(source_file = str_replace(source_file, "EC", "results")) %>% 
  dplyr::select(-`Disk Calibration date`, -Location, -`Time off`, -SpC, -pH, -Temp, -Notes, -file)

import_data = function(map.path){ 
  
  filePaths <- list.files(path = map.path, recursive = T, pattern = "\\.txt$", full.names = TRUE)
  
  filePaths <- filePaths[!grepl("cal|images|EV|readme", filePaths)]

  mapping <- lapply(filePaths, read.delim)
  
  for (i in 1:length(mapping)){mapping[[i]] <- cbind(mapping[[i]], filePaths[i])}
  
  all.txt <- 
    do.call(rbind,mapping)
}

all.txt = import_data(map.path)  

img_all <- all.txt %>% 
  mutate(`filePaths[i]` = str_remove_all(`filePaths[i]`, paste0(map.path, "/"))) %>% 
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

samples <- merge(corr.time, bind, by = "Sample_ID", all = TRUE)

time_samples <- samples %>% 
  rename(source_file = source_file.x) %>% 
  rename(disc_number = disc_number.x) %>% 
  dplyr::select(-c(source_file.y, disc_number.y, Time_HMS, Time_HM, `Time on`)) %>% 
  separate(Sample_ID, into = c("Study Code", "Kit", "Rep"), sep = "_", remove = FALSE) %>% 
  filter(!grepl("Blank", `Study Code`))

#add "0" to start of sample kit names that don't have it

for (i in 1:nrow(time_samples)){
  
  if (str_count(time_samples$Kit[i], "[0-9]") <= 2){
    
    time_samples$Kit[i] = paste0("0", time_samples$Kit[i])
    
  }
  
  else {
    
    time_samples$Kit[i] = time_samples$Kit[i]
  }
  
}

time_samples <- time_samples %>% 
  unite("Sample_ID",  c("Study Code", "Kit", "Rep"), sep = "_") %>% 
  relocate(Sample_ID, .before = source_file) %>% 
  #dplyr::select(-`calibration date`)  %>% 
  filter(Sample_ID != "EC_013_INC-W5") %>% 
  filter(Sample_ID != "EC_013_INC-D4") %>% 
  filter(!(elapsed_min < 8 & Sample_ID == "EC_014_INC-W1")) %>%
  filter(Sample_ID != "EC_014_INC-W5") %>% 
  filter(Sample_ID != "EC_027_INC-D1") %>%
  filter(Sample_ID != "EC_027_INC-D2") #%>% 
# filter(!(elapsed_min == 2 & Sample_ID == "EC_32_INC-D3")) 

#13 - overexposed samples
#27 - overexposed
#14 kit - W1 second low point not being removed in script currently
#32-D3 kit - low first point, might be partially fixed with heteroscedasticity

# start loop to remove samples ####

#generate another dataset (w/ time, DO) with everything that has been removed: high, low (median), at the end, then can plot everything on top of each other with different colors to see what has been removed

min.points = 2
#threshold = 0.99
#res.threshold = 0.25
slope.thresh = -0.006
do.thresh = 2
#do.thresh = 1.07
time.thresh = 4
fast = 5.5 #5.5 has least number of differences between theoretical and real when choosing threshold
#ymax = max(na.omit(time_samples$DO_mg_L))
#ymin = min(na.omit(time_samples$DO_mg_L))
bpfit = 1.0
#daviesfit = 0.005
high.do = 14
p.value = 0.05
#mean.conc = 7.5
high.slope.thresh = -0.04 #keep saturation point if exponential
time.same = 7


respiration <- as.data.frame(matrix(NA, ncol = 16, nrow =1))

colnames(respiration) = c("Sample_ID","slope_of_the_regression", "rate_mg_per_L_per_min", "rate_mg_per_L_per_h", "R_squared", "R_squared_adj",  "p_value", "total_incubation_time_min", "number_of_points", "no_points_rem",  "breusch_p_value","break_point",  "0_min_concentration", "2_min_concentration", "last_concentration", "theoretical")

location = c("-W1", "-W2", "-W3","-W4", "-W5","-D1", "-D2", "-D3", "-D4", "-D5")

rate = as.data.frame(matrix(NA, ncol = 13, nrow = length(unique(time_samples$Sample_ID))))

colnames(rate) = c("Sample_ID", "slope_of_the_regression", "rate_mg_per_L_per_min", "rate_mg_per_L_per_h", "R_squared", "R_squared_adj", "p_value", "total_incubation_time_min", "breusch_p_value","0_min_concentration", "2_min_concentration", "last_concentration", "theoretical")


for (i in 1:length(location)){
  
 data_location_subset = time_samples[grep(location[i],time_samples$Sample_ID),]
  
 #data_location_subset = test
  
  unique.incubations = unique(data_location_subset$Sample_ID)
  
  rate = as.data.frame(matrix(NA, ncol = 13, nrow = length(unique(time_samples$Sample_ID))))
  
  colnames(rate) = c("Sample_ID", "slope_of_the_regression", "rate_mg_per_L_per_min", "rate_mg_per_L_per_h", "R_squared", "R_squared_adj", "p_value", "total_incubation_time_min", "breusch_p_value","0_min_concentration", "2_min_concentration","last_concentration", "theoretical")
  
  
  for (j in 1:length(unique.incubations)){
    
    data_site_subset = subset(data_location_subset, data_location_subset$Sample_ID == unique.incubations[j])
    
    data_site_subset = data_site_subset[order(data_site_subset$elapsed_min, decreasing = FALSE),]
    data_site_subset$elapsed_min = as.numeric(data_site_subset$elapsed_min)
    
    data_site_subset <- data_site_subset %>% 
      fill(break_toggle, .direction = c("down")) %>% 
      mutate(break_toggle = replace_na(break_toggle, "NA"))
    
    #fit linear model to all data
    
    fit_all = lm(DO_mg_L~elapsed_min, data = data_site_subset)
   
    #remove points greater than 14 mg/L
    
    data_site_subset_low = data_site_subset %>% 
      filter(DO_mg_L < high.do) 
    
    fit_low = lm(DO_mg_L~elapsed_min, data = data_site_subset_low)
    low.slope = fit_low$coefficients[[2]]
    
    
    
     ##remove samples if at >6 minutes, they are below the DO threshold. This is trying to remove low values from the back end of curves
    
    data_site_subset_thresh = data_site_subset_low %>% 
      filter(!(elapsed_min > 4 & DO_mg_L < do.thresh)) 
    
    data_site_subset_rem = data_site_subset_thresh
    
    for(n in 1:2){
      if(data_site_subset_rem$elapsed_min[2] == 2 & data_site_subset_rem$DO_mg_L[2] <= do.thresh){
        
        data_site_subset_rem = data_site_subset_rem %>% 
          filter(!elapsed_min > 2 )
        
      }
    }
    
    #remove saturation point from kits that have flat slopes
    
    data_site_subset_fast = data_site_subset_rem
    
    for(m in 1:1){
      
      if (data_site_subset_fast$DO_mg_L[2] <= fast & data_site_subset_fast$elapsed_min[2] <= time.thresh & low.slope < high.slope.thresh){
        
        data_site_subset_fast = data_site_subset_fast
        
      } 
      
      else if ((data_site_subset_fast$DO_mg_L[2] >= fast & data_site_subset_fast$elapsed_min[2] <= time.thresh)|(data_site_subset_fast$DO_mg_L[2] >= fast & data_site_subset_fast$elapsed_min[2] == 62)){
        
        data_site_subset_fast = data_site_subset_fast[-1,]
        
      }
    }
    
    #remove 2 minute point if put on at same time as picture is taken
    data_site_subset_beg = data_site_subset_fast
    
    for (k in 1:1) {

      if (nrow(data_site_subset_beg) <= 2) {

        data_site_subset_beg = data_site_subset_beg
      }

      else if ((data_site_subset_beg$time_same[k] == "yes" | data_site_subset_beg$time_same[k] == "maybe") & ((data_site_subset_beg$DO_mg_L[k] > time.same & data_site_subset_beg$elapsed_min[k] == 2))) {

        data_site_subset_beg = data_site_subset_beg[-k,]

      }

    }
    
      
      if (data_site_subset_beg$elapsed_min[1] == 0) {
        
        data_site_subset_beg = head(data_site_subset_beg, 2)
      
        }
      
      

    # 
    #fit cleaned data (low points removed, high points removed, first point removed if put on rollers at same time as picture) 
    
    fitog = lm(DO_mg_L~elapsed_min, data = data_site_subset_beg)
   
  #use if using only break point
  if (data_site_subset_beg$break_toggle[1] == "yes") { 
    if (nrow(data_site_subset_beg) > 3) {

    segmentog = segmented(fitog, seg.Z = ~elapsed_min, psi = (((max(data_site_subset_beg$elapsed_min)-min(data_site_subset_beg$elapsed_min))/2)+min(data_site_subset_beg$elapsed_min)))
    
    fit_seg = numeric(length(data_site_subset_beg$elapsed_min)) * NA

    fit_seg[complete.cases(rowSums(cbind(data_site_subset_beg$DO_mg_L, data_site_subset_beg$elapsed_min)))] <- broken.line(segmentog)$fit

    data_seg = data.frame(DO_mg_L = data_site_subset_beg$DO_mg_L, elapsed_min = data_site_subset_beg$elapsed_min, fit = fit_seg)

    }
   
    seg = segmentog$psi
    est = seg[[2]]
    data_site_subset_beg$break_point = 2*round(seg[[2]]/2)
     
  }
    
    else {
      
      data_site_subset_beg = data_site_subset_beg
      
      est = "NA"
      
    }

    
    data_site_subset_break = data_site_subset_beg
    
    ##remove points if break point = yes
    
    if (data_site_subset_break$break_toggle[1] == "yes") {
      
      data_site_subset_break = subset(data_site_subset_break, elapsed_min <= data_site_subset_break$break_point[1])
      
    }
    
    fit_break = lm(DO_mg_L~elapsed_min, data = data_site_subset_break)
    u = fit_break$coefficients
    b = u[[1]]
    c = u[[2]]
    slope_beginning = c
    r = summary(fit_break)$r.squared
    residuals = deviance(fit_break)
    r.adj = summary(fit_break)$adj.r.squared
    p = summary(fit_break)$coefficients[4]
    pog = p
    rog = r
    resog = residuals
    bp = bptest(fit_break)[[4]]
    bpog = bptest(fit_break)[[4]]
    
    data_site_subset_fin = data_site_subset_break
    
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
    
    if (slope_beginning >= slope.thresh | ((first(data_site_subset_fin$DO_mg_L) - last(data_site_subset_fin$DO_mg_L)) < 1.5)) #all(data_site_subset_fin$DO_mg_L >= mean.conc)) 
    {
      
      data_site_subset_fin = data_site_subset_fin
      
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
    
    else {
    
      for (l in 1:60) {
        
        if (nrow(data_site_subset_fin)<= min.points){
          
          data_site_subset_fin = data_site_subset_fin
          
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
        
        else if (bp < bpfit & nrow(data_site_subset_fin) >= min.points){
        
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

    clean <- ggplot(data_site_subset_beg, aes(x = elapsed_min, y = DO_mg_L)) + coord_cartesian(ylim = c(0,15))+ geom_point(size = 2) + expand_limits(x = 0, y = 0) +
      geom_smooth(method = "lm", se=F, formula = my.formula) +
      geom_label(aes(x = 0, y = 13), size = 2.5, hjust = 0,
                 label = paste("R2 = ", signif(summary(fitog)$r.squared, 3),
                               "\nP = ", signif(summary(fitog)$coefficients[[4]], 3),
                               "\nSlope = ", signif(fitog$coefficients[[2]], 3)))+
    #   # stat_poly_eq(formula = my.formula,label.y = "top",label.x = "right", aes(label = paste( ..rr.label.., sep = "~~~"),size=1), parse = TRUE)+stat_fit_glance(data=data_site_subset_low, method = 'lm', method.args = list(formula = my.formula),geom = 'text',aes(label =paste("         p = ",signif(..p.value.., digits = 1), sep = ""),size=1),label.y = c(14.25),label.x = "left") +

      theme_bw()+theme(legend.title = element_blank(),legend.background = element_rect(fill = 'NA'), legend.text = element_text(size = 12,face="bold"))+
      labs(y = expression(Dissolved_Oxygen_mg_per_L), x = expression(Time_Elapsed_min))+ theme(axis.text.x=element_text(size = 12,face="bold"))+
      ggtitle("Cleaned" ,data_site_subset_fast$Sample_ID[1]) +
      theme(plot.title = element_text(lineheight=.8, face="bold"))+
      theme(axis.text.x=element_text(colour = c("black","black")))+
      theme(aspect.ratio=1)+
      theme(axis.text.y=element_text(size = 12,face="bold"))+
      theme(axis.title.x =element_text(size = 12,face="bold"))+
      theme(axis.title =element_text(size = 12,face="bold"))+
      theme(axis.title.y =element_text(size = 12,face="bold"))

    break_rem <- ggplot(data_site_subset_break, aes(x = elapsed_min, y = DO_mg_L)) + coord_cartesian(ylim = c(0,15))+ geom_point(size = 2) + expand_limits(x = 0, y = 0) +
      geom_smooth(method = "lm", se=F, formula = my.formula) +
      geom_label(aes(x = 0, y = 13), size = 2.5, hjust = 0,
                 label = paste("R2 = ", signif(summary(fit_break)$r.squared, 3),
                               "\nP = ", signif(summary(fit_break)$coefficients[[4]], 3),
                               "\nSlope = ", signif(fit_break$coefficients[[2]], 3)))+
    #   # stat_poly_eq(formula = my.formula,label.y = "top",label.x = "right", aes(label = paste( ..rr.label.., sep = "~~~"),size=1), parse = TRUE)+stat_fit_glance(data=data_site_subset_beg, method = 'lm', method.args = list(formula = my.formula),geom = 'text',aes(label =paste("         p = ",signif(..p.value.., digits = 1), sep = ""),size=1),label.y = c(14.25),label.x = "left") +
      theme_bw()+theme(legend.title = element_blank(),legend.background = element_rect(fill = 'NA'), legend.text = element_text(size = 12,face="bold"))+
      labs(y = expression(Dissolved_Oxygen_mg_per_L), x = expression(Time_Elapsed_min))+ theme(axis.text.x=element_text(size = 12,face="bold"))+
      ggtitle("Break Point Removed ", data_site_subset_break$Sample_ID[1]) +
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

   multi <- (final + break_rem + clean + all) +
   plot_layout(widths = c(2,2))

   ggsave(file=paste0(input.path,"/Plots/Sensitivity_Analysis/No_Breusch/DO_vs_Incubation_Time_",data_site_subset$Sample_ID[1],".pdf"), width = 7, height = 7, units = "in")
   
   # pdf(file = paste0(input.path,"/Plots/Sensitivity Analysis/Sat_6.5_removal/All_Fits.pdf"))
   # for(m in 1:1){
   #   
   #   print(multi)
   #   
   # }
   # 
   # dev.off()
   
    rate$Sample_ID[j] = as.character(data_site_subset_fin$Sample_ID[1])
    rate$slope_of_the_regression[j] = round(as.numeric((c)),3) #in mg O2/L min
    rate$rate_mg_per_L_per_min[j] = round(abs(as.numeric((c))),3) #in mg O2/L min
    rate$rate_mg_per_L_per_h[j] = round(abs(as.numeric((c))*60),3) #in mg O2/L h 
    #rate$Initial_R_squared[j] = round(abs(as.numeric((rog))),3)
    rate$R_squared[j] = round(as.numeric(abs(r)),3)
    rate$R_squared_adj[j] = round(as.numeric(abs(r.adj)),3)
    rate$p_value[j] = p
    #rate$initial_p_value[j] = pog
    #rate$residuals[j] = round(as.numeric(residuals),3)
    #rate$slope_beginning[j] = round(as.numeric(slope_beginning),3)
    
    rate$total_incubation_time_min[j] = as.numeric(tail(data_site_subset_fin$elapsed_min, n = 1) - head(data_site_subset_fin$elapsed_min, n = 1)) 

    rate$number_of_points[j] = nrow(data_site_subset_fin)
    rate$no_points_rem[j] = (as.numeric(nrow(data_site_subset)) - as.numeric(nrow(data_site_subset_fin)))
    
    #rate$removed_points_beg[j] = (as.numeric(nrow(data_site_subset_low)) - as.numeric(nrow(data_site_subset_beg)))
    
    #rate$removed_points_end[j] = (as.numeric(nrow(data_site_subset_beg)) - as.numeric(nrow(data_site_subset_fin)))
    
    rate$breusch_p_value[j] = round(as.numeric(bp),3)
    #rate$davies_p_value[j] = round(as.numeric(davies),3)
    
    #rate$flag_r2[j] = case_when(rate$Initial_R_squared[j] > rate$Final_R_squared[j] ~ "Final R2 < Initial R2")
    
    #rate$flag_pos_slope[j] = case_when(rate$slope_beginning[j] < 0 & rate$slope_of_the_regression[j] > 0~ "Final Slope positive - initial slope negative, reporting initial slope")
    
   #rate$flag_heteroscedastic[j] = case_when( rate$breusch_p_value[j] <= 0.1 ~ "Heteroscedastic but not removing data")
    
    rate$`0_min_concentration`[j] = ifelse(
      first(data_site_subset_fin$elapsed_min == 0), 
        subset(data_site_subset_fin$DO_mg_L,                             data_site_subset_fin$elapsed_min == 0),
      "NA"
      )

rate$`2_min_concentration`[j] = ifelse(
  first(data_site_subset_fin$elapsed_min == 4), 
    subset(data_site_subset_fin$DO_mg_L, 
    data_site_subset_fin$elapsed_min == 4),
  subset(data_site_subset_fin$DO_mg_L,
  data_site_subset_fin$elapsed_min == 2)) 
    
    rate$last_concentration[j] = last(data_site_subset_break$DO_mg_L)
    
    rate$break_point[j] = est
    
    rate$theoretical[j] = case_when(data_site_subset_fin$elapsed_min[1] == 0 ~ "yes")
    
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

write.csv(respiration,paste0(input.path,"/Plots/ECA_Sediment_Incubations_Respiration_Rates_No_Breusch_on_",Sys.Date(),".csv"), row.names = F)


## Pull all exported pdfs into one file
library(pdftools)

setwd(paste0(input.path,"/Plots/Sensitivity_Analysis/No_Breusch/"))

ex_pdf <- list.files(pattern = "pdf")

ex_pdf1 <- head(ex_pdf, 275)
ex_pdf2 <- tail(ex_pdf, 278)

pdf_combine(ex_pdf1, output = "First_Half.pdf")
pdf_combine(ex_pdf2, output = "Second_Half.pdf")

all_pdf <- list.files(pattern = "Half.pdf")
pdf_combine(all_pdf, output = "All_Fits.pdf")
