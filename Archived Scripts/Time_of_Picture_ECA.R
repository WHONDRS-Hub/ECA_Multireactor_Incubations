library(dplyr); library(tidyverse); library(lubridate);library(readxl);library(readr);library(hms)


rm(list=ls());graphics.off()

pnnl.user = 'laan208'

input.path = paste0("Y:/Optode_multi_reactor/Optode_multi_reactor_incubation/")


setwd(input.path)


import_data = function(input.path){ 
  
  #pull a list of files in target folder with correct pattern
  #read all files and combine
  
  map.file <-  list.files(input.path, recursive = T, pattern = "\\.xlsx$",full.names = T)
  
  #map.file <- map.file[!grepl("red|Red|EC_01_|EC_02_|EC_03_|EC_06_07|EC_10_15|EC_04_08|fast|combined|SPC|IC|QAQC", map.file)]
  
  map.file <- map.file[grepl("EV", map.file)]
  
  mapping <- lapply(map.file, read_xlsx)
  
  for (i in 1:length(mapping)){mapping[[i]] <- cbind(mapping[[i]], map.file[i])}
  
  all.map <- 
    do.call(rbind,mapping)
}

all.map = import_data(input.path)

all.map$`Time on` <- as.POSIXct(all.map$`Time on`, format = "%Y/%m/%d %H:%M:%%S")

all.map$`Time on` <- format(all.map$`Time on`, format = "%H:%M:%S")

all.map$`Time off` <- as.POSIXct(all.map$`Time off`, format = "%Y/%m/%d %H:%M:%%S")

all.map$`Time off` <- format(all.map$`Time off`, format = "%H:%M:%S")


all.map = all.map %>% 
  rename("source_file" = "map.file[i]") %>% 
  rename("disc_number" = "Disk_ID") %>% 
  mutate(source_file = str_remove_all(source_file, paste0(input.path, "/"))) %>% 
  separate(source_file, sep = "/", c("source_file", "file")) %>% 
  mutate(source_file = str_replace(source_file, "EC", "results")) %>% 
  dplyr::select(-`Disk Calibration date`, -Location, -`Time off`, -SpC, -pH, -Temp, -Notes, -file)

  # pull a list of file names in the target folder with the target pattern
  # then read all files and combine
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

#img_time <- slice(img_all, seq(1,nrow(img_all),8))

 img_time <- img_all[!grepl("Custom|type|.tif|----", img_all$TIFF.image.set.saved.with.Look.RGB.v0.1),]

img_time <- rename(img_time, Time = TIFF.image.set.saved.with.Look.RGB.v0.1)

img_time$Time <- gsub('[AMP]','', img_time$Time)

img_time$Time <- as.POSIXct(img_time$Time, format = "%m/%d/%Y %H:%M:%S")

img_time$Date <- format(as.POSIXct(img_time$Time), format = "%m/%d/%Y")
 
img_time$Time_S <- format(img_time$Time, format = "%H:%M:%S")

img_time$Time_S <- hms::as_hms(img_time$Time_S)

img_time$Time_H <- format(as.POSIXct(img_time$Time), format = "%H:%M")

img_time$Time_H <-format(img_time$Time_H, format = "%H:%M:%S") 

all.samples <- merge(img_time, all.map)

all.samples$`Time on` <- hms::as_hms(all.samples$`Time on`)

all.samples


time.remove <- all.samples %>% 
  mutate(same = `Time on` >= (Time_S - 20) & 
           `Time on` <= (Time_S) | `Time on` == Time_H)

all.samples$time.add = all.samples$Time_S + 20

time.same <- all.samples %>% 
  filter(all.samples$Time_ == all.samples$`Time on`)

time.same <- all.samples[all.samples$Time==all.samples$`Time on`,]   
