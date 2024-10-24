###### Load Library ######

library(tidyverse); library(readxl)

##### Load data ######
rm(list=ls());graphics.off()

# Set paths to data files and outputs

## INPUTS ####
pnnl.user = 'laan208'

study.code = 'EL_'
  #EC_
  #EV_
  #RS_
  #EL_

if (study.code == "EC_") {
  
  fast.rates.in = 'ec_fast_rate_calculations.xlsx'
  
  ## For .txt files for image times
  input.path = ("Y:/Optode_multi_reactor/Optode_multi_reactor_incubation/")
  
  ## For mapping files and result .csv's
  map.path = paste0("C:/Users/",pnnl.user,"/OneDrive - PNNL/Shared Documents - Core Richland and Sequim Lab-Field Team/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/")
  
  ## Where to put Final Raw DO file
  output.path = paste0("C:/Users/",pnnl.user,"/OneDrive - PNNL/Shared Documents - Core Richland and Sequim Lab-Field Team/Data Generation and Files/ECA/INC/03_ProcessedData/Raw_DO_by_Min/")
  
  
} else if (study.code == "EV_") {
  
  fast.rates.in = 'ev_fast_rate_calculations.xlsx'
  
  ## For .txt files for image times
  input.path = ("Y:/Optode_multi_reactor/Optode_multi_reactor_incubation/")
  
  ## For mapping files and result .csv's
  map.path = paste0("C:/Users/",pnnl.user,"/OneDrive - PNNL/Shared Documents - Core Richland and Sequim Lab-Field Team/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/")
  
  ## Where to put Final Raw DO file
  output.path = paste0("C:/Users/",pnnl.user,"/OneDrive - PNNL/Shared Documents - Core Richland and Sequim Lab-Field Team/Data Generation and Files/ECA/INC/03_ProcessedData/Raw_DO_by_Min/")
  
  
} else if (study.code == "RS_") {
  
  fast.rates.in = 'rs_fast_rate_calculations.xlsx'
  
  ## For .txt files for image times
  input.path = paste0("C:/Users/",pnnl.user,"/OneDrive - PNNL/Shared Documents - Core Richland and Sequim Lab-Field Team/Data Generation and Files/Cincinnati_RS/Optodes")
  
  ## For mapping files and result .csv's
  map.path = paste0("C:/Users/",pnnl.user,"/OneDrive - PNNL/Shared Documents - Core Richland and Sequim Lab-Field Team/Data Generation and Files/Cincinnati_RS/Optodes")
  
  ## Where to put Final Raw DO file
  output.path = paste0("C:/Users/",pnnl.user,"/OneDrive - PNNL/Shared Documents - Core Richland and Sequim Lab-Field Team/Data Generation and Files/Cincinnati_RS/Optodes/Rate_Processing/Raw_DO_by_Min")
  
} else if (study.code == "EL_") {
  
  fast.rates.in = 'el_fast_rate_calculations.xlsx'
  
  ## For .txt files for image times
  input.path = ("Y:/Optode_multi_reactor/Optode_multi_reactor_incubation/")
  
  ## For mapping files and result .csv's
  map.path = paste0("C:/Users/",pnnl.user,"/OneDrive - PNNL/Shared Documents - Core Richland and Sequim Lab-Field Team/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/")
  
  ## Where to put Final Raw DO file
  output.path = paste0("C:/Users/",pnnl.user,"/OneDrive - PNNL/Shared Documents - Core Richland and Sequim Lab-Field Team/Data Generation and Files/ECA/INC/03_ProcessedData/Raw_DO_by_Min/")
  
  
}


# Read in 100% saturation values for each kit based on pressure/temperature during disk calibration
fast.sat <- read_excel(paste0(map.path, "/rates/", fast.rates.in)) %>% 
  rename("DO_mg_L" = "DO_sat_mg_L") 

#read in respiration data and clean
import_data = function(map.path){
  
  # pull a list of file names in the target folder with the target pattern
  # then read all files and combine
  
  filePaths <- list.files(path = map.path, recursive = T, pattern = "\\.csv$", full.names = TRUE)
  
  filePaths <- filePaths[grepl(study.code, filePaths)]
  
  filePaths <- filePaths[grepl("results", filePaths)]
  
  # Remove ECA samples incubated in Jars. These were done at the beginning of the ECA 2022 experiment, before moving to 50 mL tubes because samples were not drying out sufficiently 
  
  if (study.code == 'EC') {
  filePaths <- filePaths[!grepl("EC_01|EC_02|EC_03|EC_04_08|EC_06_07|EC_10_15", filePaths)] }
  
    # dat <- 
  do.call(rbind, lapply(filePaths, function(map.path){
    # then add a new column `source` to denote the file name
    df <- read.csv(map.path, skip = 4, stringsAsFactors =  F)

    df[["source_file"]] <- rep(map.path, nrow(df)) # add a column for the source file
    
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
      mutate(source_file = str_remove_all(source_file, paste0(map.path, "/")))
  }
  ))
}
data = import_data(map.path)

##### Clean Data ####

if (grepl("EC|EV|EL", study.code)) {
  # Put in long form
data_long = 
  data %>% 
  mutate(disc_number = str_remove_all(disc_number, "X")) %>%
  mutate(source_file = str_remove_all(source_file, ".*Optode_multi_reactor_incubation//"),source_file = str_remove_all(source_file, "/results/results.csv")) %>%   filter(elapsed_min > 0) } else if (study.code == "RS") {
    # Put in long form
    data_long = 
      data %>% 
      mutate(disc_number = str_remove_all(disc_number, "X")) %>%
      mutate(source_file = str_remove_all(source_file, ".*Optodes/"),source_file = str_remove_all(source_file, "/results/results.csv")) %>%  #use this for RS
      filter(elapsed_min > 0) }

## Pull in mapping sheets
import_data = function(map.path){ 
  
  #pull a list of files in target folder with correct pattern
  #read all files and combine
  
  map.file <-  list.files(map.path, recursive = T, pattern = "\\.xlsx$",full.names = T)
  
  map.file <- map.file[grepl(study.code, map.file)]
  
  if (study.code == 'EC') {
  # Remove samples incubated in jars for EC samples
  map.file <- map.file[!grepl("Green|EC_01_|EC_02_|EC_03_|EC_06_07|EC_10_15|EC_04_08", map.file)] }
  
  mapping <- lapply(map.file, read_xlsx)
  
  for (i in 1:length(mapping)){mapping[[i]] <- cbind(mapping[[i]], map.file[i])}
  
  all.map <- 
    do.call(rbind,mapping)
}

all.map = import_data(map.path)

all.map$`Time on` <- as.POSIXct(all.map$`Time on`, format = "%Y/%m/%d %H:%M:%%S")

all.map$`Time on` <- format(all.map$`Time on`, format = "%H:%M")

if (grepl("EC|EV|EL", study.code)) {

all.map.clean = all.map %>% 
  rename("source_file" = "map.file[i]") %>% 
  rename("disc_number" = "Disk_ID") %>% 
  mutate(source_file = str_remove_all(source_file, ".*//")) %>% # use this for EC/EV
  separate(source_file, sep = "/", c("source_file", "file")) %>% 
  dplyr::select(-`Disk Calibration date`, -Location, -`Time off`, -SpC, -pH, -Temp, -Notes, -file) } else if(study.code == 'RS') {
    all.map.clean = all.map %>% 
    rename("source_file" = "map.file[i]") %>% 
    rename("disc_number" = "Disk_ID") %>% 
    mutate(source_file = str_remove_all(source_file, ".*Optodes/")) %>% # use this for RS
    separate(source_file, sep = "/", c("source_file", "file")) %>% 
    dplyr::select(-`Disk Calibration date`, -Location, -`Time off`, -SpC, -pH, -Temp, -Notes, -file)}

all.samples <- merge(data_long, all.map.clean)

# Fill in samples, remove last two columns for cleaner data frame
bind <- merge(all.samples, fast.sat, all = TRUE) %>% 
  filter(!grepl("Blank", Sample_ID)) 

cleaned_data <- bind %>% 
  relocate(DO_mg_L, .before = elapsed_min) %>% 
  mutate(DO_mg_per_L = DO_mg_L) %>% 
  dplyr::select(c(Sample_ID, DO_mg_per_L, LOD_mg_L, elapsed_min, `Time on`)) %>% 
  arrange(Sample_ID) %>% 
  mutate(Methods_Deviation = "N/A") %>% 
  mutate(Methods_Deviation = ifelse(elapsed_min == 0, "RATE_004", "N/A")) %>% 
  rename(Elapsed_Minute = elapsed_min) %>% 
  rename(Sample_Name = Sample_ID) 
  


# Pull in .txt files with picture times and combine
import_txt = function(input.path){ 
  
  filePaths <- list.files(path = input.path, recursive = T, pattern = "\\.txt$", full.names = TRUE)
  
  filePaths <- filePaths[grepl(study.code, filePaths)]
  filePaths <- filePaths[!grepl("cal|images", filePaths)]
  #filePaths <- filePaths[!grepl("images", filePaths)]
  
  mapping <- lapply(filePaths, read.delim)
  
  for (i in 1:length(mapping)){mapping[[i]] <- cbind(mapping[[i]], filePaths[i])}
  
  all.txt <- 
    do.call(rbind,mapping)
}

all.txt = import_txt(input.path)  

img_all <- all.txt %>% 
  mutate(`filePaths[i]` = str_remove_all(`filePaths[i]`, paste0(input.path, "/"))) %>% 
  separate(col = `filePaths[i]`, into = c("source_file", "photo"), sep = "/")

# pull out times and dates of images
img_time <- img_all[!grepl("Custom|type|.tif|----", img_all$TIFF.image.set.saved.with.Look.RGB.v0.1),] 

# rename date/time column as Time
img_time <- rename(img_time, Time = TIFF.image.set.saved.with.Look.RGB.v0.1)

# remove AM/PM in Time
img_time$Time <- gsub('[AMP]','', img_time$Time)

# make time as.POSIXct format
img_time$Time <- as.POSIXct(img_time$Time, format = "%m/%d/%Y %H:%M:%S")

# separate M/D/Y from time
img_time$Day <- format(img_time$Time, format = "%m/%d/%Y")

# make M/D/Y a date
img_time$Day <- as.Date(img_time$Day, format = "%m/%d/%Y")

# Correct for daylight savings, imageJ and mapping file times won't match during these times
img_time$Time_Corr <- ifelse(img_time$Day < "2022-11-06" | img_time$Day > "2023-03-12", img_time$Time + 3600, img_time$Time)

# correct to POSIXct
class(img_time$Time_Corr) <- c("POSIXct", "POSIXt")

# format Time as H:M:S
img_time$Time_HMS <- format(as.POSIXct(img_time$Time_Corr), format = "%H:%M:%S")

# Pull out H:M - this as the best resolution we get from mapping file
img_time$Time_HM <- format(as.POSIXct(img_time$Time_Corr), format = "%H:%M")

# pull out second picture is taken - use this if picture was taken early on in minute
img_time$Time_S <- format(as.POSIXct(img_time$Time_Corr), format = "%S")

#merge time of pictures with 'Time on' from mapping files
all.times <- merge(img_time, all.map.clean) %>% 
  drop_na(`Time on`)

# flag what the minute before would be
all.times$min_bef <- format(strptime(all.times$Time_HM, format = "%H:%M") - 60, "%H:%M")

all.times$time_same <- NA

# Flag if samples are put on at same time a picture is taken
for (i in 1:nrow(all.times)) {
  
  if (all.times$Time_HM[i] == all.times$`Time on`[i]) {
    
    all.times$time_same[i] <- "yes"
    
  }
  
  else if (all.times$Time_S[i] <= 20 & all.times$`Time on`[i] == all.times$min_bef[i]) {
    
    # if picture taken < 20 s after new minute, sample was maybe put on near start
    all.times$time_same[i] <- "maybe"
    
  }
  
  else { 
    
    all.times$time_same[i] <- "no"
    
  }
  
}

corr.time <- all.times %>% 
  dplyr::select(c(source_file, Sample_ID, disc_number, Time_HMS, Time_HM, `Time on`, time_same)) %>% 
  filter(!grepl("Blank", Sample_ID)) %>% 
  group_by(Sample_ID) %>% 
  filter(Time_HM >= `Time on`) %>% 
  arrange(Time_HM) %>% 
  filter(row_number() == 1) %>% 
  rename(Sample_Name = Sample_ID) %>% 
  ungroup() %>% 
  mutate(Elapsed_Minute = 2)

samples <- merge(corr.time, cleaned_data, by = c("Sample_Name", "Elapsed_Minute"), all = TRUE) %>% 
  mutate(Methods_Deviation = if_else(is.na(time_same), Methods_Deviation, if_else(time_same == "no", Methods_Deviation, if_else(Methods_Deviation == "N/A", "RATE_006", paste0(Methods_Deviation, "; RATE_006"))))) %>% 
  dplyr::select(-c(source_file, Time_HMS, Time_HM, disc_number, `Time on.x`, `Time on.y`, time_same)) %>% 
  separate(Sample_Name, into = c("EC", "kit", "rep"), remove = FALSE, sep = "_")

if (study.code == 'EC'){
# Add extra 0 in EC Sample Names
for (i in 1:nrow(samples)){
  
  if (str_count(samples$kit[i], "[0-9]") <= 2){
    
    samples$kit[i] = paste0("0", samples$kit[i])
    
  }
  
  else {
    
    samples$kit[i] = samples$kit[i]
  }
  
}
}

## Don't need for RS samples
#EC 011/012 - W: Add INC_008 deviation
#EC 13, 14, 27: Flag overexposed replicates
#EC 012 D5: Less sediment in replicate
# EC 72 - D5/W5: missing replicates

if (study.code == 'EC') {

missing_reps <- data.frame(
  Sample_Name = c("EC_072_INC-D5", "EC_072_INC-W5"),
  Elapsed_Minute = c(-9999, -9999),
  DO_mg_per_L = c(-9999, -9999),
  LOD_mg_L = c(-9999, -9999), 
  Methods_Deviation = c("INC_Method_001", "INC_Method_001"),
  stringsAsFactors = FALSE  # Ensure that strings are not converted to factors
)

samples_clean = samples %>% 
  unite(Sample_Name, c("EC", "kit", "rep"), sep = "_") %>%   bind_rows(missing_reps) %>% 
  arrange(Sample_Name) %>% 
  mutate(DO_mg_per_L = round(DO_mg_per_L, 2)) %>% # round to 2 for data package publishing
  relocate(Sample_Name, .before = Elapsed_Minute) %>% 
  group_by(Sample_Name) %>% 
  fill(LOD_mg_L) %>% 
  ungroup() %>% 
  mutate(Methods_Deviation = if_else(grepl("EC_011_INC-W|EC_012_INC-W", Sample_Name), if_else(Methods_Deviation == "N/A", "INC_008", paste0(Methods_Deviation, "; INC_008")), Methods_Deviation)) %>% 
  mutate(Methods_Deviation = if_else(grepl("EC_011_INC-W|EC_012_INC-W", Sample_Name), if_else(Methods_Deviation == "N/A", "INC_008", paste0(Methods_Deviation, "; INC_008")), Methods_Deviation)) %>% 
  mutate(Methods_Deviation = if_else(grepl("EC_012_INC-D5", Sample_Name), if_else(Methods_Deviation == "N/A", "INC_Method_002", paste0(Methods_Deviation, "; INC_Method_002")), Methods_Deviation)) %>% 
  filter(!(Sample_Name %in% c("EC_013_INC-D4", "EC_013_INC-W5", "EC_014_INC-W5", "EC_027_INC-D1", "EC_027_INC-D2", "EC_027_INC-W5") & Elapsed_Minute > 0)) %>% 
  mutate(Elapsed_Minute = if_else(grepl("EC_013_INC-D4|EC_013_INC-W5|EC_014_INC-W5|EC_027_INC-D1|EC_027_INC-D2|EC_027_INC-W5", Sample_Name), -9999, Elapsed_Minute)) %>% 
    mutate(DO_mg_per_L = if_else(grepl("EC_013_INC-D4|EC_013_INC-W5|EC_014_INC-W5|EC_027_INC-D1|EC_027_INC-D2|EC_027_INC-W5", Sample_Name), -9999, DO_mg_per_L)) %>%
  mutate(Methods_Deviation = if_else(grepl("EC_013_INC-D4|EC_013_INC-W5|EC_014_INC-W5|EC_027_INC-D1|EC_027_INC-D2|EC_027_INC-W5", Sample_Name), if_else(Methods_Deviation == "N/A", "INC_QA_004",paste0(Methods_Deviation, "; INC_QA_004")), Methods_Deviation)) %>% 
  mutate(Methods_Deviation = if_else(Sample_Name == "EC_091_INC-D5" & DO_mg_per_L == -9999, if_else(Methods_Deviation == "N/A", "INC_002", paste0(Methods_Deviation, "; INC_002")), Methods_Deviation)) %>% 
mutate(Methods_Deviation = if_else(grepl("EC_022_INC-D5|EC_044_INC-D4|EC_044_INC-D5", Sample_Name), if_else(Methods_Deviation == "N/A", "INC_003", paste0(Methods_Deviation, "; INC_003")), Methods_Deviation)) %>% 
  mutate(Methods_Deviation = if_else(grepl("EC_042_INC-D5|EC_064_INC-W1|EC_032_INC-D5", Sample_Name), if_else(Methods_Deviation == "N/A", "INC_004", paste0(Methods_Deviation, "; INC_004")), Methods_Deviation)) %>% 
  mutate(Methods_Deviation = if_else(grepl("EC_041_INC-W4", Sample_Name), if_else(Methods_Deviation == "N/A", "INC_005", paste0(Methods_Deviation, "; INC_005")), Methods_Deviation)) %>% 
  mutate(Methods_Deviation = if_else(grepl("EC_071", Sample_Name), if_else(Methods_Deviation == "N/A", "VI_ALL_006", paste0(Methods_Deviation, "; VI_ALL_006")), Methods_Deviation)) %>% 
  mutate(Methods_Deviation = if_else(grepl("EC_069_INC-D2", Sample_Name), if_else(Methods_Deviation == "N/A", "VI_ALL_007", paste0(Methods_Deviation, "; VI_ALL_007")), Methods_Deviation)) %>% 
  mutate(LOD_mg_L = round(LOD_mg_L, 3)) %>% 
  mutate(DO_mg_per_L_corr = as.character(DO_mg_per_L)) %>%
  mutate(LOD_mg_L = as.character(LOD_mg_L)) %>% 
  mutate(DO_mg_per_L_corr = if_else(DO_mg_per_L_corr <= 0, paste0("DO_Below_", LOD_mg_L,"_mg_per_L_LOD|Negative_Value_Raw_Not_Corrected|-9999_mg_per_L_Final_Corrected"), DO_mg_per_L_corr)) %>% 
  mutate(DO_mg_per_L_corr = if_else(DO_mg_per_L_corr <= LOD_mg_L, paste0("DO_Below_", LOD_mg_L,"_mg_per_L_LOD|", DO_mg_per_L_corr, "_mg_per_L_Raw_Not_Corrected|", DO_mg_per_L_corr,"mg_per_L_Final_Corrected"), DO_mg_per_L_corr)) %>% 
  mutate(Methods_Deviation = if_else(grepl("Raw", DO_mg_per_L_corr), if_else(Methods_Deviation == "N/A", "DTL_003", paste0(Methods_Deviation, "; DTL_003")), Methods_Deviation)) %>% 
  select(c(Sample_Name, DO_mg_per_L_corr, Elapsed_Minute, Methods_Deviation)) %>% 
  rename(DO_mg_per_L = DO_mg_per_L_corr) } else {
    
    samples_clean = samples %>% 
      unite(Sample_Name, c("EC", "kit", "rep"), sep = "_") %>%  # bind_rows(missing_reps) %>% 
      arrange(Sample_Name) %>% 
      mutate(DO_mg_per_L = round(DO_mg_per_L, 2)) %>% # round to 2 for data package publishing
      relocate(Sample_Name, .before = Elapsed_Minute) %>% 
      group_by(Sample_Name) %>% 
      fill(LOD_mg_L) %>% 
      ungroup() %>% 
      mutate(LOD_mg_L = round(LOD_mg_L, 3)) %>% 
      mutate(DO_mg_per_L_corr = as.character(DO_mg_per_L)) %>%
      mutate(LOD_mg_L = as.character(LOD_mg_L)) %>% 
      mutate(DO_mg_per_L_corr = if_else(DO_mg_per_L_corr <= 0, paste0("DO_Below_", LOD_mg_L,"_mg_per_L_LOD|Negative_Value_Raw_Not_Corrected|-9999_mg_per_L_Final_Corrected"), DO_mg_per_L_corr)) %>% 
      mutate(DO_mg_per_L_corr = if_else(DO_mg_per_L_corr <= LOD_mg_L, paste0("DO_Below_", LOD_mg_L,"_mg_per_L_LOD|", DO_mg_per_L_corr, "_mg_per_L_Raw_Not_Corrected|", DO_mg_per_L_corr,"mg_per_L_Final_Corrected"), DO_mg_per_L_corr)) %>% 
      mutate(Methods_Deviation = if_else(grepl("Raw", DO_mg_per_L_corr), if_else(Methods_Deviation == "N/A", "DTL_003", paste0(Methods_Deviation, "; DTL_003")), Methods_Deviation)) %>% 
      select(c(Sample_Name, DO_mg_per_L_corr, Elapsed_Minute, Methods_Deviation)) %>% 
      rename(DO_mg_per_L = DO_mg_per_L_corr)
    
  } 

output.name = paste0(study.code,"INC_Raw_DO_By_Min_ReadyForBoye_",Sys.Date(),".csv")

write.csv(samples_clean, file.path(output.path, output.name), quote = F, row.names = F) 
