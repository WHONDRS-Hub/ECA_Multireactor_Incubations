library(readxl);library(tidyverse)

#set working directory
pnnl.user = 'laan208'


##update this to have most recent date
effect <- read.csv(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Effect_Size_Data/Effect_Size_merged_by_laan208_on_2023-05-25.csv"))

effect = effect[,c("Site", "effect")]
effect = effect[!duplicated(effect$effect), ]

effect$Sample_ID = substr(effect$Sample_ID,1,nchar(effect$Sample_ID)-7)

effect$Site[effect$Site == "5"] <- "05"
effect$Site[effect$Site == "9"] <- "09"

sites <- read_xlsx(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Site Coordinates.xlsx"))

sites <- sites %>% 
  filter(Type != "jar") %>% 
  separate(Site, into = c("EC", "Site"), sep = "_")

join <- merge(sites, effect, by = "Site")


limits = c(-4,
           max(c(join$effect)))

join$Longitude <- as.numeric(join$Longitude)
join$Latitude <- as.numeric(join$Latitude)

png(file = paste0("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/", as.character(Sys.Date()),"_Effect_Size_Map.png"), width = 10, height = 10, units = "in", res = 300)

effect.plot = ggplot(data = join)+
  borders("state", colour = "black", size = 0.3)+ 
  geom_point(aes(x = Longitude, y = Latitude, color = effect), size = 3.5)+
   scale_color_gradient2(limits = limits, low = "firebrick2", mid = "goldenrod2",
                        high = "dodgerblue2", midpoint = (max(limits)+min(limits))/2)+
  ggtitle(paste("Effect Size")) + labs(color = var)+
  labs(color = paste0("Wet - Dry Rate"))+
  coord_fixed() + theme_bw() + theme(axis.line=element_blank(),
                                     axis.text.x=element_blank(),
                                     axis.text.y=element_blank(),
                                     axis.ticks=element_blank(),
                                     axis.title.x=element_blank(),
                                     axis.title.y=element_blank(),
                                     panel.background=element_blank(),
                                     panel.border=element_blank(),
                                     panel.grid.major=element_blank(),
                                     panel.grid.minor=element_blank(),
                                     plot.background=element_blank())+
  theme(aspect.ratio=6/10)

print(effect.plot)

dev.off()

## Mud Map ####

grain <- read_csv(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ICON_ModEx_SSS/09_Grain_Size/03_ProcessedData/20221230_Grain_Size_SBR_RC4_CM_1-42/20221230_Data_Processed_Grain_Size_SBR_RC4_CM_1-42.csv"))

grn <- grain %>% 
  mutate(Percent_Mud = Percent_Silt + Percent_Clay) %>% 
  separate(Sample_ID, into = c("CM", "Site", "GRN"), sep = "_")

sub <- 2
grn$Site <- substr(grn$Site, nchar(grn$Site) - sub + 1, nchar(grn$Site))


grn_map <- merge(sites, grn, by = "Site")

limits = c(0,
           max(c(grn_map$Percent_Mud)))

grn_map$Longitude <- as.numeric(grn_map$Longitude)
grn_map$Latitude <- as.numeric(grn_map$Latitude)

png(file = paste0("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/Figures/", as.character(Sys.Date()),"_Percent_Mud_Map.png"), width = 10, height = 10, units = "in", res = 300)

ggplot(data = grn_map)+
  borders("state", colour = "black", size = 0.3)+ 
  geom_point(aes(x = Longitude, y = Latitude, color = Percent_Mud), size = 3.5)+
  scale_color_gradient2(limits = limits, low = "dodgerblue2", mid = "goldenrod2",
                        high = "firebrick2", midpoint = (max(limits)+min(limits))/2)+
  ggtitle(paste("Percent Mud")) + labs(color = var)+
  labs(color = paste0("Percent_Mud"))+
  coord_fixed() + theme_bw() + theme(axis.line=element_blank(),
                                     axis.text.x=element_blank(),
                                     axis.text.y=element_blank(),
                                     axis.ticks=element_blank(),
                                     axis.title.x=element_blank(),
                                     axis.title.y=element_blank(),
                                     panel.background=element_blank(),
                                     panel.border=element_blank(),
                                     panel.grid.major=element_blank(),
                                     panel.grid.minor=element_blank(),
                                     plot.background=element_blank())+
  theme(aspect.ratio=6/10)

dev.off()

