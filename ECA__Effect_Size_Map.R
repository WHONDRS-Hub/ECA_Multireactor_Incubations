library(readxl);library(tidyverse)

#set working directory
pnnl.user = 'laan208'


##update this to have most recent date
effect <- read.csv(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/rates/Effect_Size_Merged_by_laan208_on_2023-04-26.csv"))

effect = effect[,c("Sample_ID", "pos_effect")]
effect = effect[!duplicated(effect$pos_effect), ]

effect$Sample_ID = substr(effect$Sample_ID,1,nchar(effect$Sample_ID)-7)

sites <- read_xlsx(paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Site Coordinates.xlsx"))

sites <- sites %>% 
  filter(Type != "jar") %>% 
  filter(Site != "EC_27") %>%
  filter(Site != "EC_39")

colnames(sites)[1] = "Sample_ID"

join <- merge(sites, effect, by = "Sample_ID")


limits = c(-4,
           max(c(join$pos_effect)))


png(file = paste0("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/Optode multi reactor/Optode_multi_reactor_incubation/effect size/ESS-PI_EGU/", as.character(Sys.Date()),"_Effect_Size_Map.png"), width = 10, height = 10, units = "in", res = 300)

effect.plot = ggplot(data = join)+
  borders("state", colour = "black", size = 0.3)+ 
  geom_point(aes(x = Longitude, y = Latitude, color = pos_effect), size = 3.5)+
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
