library(patchwork)


pnnl.user = 'laan208'

output.path = paste0("C:/Users/",pnnl.user,"/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/INC/03_ProcessedData/")

drying <- read.csv("C:/Users/laan208/PNNL/Core Richland and Sequim Lab-Field Team - Documents/Data Generation and Files/ECA/INC/03_ProcessedData/ECA_Drying_Masses_merged_by_laan208_on_2023-09-15.csv")

drying <- drying %>% 
  filter(!grepl("EC_023", Sample_Name))

unique.incubations = unique(drying$Sample_Name)

my.formula = y ~ x


for (i in 1:length(unique.incubations)) {
  
  data_subset = subset(drying, drying$Sample_Name == unique.incubations[i])
  
  data_subset$grav_moisture = data_subset$Water_total_g/data_subset$Dry_weight_sed_g
  
  data_subset$Date <- as.Date(data_subset$Date, format='%Y-%m-%d')
  
  fit_grav = lm(data_subset$grav_moisture~data_subset$Date)
  u = fit_grav$coefficients
  b = u[[1]] #Intercept
  c = u[[2]] #slope
  r = summary(fit_grav)$r.squared
  r.adj = summary(fit_grav)$adj.r.squared
  residuals = deviance(fit_grav)
  p = summary(fit_grav)$coefficients[4]
  r2 = r
  
  fit_wat = lm(data_subset$Water_total_g~data_subset$Date)
  u = fit_wat$coefficients
  b = u[[1]] #Intercept
  c = u[[2]] #slope
  r = summary(fit_wat)$r.squared
  r.adj = summary(fit_wat)$adj.r.squared
  residuals = deviance(fit_wat)
  p = summary(fit_wat)$coefficients[4]
  r2 = r
  
  grav <- ggplot(data_subset, aes(x = Date, y = grav_moisture))+
    geom_point()+
    coord_cartesian(ylim = c(0,1.5))+ 
    geom_smooth(method = "lm", se=F, formula = my.formula)+
    geom_label(aes(x = as.Date(first(Date)), y = 1), size = 2.5, hjust = 0,
               label = paste("R2 = ", signif(summary(fit_grav)$r.squared, 3),
                             "\nP = ", signif(summary(fit_grav)$coefficients[[4]], 3),
                             "\nSlope = ", signif(fit_grav$coefficients[[2]], 3)))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    xlab("Gravimetric Moisture")+
    ggtitle(data_subset$Sample_Name[1])
    
  
  water <- ggplot(data_subset, aes(x = Date, y = Water_total_g))+
    geom_point()+
    coord_cartesian(ylim = c(0,15))+ 
    geom_smooth(method = "lm", se=F, formula = my.formula)+
    geom_label(aes(x = as.Date(first(Date)), y = 12), size = 2.5, hjust = 0,
               label = paste("R2 = ", signif(summary(fit_grav)$r.squared, 3),
                             "\nP = ", signif(summary(fit_grav)$coefficients[[4]], 3),
                             "\nSlope = ", signif(fit_grav$coefficients[[2]], 3)))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    xlab("Total Water (g)")+
    ggtitle(data_subset$Sample_Name[1])
    
  
  multi <- (grav + water) +
    plot_layout(widths = c(1,1))
  
  ggsave(file=paste0(output.path,"Plots/ECA_dry_down_",data_subset$Sample_Name[1],".pdf"), width = 7, height = 7, units = "in")
  
}

library(pdftools)

setwd(paste0(output.path,"/Plots/"))

ex_pdf <- list.files(pattern = "pdf")

ex_pdf1 <- head(ex_pdf, 275)
ex_pdf2 <- tail(ex_pdf, 278)

pdf_combine(ex_pdf1, output = "First_Half.pdf")
pdf_combine(ex_pdf2, output = "Second_Half.pdf")

all_pdf <- list.files(pattern = "Half.pdf")
pdf_combine(all_pdf, output = "Dry_Down_Fits.pdf")



