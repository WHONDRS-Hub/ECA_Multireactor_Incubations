## combine plots from Sensitivity Analysis Outputs
library(pdftools)

setwd("C:/Users/laan208/OneDrive - PNNL/Shared Documents - Core Richland and Sequim Lab-Field Team/Data Generation and Files/ECA/INC/03_ProcessedData/Plots/")

ex_pdf = list.files(pattern = "pdf")


ex_pdf1 <- head(ex_pdf, 275)
ex_pdf2 <- tail(ex_pdf, 278)

pdf_combine(ex_pdf1, output = "First_Half.pdf")

pdf_combine(ex_pdf2,  output = "Second_Half.pdf")
 
all_pdf <- list.files(pattern = "Half.pdf")
pdf_combine(all_pdf, output = "All_Fits.pdf")
