## combine plots from Sensitivity Analysis Outputs
library(pdftools)

setwd("C:/GitHub/Cincinnati_Multireactor_Respiration/Data/Plots/Linear_Fits_Rules/")

ex_pdf = list.files(pattern = "pdf")


ex_pdf1 <- head(ex_pdf, 275)
ex_pdf2 <- tail(ex_pdf, 278)

pdf_combine(ex_pdf1, output = "First_Half.pdf")

pdf_combine(ex_pdf2,  output = "Second_Half.pdf")
 
all_pdf <- list.files(pattern = "Half.pdf")
pdf_combine(ex_pdf, output = "RS_All_Combined_Fits.pdf")
