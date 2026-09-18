#------------------------------------------------------------------------------#
#publish.R
#
#Author: 
#Sindre Vatnehol
#Institute of marine Research
#Norway
#
#In this script we produce the different presentation for the assessment and 
#forecast. 
#
#
#------------------------------------------------------------------------------#





#------------------------------------------------------------------------------#
# Remove everything from memory to avoid conflicts
#------------------------------------------------------------------------------#
rm(list=ls())





#------------------------------------------------------------------------------#
#Load library, and source files
#------------------------------------------------------------------------------#
library(icesTAF)
library(ggplot2)
library('stockassessment')
library(quarto)
library(ggplotify)   
library(patchwork)  
library(ggrepel)
library(gridExtra)




#------------------------------------------------------------------------------#
#make the folder structure
#------------------------------------------------------------------------------#
mkdir('publish')
mkdir('publish/presentation')




#------------------------------------------------------------------------------#
#Render the assessment presentation. 
#------------------------------------------------------------------------------#
# reproducible <- Sys.getenv("REPRODUCIBLE", unset = "true")


file.copy("boot/quarto/assessment_presentation.qmd", 
          "assessment_presentation.qmd",overwrite = T )
quarto::quarto_render(
  input       = "assessment_presentation.qmd",
  output_file = "assessment_presentation.pdf"
)
unlink("assessment_presentation.qmd")
file.rename("assessment_presentation.pdf", 
            "publish/presentation/assessment_presentation.pdf")


#------------------------------------------------------------------------------#
#Render the forecast presentation.
#------------------------------------------------------------------------------#
file.copy("boot/quarto/forecast_presentation.qmd",
          "forecast_presentation.qmd",overwrite = T)
quarto::quarto_render(
  input       = "forecast_presentation.qmd",
  output_file = "forecast_presentation.pdf"
)
unlink("forecast_presentation.qmd")
file.rename("forecast_presentation.pdf",
            "publish/presentation/forecast_presentation.pdf")
