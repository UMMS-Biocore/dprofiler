# Dprofiler Shiny App

# library
library(shiny)
library(shinyjs)
library(shinydashboard)
library(debrowser)
library(cluster)
library(xbioc)
library(MuSiC)
library(VennDiagram)
library(SignallingSingleCell)
library(heatmaply)
library(gplots)

# source
source("ui.R")
source("server.R")
source("functions/dataLoad.R")
source("functions/condSelect.R") 
source("functions/iterdeprogs.R")
source("functions/deprogs.R")
source("functions/scoring.R")
source("functions/funcs.R")
source("functions/deconvolute.R")
source("functions/mainScatter.R")
source("functions/heatmap.R")
source("functions/Bar.R")

# Run the application 
shinyApp(ui = shinyUI(DprofilerUI),
         server = shinyServer(DprofilerServer))

