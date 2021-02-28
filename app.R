# Dprofiler Shiny App

# library
library(shiny)
library(shinyjs)
library(shinydashboard)
library(debrowser)
library(cluster)
library(MuSiC)
library(VennDiagram)

# source
source("ui.R")
source("server.R")
source("functions/dataLoad.R")
# source("functions/dataLCFUI.R")
# source("functions/batcheffect.R")
source("functions/condSelect.R") 
source("functions/iterdeprogs.R")
# source("functions/funcs.R") 
source("functions/deprogs.R")
source("functions/silhouette.R")

# Run the application 
shinyApp(ui = shinyUI(DprofilerUI),
         server = shinyServer(DprofilerServer))
  