# Dconnect Shiny App

# library
library(shiny)
library(shinyjs)
library(shinydashboard)
library(debrowser)

# source
source("ui.R")
source("server.R")
source("functions/dataLoad.R")
source("functions/dataLCFUI.R")
source("functions/batcheffect.R")
source("functions/funcs.R")

# Run the application 
shinyApp(ui = shinyUI(dconnectUI), 
         server = shinyServer(dconnectServer))
