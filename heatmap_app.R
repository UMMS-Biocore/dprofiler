# Dprofiler Shiny App

# library
library(shiny)
library(shinyjs)
library(shinydashboard)
library(debrowser)
library(heatmaply)
library(gplots)

source("functions/heatmap.R")
source("functions/heatmapfunctions.R")

shinyApp(ui = shinyUI(heatmapUI), server = shinyServer(heatmapServer))