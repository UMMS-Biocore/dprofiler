#' deUI
#'
#' Creates a shinyUI to be able to run DEBrowser interactively.
#'
#' @note \code{deUI}
#' @return the panel for main plots;
#'
#' @examples
#'     x<-deUI()
#'
#' @export
#'

dprofilerUI <- function() {
  dbHeader <- shinydashboard::dashboardHeader(titleWidth = 250)
  dbHeader$children[[2]]$children <- tags$a(style='color: white;',
                                            id="top_logo" , "Dprofiler")
  addResourcePath(prefix = "www", directoryPath = "www/")
  library("debrowser")
  dprofiler <- (fluidPage(
    shinyjs::useShinyjs(),
    tags$head(tags$title("Dprofiler"),
              tags$link(rel = "stylesheet", type = "text/css",
                        href = "www/shinydashboard_additional.css")
    ),
    dashboardPage(
      dbHeader,
      
      # Shiny dashboard Side
      dashboardSidebar(
        width = 250,
        getJSLine(),
        add_busy_spinner(spin = "fading-circle", position = "bottom-right"),
        use_waiter(),
        tabsetPanel(id = "menutabs", type = "tabs",
                    tabPanel(title = "Menu", value = "dataprep", id="dataprep",
                             sidebarMenu(id="MenuItems",
                                         menuItem("Quick Start Guide", icon = icon("user"),
                                                  menuSubItem("Intro. and Data Upload", tabName = "Intro"),
                                                  menuSubItem("Data Preprocessing", tabName = "assesment"),
                                                  menuSubItem("Diff. Hetero. Analysis", tabName = "heteroanalysis"),
                                                  menuSubItem("FAQ", tabName ="FAQ")
                                         ),
                                         menuItem("Upload", icon = icon("upload"), tabName = "Upload"),
                                         menuItem("Data Processing", icon = icon("filter"), tabName = "DataProcessing"),
                                         menuItem("Cond. Select", icon = icon("adjust"), tabName = "CondSelect"),
                                         menuItem("Diff. Hetero. Analysis", icon = icon("adjust"), tabName = "DEAnalysis"),
                                         menuItem("Cellular Comp.", icon = icon("adjust"), tabName = "CellComp"),
                                         menuItem("Profiling", icon = icon("adjust"), tabName = "Profile"),
                                         menuItem("DEFilter",  icon = icon("code"), tabName = "CondSelect",  startExpanded = TRUE,
                                                  uiOutput("cutOffUI")),
                                         menuItem("ScoreFilter",  icon = icon("code"), tabName = "CondSelect",  startExpanded = TRUE,
                                                  uiOutput("ScoreCutOffUI"))
                             ),
                             helpText("Developed by ", a("UMMS Biocore.", href="https://www.umassmed.edu/biocore/", target = "_blank"))
                             ),
                    tabPanel(title = "Discover", value = "discover", id="discover",
                             mainPlotControlsUI("deresults"),
                             # uiOutput('cutoffSelection'),
                             shinydashboard::menuItem("DE Heatmaps", 
                                                      heatmapControlsUI("deresults")),
                             shinydashboard::menuItem("Deconvolution Heatmaps", 
                                                      heatmapControlsUI("deconvolute"))
                             )
        )
      ),
      
      # Shiny dashboard Body
      dashboardBody(
          
          tabItems(#id = "methodtabs", type = "tabs",
                    
                   # Help Tab
                   tabItem(tabName="Intro", getIntroText()),
                   tabItem(tabName="assesment", getDataAssesmentText()),
                   tabItem(tabName="heteroanalysis", getHeteroAnalysisText()),
                   tabItem(tabName="FAQ",  getQAText()),
                   
                   # Upload Tab
                   tabItem(tabName="Upload",
                           tabBox(id = "UploadBox", 
                                  width = NULL,
                                  tabPanel(title = "Upload Data",
                                           dataLoadUI("load")
                                  ),
                                  tabPanel(title = "Upload Summary",
                                           dataSummaryUI("load"), 
                                           value = "uploadsummary"
                                  )
                           )
                   ),

                   # Upload Tab
                   tabItem(tabName="DataProcessing", 
                           tabBox(id = "DataProcessingBox", 
                                  width = NULL,
                                  tabPanel(title = "Filter",
                                           conditionalPanel((condition <- "input.Filter"),
                                                            dataLCFUI("lcf")),
                                           value = "filter"
                                  ),
                                  tabPanel(title = "BatchEffect",
                                           conditionalPanel((condition <- "input.Batch"),
                                                            batchEffectUI("batcheffect")),
                                           value = "batcheffect"
                                  )
                           )
                   ),
                   
                   # Condition Selection Tab
                   tabItem(tabName="CondSelect", 
                           condSelectUI()
                   ),
                   
                   # DE Analysis Tab
                   tabItem(tabName="DEAnalysis", 
                           uiOutput("deresUI")
                   ),
                   
                   # Deconvolution Tab
                   tabItem(tabName="CellComp", 
                           uiOutput("cellcompUI")
                   ),
                   
                   # Profiling Tab
                   tabItem(tabName="Profile", 
                           uiOutput("ProfilingUI")
                   )
          )
          
        
        # hide and show tabs as analysis progresses
        # change 'hide' part to manage order of popping up sidebar menus
        # getTabUpdateJS() 
      ))
  )
  )
  dprofiler
}