#' dprofilerUI
#'
#' Creates a shinyUI to be able to run DEBrowser interactively.
#'
#' @note \code{dprofilerUI}
#' @return the panel for main plots;
#'
#' @examples
#'     x<-dprofilerUI()
#'     
#' @export
#' 
dprofilerUI <- function() {
  dbHeader <- shinydashboard::dashboardHeader(titleWidth = 250)
  dbHeader$children[[2]]$children <- tags$a(style='color: white;',
                                            id="top_logo" , "Dprofiler v1.0.0")
  addResourcePath(prefix = "www", directoryPath = system.file("extdata",
                                                              "www", package = "dprofiler"))
  dprofiler <- (fluidPage(
    shinyjs::useShinyjs(),
    tags$head(tags$title("Dprofiler v1.0.0"),
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
                                                  menuSubItem("What is Dprofiler?", tabName = "dprofilerintro"),
                                                  menuSubItem("Data Upload", tabName = "Intro"),
                                                  menuSubItem("Data Preprocessing", tabName = "assesment"),
                                                  menuSubItem("Computational Profiling", tabName = "heteroanalysis"),
                                                  menuSubItem("Cellular Comp. Analysis", tabName = "cellularanalysis"),
                                                  menuSubItem("Comparative Profiling", tabName = "comparativeanalysis"),
                                                  menuSubItem("FAQ", tabName ="FAQ")
                                         ),
                                         menuItem("Data Upload", icon = icon("upload"), tabName = "Upload", selected = TRUE),
                                         menuItem("Data Processing", icon = icon("filter"), tabName = "DataProcessing"),
                                         menuItem("Computational Profiling", icon = icon("dna"), tabName = "DEAnalysis"),
                                         menuItem("Cellular Comp. Analysis", icon = icon("chart-bar"), tabName = "CellComp"),
                                         menuItem("Comparative Profiling", icon = icon("project-diagram"), tabName = "Profile")
                                         # menuItem("DEFilter",  icon = icon("code"), tabName = "CondSelect",  startExpanded = TRUE,
                                         #          uiOutput("cutOffUI")),
                                         # menuItem("ScoreFilter",  icon = icon("code"), tabName = "CondSelect",  startExpanded = TRUE,
                                         #          uiOutput("ScoreCutOffUI"))
                             ),
                             helpText("Developed by ", a("UMMS Biocore.", href="https://www.umassmed.edu/biocore/", target = "_blank"))
                    ),
                    tabPanel(title = "Discover", value = "discover", id="discover",
                             mainPlotControlsUI("deresults"),
                             menuItem("DEFilter",  startExpanded = TRUE,
                                      # uiOutput("cutOffUI"),
                                      cutOffSelectionUI("deresults")),
                             menuItem("ScoreFilter", startExpanded = TRUE,
                                      # uiOutput("ScoreCutOffUI"),
                                      ScoreCutOffSelectionUI("deresults")
                                      ),
                             shinydashboard::menuItem("DE Heatmaps", 
                                                      heatmapControlsUI("deresults")),
                             shinydashboard::menuItem("Deconvolution Heatmaps", 
                                                      heatmapControlsUI("deconvolute"))
                    )
        )
      ),
      
      # Shiny dashboard Body
      dashboardBody(
          
          tabItems(
                    
                   # Help Tab
                   tabItem(tabName="dprofilerintro", getDprofilerText()),
                   tabItem(tabName="Intro", getIntroText()),
                   tabItem(tabName="assesment", getDataAssesmentText()),
                   tabItem(tabName="heteroanalysis", getCompProfilingText()),
                   tabItem(tabName="cellularanalysis", getCompCellularText()),
                   tabItem(tabName="comparativeanalysis", getComparativeProfText()),
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
                                  ),
                                  tabPanel(title = "Upload Summary (Reference Bulk)",
                                           dataProfileSummaryUI("load"),
                                           value = "uploadprofilesummary"
                                  )
                           )
                   ),

                   # Upload Tab
                   tabItem(tabName="DataProcessing", 
                           tabBox(id = "DataProcessingBox", 
                                  width = NULL,
                                  tabPanel(title = "Filtering",
                                           conditionalPanel((condition <- "input.Filter"),
                                                            dataLCFUI("lcf")),
                                           value = "filter"
                                  ),
                                  tabPanel(title = "Batch Effect Correction",
                                           conditionalPanel((condition <- "input.Batch"),
                                                            batchEffectUI("batcheffect")),
                                           value = "batcheffect"
                                  )
                           )
                   ),
                   
                   # # Condition Selection Tab
                   # tabItem(tabName="CondSelect", 
                   #         condSelectUI()
                   # ),
                   
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