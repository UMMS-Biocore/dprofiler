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
  addResourcePath(prefix = "www", directoryPath = system.file("extdata", "www", package = "dprofiler"))
  dprofiler <- (fluidPage(
    
    # use Shinyjs
    shinyjs::useShinyjs(),
    
    # Dprofilier head and tag
    tags$head(tags$title("Dprofiler v1.0.0"),
              tags$link(rel = "stylesheet", type = "text/css",
                        href = "www/shinydashboard_additional.css")
    ),
    
    # Main Dashboard Page
    dashboardPage(
      dbHeader,
      
      # Shiny dashboard Side
      dashboardSidebar(
        width = 250,
        getJSLine(),
        add_busy_spinner(spin = "fading-circle", position = "bottom-right"),
        use_waiter(),
        tabsetPanel(id = "menutabs", type = "tabs",
                    # Mongo Server Sidebar
                    tabPanel(title = "Server", value = "dataserver", id="dataserver",
                             mongoServerUI("mongoserver")
                    ),
                    # Main Navigation Sidebar
                    tabPanel(title = "Menu", value = "dataprep", id="dataprep",
                             sidebarMenu(id="MenuItems",
                                         menuItem("Quick Start Guide", icon = icon("user", verify_fa = FALSE),
                                                  menuSubItem("What is Dprofiler?", tabName = "dprofilerintro"),
                                                  menuSubItem("Data Upload and Summary", tabName = "Intro"),
                                                  menuSubItem("Data Preprocessing", tabName = "assesment"),
                                                  menuSubItem("Comp. Phenotypic Profiling", tabName = "heteroanalysis"),
                                                  menuSubItem("Compositional Profiling", tabName = "cellularanalysis"),
                                                  menuSubItem("Comparative Profiling", tabName = "comparativeanalysis"),
                                                  menuSubItem("FAQ", tabName ="FAQ")
                                         ),
                                         menuItem("Data Upload and Summary", icon = icon("upload", verify_fa = FALSE), tabName = "Upload", selected = TRUE),
                                         menuItem("Data Processing", icon = icon("filter", verify_fa = FALSE), tabName = "DataProcessing"),
                                         menuItem("Comp. Phenotypic Profiling", icon = icon("dna", verify_fa = FALSE), tabName = "DEAnalysis"),
                                         menuItem("Compositional Profiling", icon = icon("chart-bar", verify_fa = FALSE), tabName = "CellComp"),
                                         menuItem("Comparative Profiling", icon = icon("project-diagram", verify_fa = FALSE), tabName = "Profile"),
                                         menuItem("Summary", icon = icon("project-diagram", verify_fa = FALSE), tabName = "Summary")
                             ),
                             helpText("Developed by ", a("UMMS Biocore.", href="https://www.umassmed.edu/biocore/", target = "_blank"))
                    ),
                    # Discover Auxiliary Tools Sidebar
                    tabPanel(title = "Discover", value = "discover", id="discover",
                             mainPlotControlsUI("deresults"),
                             menuItem("DEFilter",  startExpanded = FALSE,
                                      cutOffSelectionUI("deresults")),
                             menuItem("ScoreFilter", startExpanded = FALSE,
                                      ScoreCutOffSelectionUI("deresults")
                             ),
                             shinydashboard::menuItem("Computational Prof. Heatmaps", 
                                                      heatmapControlsUI("deresults")),
                             shinydashboard::menuItem("Compositional Prof. Heatmaps", 
                                                      heatmapControlsUI("deconvolute")),
                             shinydashboard::menuItem("Comparative Prof. Heatmaps", 
                                                      heatmapControlsUI("profiling"))
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
                   
                   # Data Upload and Summary Tab
                   tabItem(tabName="Upload",
                           MainLoadPage()
                           # tabBox(id = "AnalysisBox",
                           #        width = NULL,
                           #        tabPanel(title = "Import Bulk RNA-Seq",
                           #                 actionButtonDE("art", "art")
                           #        )
                           # ),
                           # tabBox(id = "UploadBox", 
                           #        width = NULL,
                           #        tabPanel(title = "Import Bulk RNA-Seq",
                           #                 dataUI("load")
                           #        ),
                           #        tabPanel(title = "Reference scRNA-Seq",
                           #                 dataSCUI("load"),
                           #                 value = "uploadsummary"
                           #        ),
                           #        tabPanel(title = "Reference Bulk RNA-Seq",
                           #                 dataProfileUI("load"),
                           #                 value = "uploadprofilesummary"
                           #        ),
                           #        tabPanel(title = "Dprofiler Database",
                           #                 DBUI("mongoserver"), 
                           #                 value = "dprofilerdatabase"
                           #        )
                           # )
                   ),

                   # Data Processing Tab
                   tabItem(tabName="DataProcessing", 
                           tabBox(id = "DataProcessingBox", 
                                  width = NULL,
                                  tabPanel(title = "Filtering",
                                           conditionalPanel((condition <- "input.load-Filter"),
                                                            dataLCFUI("lcf")),
                                           value = "filter"
                                  ),
                                  tabPanel(title = "Batch Effect Correction",
                                           conditionalPanel((condition <- "input.Batch"),
                                                            batchEffectUI("batcheffect")),
                                           value = "batcheffect", 
                                  )
                           )
                   ), 
                   
                   # DE Analysis Tab
                   tabItem(tabName="DEAnalysis", 
                           getDEResultsUI("deresults")
                   ),
                   
                   # Deconvolution Tab
                   tabItem(tabName="CellComp", 
                           getDeconvoluteUI("deconvolute")
                   ),
                   
                   # Profiling Tab
                   tabItem(tabName="Profile", 
                           getProfilingUI("profiling")
                   ),
                   
                   # Summary Tab
                   tabItem(tabName="Summary", 
                           SummaryUI("summariser")
                   )
          ),
          
        
        # hide and show tabs as analysis progresses
        # change 'hide' part to manage order of popping up sidebar menus
        getTabUpdateJS() 
      ))
  )
  )
  dprofiler
}