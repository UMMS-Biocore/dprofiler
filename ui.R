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
        uiOutput("loading"),
        tabsetPanel(id = "menutabs", type = "tabs",
                    tabPanel(title = "Data Prep", value = "dataprep", id="dataprep",
                             sidebarMenu(id="MenuItems",
                                         menuItem("Upload", icon = icon("upload"), tabName = "Upload"),
                                         menuItem("Data Processing", icon = icon("filter"), tabName = "DataProcessing"),
                                         # menuItem("Filter", icon = icon("filter"), tabName = "Filter"),
                                         # menuItem("BatchEffect",  icon = icon("align-left"), tabName = "BatchEffect"),
                                         menuItem("Cond. Select", icon = icon("adjust"), tabName = "CondSelect"),
                                         menuItem("DE Analysis", icon = icon("adjust"), tabName = "DEAnalysis"),
                                         # menuItem("Iter. DEAnalysis", icon = icon("adjust"), tabName = "IterDEAnalysis"),
                                         menuItem("DEFilter",  icon = icon("code"), tabName = "CondSelect",  startExpanded = TRUE,
                                                  uiOutput("cutOffUI"),
                                                  uiOutput("compselectUI"))
                             ),
                             helpText("Developed by ", a("UMMS Biocore.", href="https://www.umassmed.edu/biocore/", target = "_blank")))
        )
      ),
      
      # Shiny dashboard Body
      dashboardBody(
          
          tabItems(#id = "methodtabs", type = "tabs",

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