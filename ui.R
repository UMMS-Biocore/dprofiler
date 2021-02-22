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

dconnectUI <- function() {
  dbHeader <- shinydashboard::dashboardHeader(titleWidth = 250)
  dbHeader$children[[2]]$children <- tags$a(style='color: white;',
                                            id="top_logo" , "Dconnect")
  addResourcePath(prefix = "www", directoryPath = "www/")
  library("debrowser")
  dconnect <- (fluidPage(
    shinyjs::useShinyjs(),
    tags$head(tags$title("Dconnect"),
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
                                         menuItem("Filter", icon = icon("filter"), tabName = "Filter"),
                                         menuItem("BatchEffect",  icon = icon("align-left"), tabName = "BatchEffect"),
                                         menuItem("Cond. Select", icon = icon("adjust"), tabName = "CondSelect"),
                                         menuItem("DEAnalysis", icon = icon("adjust"), tabName = "DEAnalysis"),
                                         menuItem("Iter. DEAnalysis", icon = icon("adjust"), tabName = "IterDEAnalysis"),
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

                   # Filter Tab
                   tabItem(tabName="Filter",
                           conditionalPanel((condition <- "input.Filter"),
                                            dataLCFUI("lcf"))
                   ),

                   # Batch Effect Tab
                   tabItem(tabName="BatchEffect",
                           conditionalPanel((condition <- "input.Batch"),
                                            batchEffectUI("batcheffect"))
                   ), 
                   
                   # Condition Selection Tab
                   tabItem(tabName="CondSelect", 
                           condSelectUI()
                   ),
                   
                   # Upload Tab
                   tabItem(tabName="DEAnalysis", 
                           tabBox(id = "DEAnalysisBox", 
                                  width = NULL,
                                  tabPanel(title = "DE Analysis",
                                           uiOutput("deresUI"), 
                                           value = "deresults"
                                  )
                           )
                   ),
                   
                   # Upload Tab
                   tabItem(tabName="IterDEAnalysis", 
                           tabBox(id = "IterDEAnalysisBox", 
                                  width = NULL,
                                  tabPanel(title = "Iterative DE Analysis",
                                           uiOutput("iterderesUI"), 
                                           value = "membershipscore"
                                  )
                           )
                   )
          )
          
        
        # hide and show tabs as analysis progresses
        # change 'hide' part to manage order of popping up sidebar menus
        # getTabUpdateJS() 
      ))
  )
  )
  dconnect
}