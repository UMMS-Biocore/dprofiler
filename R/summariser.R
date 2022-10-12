#' summaryServer
#' 
#' a module server to summarise all membership scores
#'
#' @param input input
#' @param output output
#' @param session session
#' @param comppheno results from computational phenotypic profiling
#' @param composprof results from cmopositional profiling
#' @param comparprof results from comparative profiling
#' 
#' @examples
#'     x <- summaryServer()
#'     
#' @export
#' 
dprofilersummariser <- function (input = NULL, output = NULL, session = NULL, parent_session = NULL, comppheno = NULL, composprof = NULL, comparprof = NULL) {
  
  ###
  # Reactive Values for Summarizing ####
  ###
  
  SummaryScoreTables <- reactiveValues(comppheno = NULL, comparprof = NULL, composprof = NULL)
  InsertMongo <- reactive(input$recordscores)
  
  # get scores from the comp pheno profiling module
  observeEvent(comppheno$Summarise(), {
    SummaryScoreTables$comppheno <- comppheno$ScoreTable()
    updateTabItems(parent_session, "MenuItems", "Summary")
  })
  
  # get scores from the compositional profiling module
  observeEvent(composprof$Summarise(), {
    SummaryScoreTables$composprof <- composprof$ScoreTable()
    updateTabItems(parent_session, "MenuItems", "Summary")
  })
  
  # get scores from the comparative profiling module
  observeEvent(comparprof$Summarise(), {
    SummaryScoreTables$comparprof <- comparprof$ScoreTable()
    updateTabItems(parent_session, "MenuItems", "Summary")
  })
  
  ScoreTable <- reactive({
    if(all(sapply(reactiveValuesToList(SummaryScoreTables), is.null))) return(NULL)
    
    # parse tables
    summary_table_list <- reactiveValuesToList(SummaryScoreTables)
    summary_table_list <- summary_table_list[c("comppheno","comparprof","composprof")]
    summary_table_list <- summary_table_list[!sapply(summary_table_list, is.null)]
    print(summary_table_list)
    
    # merge tables
    summary_table <- NULL
    for(i in 1:length(summary_table_list)){
      if(i==1){
        summary_table <- summary_table_list[[i]]
      } else {
        summary_table <- as_tibble(summary_table) %>% dplyr::left_join(summary_table_list[[i]])
      }
    }
    # colnames_summary_table <- colnames(summary_table)
    # summary_table <- summary_table[,c(colnames_summary_table[grepl("Score|Sample|Cross",colnames_summary_table)], 
    #                                   colnames_summary_table[!grepl("Score|Sample|Cross",colnames_summary_table)])]
    
    return(summary_table)
  })
  
  ###
  # Main Observable ####
  ###
  
  observe({
    
    # get scores
    getAllScoreDetails(output, session, "MembershipScores", ScoreTable(), modal = FALSE, highlight = TRUE)
    
  })
  
  return(list(ScoreTable = ScoreTable, InsertMongo = InsertMongo))
}

#' SummaryUI
#' 
#' Creates the Mongo Database Pane
#'
#' @param id namespace id
#' 
#' @examples
#'     x <- SummaryUI("summariser")
#'     
#' @export
#' 
SummaryUI <- function (id) {
  ns <- NS(id)
  list(
    tabBox(id = "AllResults",
           width = NULL,
           tabPanel(title = "Summary",
                    fluidRow(
                      shinydashboard::box(title = "Membership Scores",
                                          solidHeader = T, status = "info",  width = 12, collapsible = TRUE,
                                          DT::dataTableOutput(ns("MembershipScores")),
                                          column(2,actionButtonDE(ns("recordscores"), "Record Scores", styleclass = "primary", 
                                                                  style = "width: 100%; margin-top: 25px; margin-left: 0px"))
                      )
                    ),
                    value = "allscores"
           )
    )
  )
}

#' getAllScoreDetails
#' 
#' Details and scores of all profiling methods at the final Summary Pane
#'
#' @param output, output
#' @param session, session
#' @param tablename, table name
#' @param data, matrix data
#' @param modal, if it is true, the matrix is going to be in a modal
#' @param highlight if it is true, numerical columns are highlighted
#'
#' @examples
#'     x <- getAllScoreDetails()
#'     
#' @export
#'     
getAllScoreDetails <- function(output = NULL, session = NULL, tablename = NULL, data = NULL, 
                                   modal = NULL, highlight = FALSE){
  if (is.null(data)) return(NULL)
  output[[tablename]] <- DT::renderDataTable({
    if (!is.null(data)){
      dttable <- DT::datatable(data, extensions = 'Buttons',
                               rownames = FALSE,
                               options = list( server = TRUE,
                                               dom = "Blfrtip",
                                               buttons = 
                                                 list("copy", list(
                                                   extend = "collection"
                                                   , buttons = c("csv", "excel", "pdf")
                                                   , text = "Download"
                                                 ) ), # end of buttons customization
                                               
                                               # customize the length menu
                                               lengthMenu = list( c(10, 20,  50, -1) # declare values
                                                                  , c(10, 20, 50, "All") # declare titles
                                               ), # end of length Menu customization
                                               pageLength = 10))
      numeric_names <- colnames(data[,sapply(as.data.frame(data), is.numeric), drop = FALSE])
      numeric_names <- numeric_names[!numeric_names %in% "Reads"]
      dttable <- dttable %>% DT::formatRound(numeric_names, digits=3)
      if(highlight){
        colours <- rainbow(length(numeric_names))
        for(i in 1:length(numeric_names)){
          dttable <-  dttable %>% DT::formatStyle(numeric_names[i],
                                                  background = DT::styleColorBar(c(0,1), colours[i]))
        }
      } 
      dttable
    }
  })
}