#' dprofilerdeconvolute
#'
#' Module to perform and visualize DE results.
#'
#' @param input, input variables
#' @param output, output objects
#' @param session, session
#' @param data, a matrix that includes expression values
#' @param columns, columns
#' @param conds, conditions
#' @param params, de parameters
#' @return DE panel
#' @export
#'
#' @examples
#'     x <- dprofilerdeconvolute()
#'
dprofilerdeconvolute <- function(input = NULL, output = NULL, session = NULL, scores = NULL) {
  if(is.null(scores)) return(NULL)

  # Observe for Tables and Plots
  observe({
    getScoreTableDetails(output, session, "MembershipScoresIterDE", scores$DEscore, modal = FALSE,
                         highlight = TRUE)
  })
  list(dat = scores)
}

deconvolute <- function(data, deres, columns, conds, scdata){
  
  data <- data[,columns]
  data_de <- data[deres$IterDEgenes,]
  Vit_BulkRNAseq <- ExpressionSet(assayData=as.matrix(data_de))
  pData(Vit_BulkRNAseq) <- data.frame(row.names = columns,
                                      conds = conds)
  
  NLandL.prop = music_prop(bulk.eset = Vit_BulkRNAseq, 
                           sc.eset = scdata, 
                           clusters = 'CellType',
                           samples = 'Patient', verbose = T)
  
  return(NLandL.prop$Est.prop.weighted)
}


#' getDeconvoluteUI
#' Creates a panel to visualize DE results
#'
#' @param id, namespace id
#' @return panel
#' @examples
#'     x <- getDeconvoluteUI("batcheffect")
#'
#' @export
#'
getDeconvoluteUI<- function (id) {
  ns <- NS(id)
  list(
    tabBox(id = "DeconvoluteBox",
           width = NULL,
           tabPanel(title = "Cellular Composition",
                    fluidRow(
                      column(12,
                             shinydashboard::box(title = "Cellular Heterogeneity Analysis",
                                                 solidHeader = T, status = "info",  width = 12, collapsible = TRUE,
                                                 DT::dataTableOutput(ns("MembershipScoresIterDE"))
                             )
                      )
                    ),
                    value = "homogeneity"
           )
    )
  )
}

#' getScoreTableDetails
#' 
#' get table details
#' To be able to put a table into two lines are necessary;
#' into the server part;
#' getScoreTableDetails(output, session, "dataname", data, modal=TRUE)
#' into the ui part;
#' uiOutput(ns("dataname"))
#'   
#' @param output, output
#' @param session, session
#' @param tablename, table name
#' @param data, matrix data
#' @param modal, if it is true, the matrix is going to be in a modal
#' @return panel
#' @examples
#'     x <- getScoreTableDetails()
#'
#' @export
#'
getScoreTableDetails <- function(output  = NULL, session  = NULL, tablename  = NULL, data = NULL, 
                                 modal = NULL, highlight = FALSE){
  if (is.null(data)) return(NULL)
  output[[tablename]] <- DT::renderDataTable({
    if (!is.null(data)){
      dttable <- DT::datatable(data, extensions = 'Buttons',
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
                                               ), # end of lengthMenu customization
                                               pageLength = 10))
      numeric_names <- colnames(data[,sapply(data, is.numeric), drop = FALSE])
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