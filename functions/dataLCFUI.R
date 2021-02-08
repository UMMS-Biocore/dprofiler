#' dataLCFUI
#'
#' @param id 
#'
#' @return
#' @export
#'
#' @examples
dataLCFUI<- function (id) {
  ns <- NS(id)
  list(
    fluidRow(
      shinydashboard::box(title = "Low Count Filtering",
                          solidHeader = TRUE, status = "info",  width = 12, 
                          fluidRow(
                            column(5,div(style = 'overflow: scroll',
                                         tableOutput(ns("uploadSummary")),
                                         DT::dataTableOutput(ns("sampleDetails"))),
                                   uiOutput(ns("loadedtable"))
                            ),
                            column(2,
                                   shinydashboard::box(title = "Filtering Methods",
                                                       solidHeader = TRUE, status = "info",
                                                       width = 12, 
                                                       lcfMetRadio(id),
                                                       uiOutput(ns("cutoffLCFMet")),
                                                       actionButtonDE(ns("submitLCF"), label = "Filter", styleclass = "primary")
                                   )
                            ),
                            column(5,div(style = 'overflow: scroll',
                                         tableOutput(ns("filteredSummary")),
                                         DT::dataTableOutput(ns("filteredDetails"))),
                                   uiOutput(ns("filteredtable"))
                            )
                          ),
                          conditionalPanel(condition = paste0("input['", ns("submitLCF"),"']"),
                                           actionButtonDE("Batch", label = "Batch Effect Correction", styleclass = "primary"),
                                           conditionalPanel(condition = "!(input.Batch)",
                                                            actionButtonDE("goDEFromFilter", "Iterative DE Analysis", styleclass = "primary")))
      ),
      shinydashboard::box(title = "Histograms",
                          solidHeader = TRUE, status = "info",  width = 12, 
                          fluidRow(
                            column(6,histogramControlsUI(ns("beforeFiltering")),
                                   getHistogramUI(ns("beforeFiltering"))),
                            column(6,histogramControlsUI(ns("afterFiltering")),
                                   getHistogramUI(ns("afterFiltering")))
                          ))
    ))
}