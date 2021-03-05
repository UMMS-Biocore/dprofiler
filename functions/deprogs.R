#' dprofilerdeanalysis
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
#'     x <- dprofilerdeanalysis()
#'
dprofilerdeanalysis <- function(input = NULL, output = NULL, session = NULL, 
                                data = NULL, columns = NULL, conds = NULL, params = NULL) {
    if(is.null(data)) return(NULL)
    
    # Iterative DE Algorithm
    deres <- reactive({
        runIterDE(data, columns, conds, params)
    })
    
    # Apply Filters for DE and Iter DE Results
    prepDat <- reactive({
        applyFiltersNew(addDataCols(data, deres()$DEResults, columns, conds), input)
    })
    iterprepDat <- reactive({
        applyFiltersNew(addDataCols(data, deres()$IterDEResults, columns, conds), input)
    })
    
    # # get Final Scores
    # finalScore <- reactive({
    #     getFinalScores(deres(), data, columns, conds, params, 
    #                    ManuelDEgenes = input$manuelgenes, TopStat = input$topstat)
    # })
    
    # Observe for Tables and Plots
    observe({
        finalScore <- getFinalScores(deres(), data, columns, conds, params, 
                                     ManuelDEgenes = input$manuelgenes, TopStat = input$topstat)
        dat <-  prepDat()[prepDat()$Legend == input$legendradio,]
        dat2 <- removeCols(c("ID", "x", "y","Legend", "Size"), dat)
        iterdat <-  iterprepDat()[iterprepDat()$Legend == input$legendradio,]
        iterdat2 <- removeCols(c("ID", "x", "y","Legend", "Size"), iterdat)
        
        # DE Results
        getTableDetails(output, session, "DEResults", dat2, modal=FALSE)
        getTableDetails(output, session, "IterDEResults", iterdat2, modal = FALSE)

        # Membership Scores
        getIterDESummary(output, session, "HomogeneityVenn", "HomogeneitySummary", deres())
        getScoreDetails(output, session, "HomogeneityScores", finalScore$DEscore, finalScore$IterDEscore)
        # getScoreTableDetails(output, session, "MembershipScoresIterDE", finalScore$IterDEscore, 
        #                      modal = FALSE, highlight = TRUE)
        
    })
    list(dat = prepDat)
}


#' getDEResultsUI
#' Creates a panel to visualize DE results
#'
#' @param id, namespace id
#' @return panel
#' @examples
#'     x <- getDEResultsUI("batcheffect")
#'
#' @export
#'
getDEResultsUI<- function (id) {
    ns <- NS(id)
    list(
        tabBox(id = "DEAnalysisBox",
               width = NULL,
               tabPanel(title = "Homogeneity Detection",
                        fluidRow(
                            column(12,
                                shinydashboard::box(title = "Summary",
                                                solidHeader = T, status = "info",  width = 3, collapsible = TRUE,
                                                tableOutput(ns("HomogeneitySummary"))
                                )
                            ),
                            column(12,
                                shinydashboard::box(title = "Venn diagram",
                                                    solidHeader = T, status = "info",  width = 6, collapsible = TRUE,
                                                    plotOutput(ns("HomogeneityVenn")),
                                                    actionButtonDE("deconvolute", "Go to Cellular Composition", styleclass = "primary")
                                ),
                                shinydashboard::box(title = "Membership Scores",
                                                    solidHeader = T, status = "info",  width = 6, collapsible = TRUE,
                                                    plotlyOutput(ns("HomogeneityScores"))
                                ),
                                # shinydashboard::box(title = "Cellular Heterogeneity Analysis",
                                #                     solidHeader = T, status = "info",  width = 12, collapsible = TRUE,
                                #                     # uiOutput(ns("MembershipScoresIterDE"))
                                #                     DT::dataTableOutput(ns("MembershipScoresIterDE"))
                                # )
                            )
                        ),
                        value = "homogeneity"
               ),
               tabPanel(title = "Heterogeneous Conditions",
                        fluidRow(
                            shinydashboard::box(title = "Differentially Expressed Genes",
                                                solidHeader = T, status = "info",  width = 12, collapsible = TRUE,
                                                fluidRow(
                                                    uiOutput(ns("IterDEResults"))
                                                )
                            )
                        ),
                        value = "iterderesults"
               ),
               tabPanel(title = "Homogeneous Conditions",
                        fluidRow(
                            shinydashboard::box(title = "Differentially Expressed Genes",
                                                solidHeader = T, status = "info",  width = 12, collapsible = TRUE,
                                                fluidRow(
                                                    uiOutput(ns("DEResults")),
                                                    actionButtonDE("goMain", "Go to Main Plots", styleclass = "primary")
                                                )
                            )
                        ),
                        value = "deresults"
               )
        )
    )
}

#' cutOffSelectionUI
#'
#' Gathers the cut off selection for DE analysis
#'
#' @param id, namespace id
#' @note \code{cutOffSelectionUI}
#' @return returns the left menu according to the selected tab;
#' @examples
#'     x <- cutOffSelectionUI("cutoff")
#' @export
#'
cutOffSelectionUI <- function(id){
    ns <- NS(id)
    list(
        getLegendRadio(id),
        textInput(ns("padj"), "padj value cut off", value = "0.01" ),
        textInput(ns("foldChange"), "foldChange", value = "2" ),
        textInput(ns("topstat"), "Top Stat", value = "" ), 
        fileInput(ns("manuelgenes"), "Manuel DEgenes")
    )
}