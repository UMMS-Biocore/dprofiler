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
                                data = NULL, columns = NULL, conds = NULL, params = NULL,
                                scdata = NULL) {
    if(is.null(data)) return(NULL)
    
    # Iterative DE Algorithm
    deres <- reactive({
        runIterDE(data, columns, conds, params)
    })
    
    # # Deconvolution
    # mixtures <- reactive({
    #     deconvolute(data, deres(), columns, conds, scdata)
    # })
    
    # Apply Filters for DE and Iter DE Results
    prepDat <- reactive({
        applyFiltersNew(addDataCols(data, deres()$DEResults, columns, conds), input)
    })
    iterprepDat <- reactive({
        remaining_columns <- (columns != deres()$cleaned_columns)
        if(length(remaining_columns) > 0){
            columns <- columns[remaining_columns]
            conds <- conds[remaining_columns]
        }
        applyFiltersNew(addDataCols(data, deres()$IterDEResults, columns, conds), input)
    })
    
    # Observe for Tables and Plots
    observe({
        finalScore <- getFinalScores(deres(), data, columns, conds, params, 
                                     ManuelDEgenes = input$manualgenes, TopStat = input$topstat)
        dat <-  prepDat()[prepDat()$Legend == input$legendradio,]
        dat2 <- removeCols(c("ID", "x", "y","Legend", "Size"), dat)
        iterdat <-  iterprepDat()[iterprepDat()$Legend == input$legendradio,]
        iterdat2 <- removeCols(c("ID", "x", "y","Legend", "Size"), iterdat)
        
        # DE Results
        getTableDetails(output, session, "DEResults", dat2, modal=FALSE)
        getTableDetails(output, session, "IterDEResults", iterdat2, modal = FALSE)

        # Membership Scores
        getIterDESummary(output, session, "HomogeneityVenn", "HomogeneitySummary", deres(), params)
        getScoreDetails(output, session, "HomogeneityScores", finalScore$DEscore, finalScore$IterDEscore)
        getScoreTableDetails(output, session, "MembershipScoresIterDE", 
                             # cbind(finalScore$IterDEscore, mixtures()), 
                             cbind(finalScore$IterDEscore), 
                             modal = FALSE, highlight = TRUE)
        
        # download handlers
        output$downloadBeforeGenes <- downloadHandler(
          filename = function() {paste('initial_degenes.txt')},
          content = function(con) {write(setdiff(deres()$DEgenes,deres()$IterDEgenes), con)}
        )
        output$downloadOverlapGenes <- downloadHandler(
            filename = function() {paste('overlapping_degenes.txt')},
            content = function(con) {write(intersect(deres()$DEgenes,deres()$IterDEgenes), con)}
        )
        output$downloadAfterGenes <- downloadHandler(
            filename = function() { paste('final_degenes.txt')},
            content = function(con) {write(setdiff(deres()$IterDEgenes,deres()$DEgenes), con)}
        )
    })
    list(dat = prepDat, DEgenes = deres()$DEgenes, iterdat = iterprepDat, IterDEgenes = deres()$IterDEgenes)
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
                            shinydashboard::box(title = "Summary", height = 260,
                                                solidHeader = T, status = "info",  width = 12, collapsible = TRUE,
                                                column(3,tableOutput(ns("HomogeneitySummary"))),
                                                column(3,tableOutput(ns("HomogeneitySummaryIter"))),
                                                column(3,tableOutput(ns("HomogeneitySummaryDE")))
                            ),
                            shinydashboard::box(title = "# of DE Genes",
                                                solidHeader = T, status = "info",  width = 6, collapsible = TRUE,
                                                plotOutput(ns("HomogeneityVenn")),
                                                column(12,
                                                       downloadButton(ns("downloadBeforeGenes"), label = "Initial DEgenes"),
                                                       downloadButton(ns("downloadOverlapGenes"), label = "Overlapping DE genes"),
                                                       downloadButton(ns("downloadAfterGenes"), label = "Final DE genes"),
                                                )
                            ),
                            shinydashboard::box(title = "Membership Scores",
                                                solidHeader = T, status = "info",  width = 6, collapsible = TRUE,
                                                plotlyOutput(ns("HomogeneityScores")),
                                                actionButtonDE("deconvolute", "Go to Cellular Composition", styleclass = "primary")
                            ),
                            uiOutput(ns("maininitialplot")),                                  
                            uiOutput(ns("mainoverlapplot")),
                            uiOutput(ns("mainfinalplot")),
                            uiOutput(ns("BarMainUI")),
                            uiOutput(ns("BoxMainUI")),
                            uiOutput(ns("heatmapUI"))
                        ),
                        value = "homogeneity"
               ),
               tabPanel(title = "Impure (Heterogeneous) Conditions",
                        fluidRow(
                            shinydashboard::box(title = "Differentially Expressed Genes",
                                                solidHeader = T, status = "info",  width = 12, collapsible = TRUE,
                                                uiOutput(ns("IterDEResults"))
                                                
                            ),
                            uiOutput(ns("maindeplot"))
                        ),
                        value = "deresults"
               ),
               tabPanel(title = "Pure (Homogeneous) Conditions",
                        fluidRow(
                            shinydashboard::box(title = "Differentially Expressed Genes",
                                                solidHeader = T, status = "info",  width = 12, collapsible = TRUE,
                                                uiOutput(ns("DEResults"))
                            ),
                            uiOutput(ns("mainiterdeplot"))
                        ),
                        value = "iterderesults"
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
        fileInput(ns("manualgenes"), "Manual DEgenes")
    )
}