#' dprofilerdeanalysis
#'
#' Module to perform and visualize Iterative DE results.
#' 
#' @param input, input variables
#' @param output, output objects
#' @param session, session 
#' @param data, a matrix that includes expression values
#' @param columns, columns
#' @param conds, conditions
#' @param params, de parameters
#' 
#' @return DE panel 
#' @export
#'
#' @examples
#'     x <- dprofilerdeanalysis()
#'
dprofilerdeanalysis <- function(input = NULL, output = NULL, session = NULL, 
                                data = NULL, columns = NULL, conds = NULL, params = NULL){
    if(is.null(data)) return(NULL)
    
    # Iterative DE Algorithm
    deres <- reactive({
        runIterDE(data, columns, conds, params, session)
    })
    
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
    
    # Create a reactive scoring table as well
    Scores <- reactiveVal()
    
    # Observe for Tables and Plots
    observe({
        
        # get scpres
        Scores(getFinalScores(deres(), data, columns, conds, params, 
                              ManualDEgenes = input$manualgenes, TopStat = input$topstat))
        
        # prepare DE tables
        dat <-  prepDat()[prepDat()$Legend == input$legendradio,]
        dat2 <- removeCols(c("ID", "x", "y","Legend", "Size"), dat)
        iterdat <-  iterprepDat()[iterprepDat()$Legend == input$legendradio,]
        iterdat2 <- removeCols(c("ID", "x", "y","Legend", "Size"), iterdat)
        
        # DE Results
        getTableDetails(output, session, "DEResults", dat2, modal=FALSE)
        getTableDetails(output, session, "IterDEResults", iterdat2, modal = FALSE)

        # Membership Scores
        getIterDESummary(output, session, "HomogeneityVenn", "HomogeneitySummary", deres(), params)
        getScoreDetails(output, session, "HomogeneityScores", Scores()$DEscore, Scores()$IterDEscore)
        
        # download handler for DE genes
        getDEgenesDownloadButtons(output, session, deres()$DEgenes, deres()$IterDEgenes)
    })
    list(dat = prepDat, DEgenes = deres()$DEgenes, 
         iterdat = iterprepDat, IterDEgenes = deres()$IterDEgenes,
         score = Scores)
}

#' getDEResultsUI
#' 
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
               tabPanel(title = "Differential Hetergeneity Detection",
                        fluidRow(
                            shinydashboard::box(title = "Summary", height = 260,
                                                solidHeader = T, status = "info",  width = 12, collapsible = TRUE,
                                                column(2,htmlOutput(ns("HomogeneitySummary"))),
                                                column(2,htmlOutput(ns("HomogeneitySummaryIter"))),
                                                column(2,htmlOutput(ns("HomogeneitySummaryDE")))
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
                        )               
                        ),
               tabPanel(title = "Pure (Homogeneous) Conditions",
                        fluidRow(
                            shinydashboard::box(title = "Differentially Expressed Genes",
                                                solidHeader = T, status = "info",  width = 12, collapsible = TRUE,
                                                uiOutput(ns("DEResults"))
                            ),
                            uiOutput(ns("mainiterdeplot"))
                        )
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
        textInput(ns("foldChange"), "foldChange", value = "2" )
    )
}

#' ScoreCutOffSelectionUI
#'
#' Gathers the cut off selection for Scoring 
#'
#' @param id, namespace id
#' @note \code{ScoreCutOffSelectionUI}
#' @return returns the left menu according to the selected tab;
#' @examples
#'     x <- ScoreCutOffSelectionUI("cutoff")
#' @export
#'
ScoreCutOffSelectionUI <- function(id){
    ns <- NS(id)
    list(
        textInput(ns("topstat"), "Top Stat", value = "" ),
        fileInput(ns("manualgenes"), "Manual DEgenes")
    )
}

#' getDEgenesDownloadButtons
#'
#' Buttons for downloading Initial, overlapping and DE genes
#' 
#' @param output output 
#' @param session session
#' @param DEgenes Initial DE genes
#' @param IterDEgenes Final DE genes
#'
#' @return
#' @export
#'
#' @examples
getDEgenesDownloadButtons <- function(output, session,  DEgenes, IterDEgenes){
    
    genes <- setdiff(DEgenes,IterDEgenes)
    if(length(genes) == 0) genes <- DEgenes
    output$downloadBeforeGenes <- downloadHandler(
        filename = function() {paste('initial_degenes.txt')},
        content = function(con) {write(genes, con)}
    )
    
    genes <- intersect(DEgenes,IterDEgenes)
    if(length(genes) == 0) genes <- IterDEgenes
    output$downloadOverlapGenes <- downloadHandler(
        filename = function() {paste('overlapping_degenes.txt')},
        content = function(con) {write(genes, con)}
    )
    
    genes <- setdiff(IterDEgenes,DEgenes)
    if(length(genes) == 0) genes <- IterDEgenes
    output$downloadAfterGenes <- downloadHandler(
        filename = function() { paste('final_degenes.txt')},
        content = function(con) {write(genes, con)}
    )
}