#' prepDeconvolute
#' 
#' Prepares the container for deconvolution methods
#'
#' @param dc reactive object of DE analysis
#' @param scdata single cell ExpressionSet Object
#'
#' @return
#' @export
#'
#' @examples
prepDeconvolute <- function(dc = NULL, scdata = NULL){
  if (is.null(dc)) return(NULL)
  
  withProgress(message = 'Running RNA Deconvolution', value = 0, {
    mixtures <- callModule(dprofilerdeconvolute, "deconvolute", dc, scdata)
    mix <- mixtures()
  })
  
  return(mix)
}


#' dprofilerdeconvolute
#'
#' Module to perform and visualize deconvolution results.
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
dprofilerdeconvolute <- function(input = NULL, output = NULL, session = NULL, dc = NULL, scdata = NULL) {
  if(is.null(dc)) return(NULL)

  # Deconvolution
  mixtures <- reactive({
      deconvolute(dc()$init_dedata, dc()$IterDEgenes, dc()$cols, scdata)
  })
  
  # Choose Cell Types and top markers
  output$heatmap_selection <- renderUI({
    list(
      column(3,
             selectInput(session$ns("select_celltype"), label = "Select Celltype", choices = unique(pData(scdata)$CellType))),
      column(3,
            textInput(session$ns("select_top_markers"), label = "Top n Markers", value = "10"))
    )
  })
  
  # heatmap of marker genes
  output$heatmap <- renderPlot({
    featuresData <- fData(scdata)[rownames(fData(scdata)) %in% dc()$IterDEgenes,]
    gene_scores <- featuresData[,paste0(input$select_celltype,"_marker_score_CellType")]
    top_n_markers <- as.numeric(input$select_top_markers)
    top_n_markers <- ifelse(is.na(top_n_markers), length(gene_scores),
                            ifelse(top_n_markers > length(gene_scores), length(gene_scores), top_n_markers))
    marker_genes <- rownames(featuresData)[order(gene_scores, decreasing = TRUE)[1:top_n_markers]]
    data_de_tmm <- getNormalizedMatrix(dc()$init_dedata[marker_genes,dc()$cols], method = "TMM")
    runHeatmap2(input, session, expdata = as.matrix(data_de_tmm))
  })
  
  # Observe for Tables and Plots
  observe({
    
    # Score and deconvolution paper 
    ScoreTable <- cbind(dc()$score()$IterDEscore, mixtures())
    getScoreTableDetails(output, session, "MembershipScoresIterDE", ScoreTable, 
                         modal = FALSE, highlight = TRUE)
    
  })
  
  return(mixtures = mixtures)
}


#' deconvolute
#' 
#' the deconvolution function based on MuSiC algorithm
#'
#' @param data Bulk expression data set
#' @param DEgenes DE genes for limiting the genes of scRNA and Bulk RNA data sets
#' @param columns samples that are deconvoluted
#' @param scdata single cell ExpressionSet Object
#'
#' @return
#' @export
#'
#' @examples
deconvolute <- function(data, DEgenes, columns, scdata){
  
  data <- data[,columns]
  data_de <- data[DEgenes,]
  Vit_BulkRNAseq <- ExpressionSet(assayData=as.matrix(data_de))
  pData(Vit_BulkRNAseq) <- data.frame(row.names = columns, columns = columns)

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
                      shinydashboard::box(title = "Cellular Heterogeneity Analysis",
                                          solidHeader = T, status = "info",  width = 12, collapsible = TRUE,
                                          DT::dataTableOutput(ns("MembershipScoresIterDE"))
                      ),
                      shinydashboard::box(
                        collapsible = TRUE, title = "Markers", status = "primary", 
                        solidHeader = TRUE, width = 6,
                        draggable = TRUE, 
                        column(12,
                               plotOutput(ns("heatmap"))
                        ),
                        column(8,
                               uiOutput(ns("heatmap_selection"))
                        )
                      )
                    )
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
#' @param highlight if it is true, numerical columns are highlighted
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
  
  