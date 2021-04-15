#' prepDeconvolute
#' 
#' Prepares the container for deconvolution methods
#'
#' @param dc reactive object of DE analysis
#' @param scdata single cell ExpressionSet Object
#' @param parent_session main session
#'
#' @return
#' @export
#'
#' @examples
prepDeconvolute <- function(dc = NULL, scdata = NULL, parent_session = NULL){
  if (is.null(dc)) return(NULL)
  
  waiter_show(html = waiting_screen, color = transparent(.5))
  withProgress(message = 'Running RNA Deconvolution', value = 0, {
    mixtures <- callModule(dprofilerdeconvolute, "deconvolute", dc, scdata, parent_session)
    mix <- mixtures()
  })
  waiter_hide()
  
  return(mix)
}


#' dprofilerdeconvolute
#'
#' Module to perform and visualize deconvolution results.
#'
#' @param input, input variables
#' @param output, output objects
#' @param session, session
#' @param dc, de results
#' @param scdata, single cell data
#' @param parent_session main session
#' 
#' @return DE panel
#' @export
#'
#' @examples
#'     x <- dprofilerdeconvolute()
#'
dprofilerdeconvolute <- function(input = NULL, output = NULL, session = NULL, dc = NULL, 
                                 scdata = NULL, parent_session = NULL) {
  if(is.null(dc)) return(NULL)

  # Deconvolution
  mixtures <- reactive({
      mixture <- deconvolute(dc()$init_dedata, isolate(dc()$deconvolute_genes()), dc()$cols, scdata)
      updateTabsetPanel(session = parent_session, "DeconvoluteBox", "deconvoluteresults")
      mixture
  })
  
  # prepare heat data
  data_de_tmm <- reactive({
    marker_genes <- getMarkerGenes(scdata, dc()$IterDEgenes, input)
    heatdata <- prepHeatData(dc()$init_dedata[,dc()$cols], input)
    as.matrix(heatdata[marker_genes,])
  })
  
  output$heatmap <- renderPlotly({
    if(!is.null(data_de_tmm())){
      withProgress(message = 'Drawing Heatmap', detail = "interactive", value = 0, {
        runHeatmap(input, session, data_de_tmm())
      })
    }
  })
  
  output$heatmap2 <- renderPlot({
    if(!is.null(data_de_tmm())){
      withProgress(message = 'Drawing Heatmap', detail = "non-interactive", value = 0, {
        runHeatmap2(input, session, data_de_tmm())
      })
    }
  })
  
  output$heatmapUI <- renderUI({
    if (is.null(input$interactive)) return(NULL)
    column(6,
           shinydashboard::box(
             collapsible = TRUE, title = session$ns("Markers"), status = "primary", 
             solidHeader = TRUE, width = NULL,
             draggable = TRUE,
             column(12,getPlotArea(input, session)),
             column(8,uiOutput(session$ns("heatmap_selection")))
           )
    )
  })
  
  # Choose Cell Types and top markers on heatmap
  output$heatmap_selection <- renderUI({
    list(
      column(6,
             selectInput(session$ns("select_celltype"), label = "Select Celltype", choices = unique(pData(scdata)$CellType))),
      column(6,
             textInput(session$ns("select_top_markers"), label = "Top n Markers", value = "10"))
    )
  })
  
  # Observe for Tables and Plots
  observe({
    
    # Score and deconvolution paper 
    ScoreTable <- cbind(dc()$score()$IterDEscore, mixtures())
    getDeconvoluteTableDetails(output, session, "MembershipScoresIterDE", ScoreTable, 
                         modal = FALSE, highlight = TRUE)
    
  })
  
  return(mixtures = mixtures)
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
           tabPanel(title = "Conditions",
                    fluidRow(
                      shinydashboard::box(title = "Select Conditions",
                                          solidHeader = T, status = "info",  width = 12, collapsible = TRUE,
                                          uiOutput(ns("conditionSelector")),
                                          column(4,actionButtonDE("deconvolute", "Start", 
                                                                  styleclass = "primary", style = 'margin-top:21px'))
                                          
                      )
                    ),
                    value = "deconvoluteconditions"
           ),
           tabPanel(title = "Cellular Composition",
                    fluidRow(
                      shinydashboard::box(title = "Cellular Heterogeneity Analysis",
                                          solidHeader = T, status = "info",  width = 12, collapsible = TRUE,
                                          DT::dataTableOutput(ns("MembershipScoresIterDE"))
                      ),
                      uiOutput(ns("heatmapUI"))
                    ),
                    value = "deconvoluteresults"
           )
    )
  )
}

#' getDeconvoluteTableDetails
#' 
#' get table details
#' To be able to put a table into two lines are necessary;
#' into the server part;
#' getDeconvoluteTableDetails(output, session, "dataname", data, modal=TRUE)
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
#'     x <- getDeconvoluteTableDetails()
#'
#' @export
#'
getDeconvoluteTableDetails <- function(output  = NULL, session  = NULL, tablename  = NULL, data = NULL, 
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
  
  data <- getNormalizedMatrix(data,method = "TMM")
  
  data_de <- data[DEgenes,]
  Vit_BulkRNAseq <- ExpressionSet(assayData=as.matrix(data_de))
  pData(Vit_BulkRNAseq) <- data.frame(row.names = columns, columns = columns)
  
  
  NLandL.prop = music_prop(bulk.eset = Vit_BulkRNAseq, 
                           sc.eset = scdata, 
                           clusters = 'CellType',
                           samples = 'Patient', verbose = T)
  
  return(NLandL.prop$Est.prop.weighted)
}


#' getMarkerGenes
#'
#' @param scdata single cell data
#' @param IterDEgenes DE genes from bulk data
#' @param input input 
#'
#' @examples
getMarkerGenes <- function(scdata, IterDEgenes, input){
  
  if(!is.null(IterDEgenes)){
    featuresData <- fData(scdata)[rownames(fData(scdata)) %in% IterDEgenes,]
  } else {
    featuresData <- fData(scdata)
  }
  gene_scores <- featuresData[,paste0(input$select_celltype,"_marker_score_CellType")]
  top_n_markers <- as.numeric(input$select_top_markers)
  top_n_markers <- ifelse(is.na(top_n_markers), length(gene_scores),
                          ifelse(top_n_markers > length(gene_scores), length(gene_scores), top_n_markers))
  marker_genes <- rownames(featuresData)[order(gene_scores, decreasing = TRUE)[1:top_n_markers]]
 
  return(marker_genes)
}
  