#' prepDeconvolute
#' 
#' Prepares the container for deconvolution methods
#'
#' @param dc reactive object of DE analysis
#' @param scdata single cell ExpressionSet Object
#' @param parent_session parent session
#'
#' @examples
#'        x <- prepDeconvolute()
#'        
prepDeconvolute <- function(dc = NULL, scdata = NULL, parent_session = NULL){
  if (is.null(dc)) return(NULL)
  waiter_show(html = spin_ring(), color = transparent(.5))
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
#' @param input input variables
#' @param output output objects
#' @param session session
#' @param dc de results
#' @param scdata single cell data
#' @param parent_session parent session
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
      if(isolate(input$deconvolute_genes) == "DE Genes After Prof."){
        degenes <- dc()$IterDEgenes
      } else if(isolate(input$deconvolute_genes) == "DE Genes Before Prof."){
        degenes <- dc()$DEgenes
      } else {
        degenes <- getAllMarkerGenes(scdata, dc()$init_dedata, input)
      }
      mixture <- deconvolute(dc()$init_dedata, degenes, dc()$cols, scdata, input)
      updateTabsetPanel(session = parent_session, "DeconvoluteBox", "deconvoluteresults")
      mixture
  })
  
  # prepare heat data
  data_de_tmm <- reactive({
    marker_genes <- getMarkerGenes(scdata, dc()$IterDEgenes, input)
    if(is.null(marker_genes)) return(NULL)
    # heatdata <- getNormalizedMatrix(dc()$init_dedata[,dc()$cols], method = "TMM")
    # heatdata <- prepHeatData(heatdata, input)
    heatdata <- prepHeatData(dc()$init_dedata[,dc()$cols], input)
    as.matrix(heatdata[marker_genes,])
  })
  
  output$heatmap <- renderPlotly({
    if(!is.null(data_de_tmm())){
      if(nrow(data_de_tmm()) > 1 & ncol(data_de_tmm()) > 1)
        withProgress(message = 'Drawing Heatmap', detail = "interactive", value = 0, {
          runHeatmap(input, session, data_de_tmm())
        })
    }
  })
  
  output$heatmap2 <- renderPlot({
    if(!is.null(data_de_tmm())){
      if(nrow(data_de_tmm()) > 1 & ncol(data_de_tmm()) > 1)
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
             selectInput(session$ns("select_celltype"), label = "Select Celltype", choices = isolate(input$condition))),
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
#' 
#' Creates a panel to visualize deconvolution results
#'
#' @param id, namespace id
#' @return panel
#' @examples
#'     x <- getDeconvoluteUI("batcheffect")
#'
getDeconvoluteUI<- function (id) {
  ns <- NS(id)
  list(
    tabBox(id = "DeconvoluteBox",
           width = NULL,
           tabPanel(title = "Mixture Conditions",
                    fluidRow(
                      shinydashboard::box(title = "Select Annotations",
                                          solidHeader = T, status = "info",  width = 12, collapsible = TRUE,
                                          uiOutput(ns("conditionSelector")),
                                          column(4,actionButtonDE("deconvolute", "Start", 
                                                                  styleclass = "primary", style = 'margin-top:21px'))
                                          
                      )
                    ),
                    value = "deconvoluteconditions"
           ),
           tabPanel(title = "Cellular Compositions",
                    fluidRow(
                      shinydashboard::box(title = "RNA Deconvolution",
                                          solidHeader = T, status = "info",  width = 12, collapsible = TRUE,
                                          p(strong("Note:"), " P65_NL has a low membership score and ", strong("estimated melanocyte proportion of P65_NL is lower than other non-lesional samples"), 
                                            " suggesting that profile of P65_NL might be similar to lesional samples due to ", strong("its melanocytes being low in number.")),
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
#' get deconvolution table details
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
getDeconvoluteTableDetails <- function(output  = NULL, session  = NULL, tablename  = NULL, data = NULL, 
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
#' @param input input
#'
#' @examples
#'     x <- deconvolute()
#'     
deconvolute <- function(data = NULL, DEgenes = NULL, columns = NULL, scdata = NULL, input = NULL){
  if (is.null(data)) return(NULL)
  
  # parameters
  celltypes <- isolate(input$condition)
  celltype_label <- isolate(input$conditions_from_meta0)
  samples <- isolate(input$deconvolute_samples)
  
  # Bulk Data
  data <- data[,columns]
  # data <- getNormalizedMatrix(data,method = "TMM")
  genes_data <- rownames(data)
  DEgenes <- intersect(genes_data, DEgenes)
  data_de <- data[DEgenes,]
  Vit_BulkRNAseq <- ExpressionSet(assayData=as.matrix(data_de))
  pData(Vit_BulkRNAseq) <- data.frame(row.names = columns, columns = columns)
  
  # Single cell data
  metadata <- pData(scdata)
  scdata <- scdata[,metadata[,celltype_label] %in% celltypes] 
  metadata <- metadata[metadata[,celltype_label] %in% celltypes,]
    
  # normalized integrated library sizes
  # exprs_srt <- exprs(scdata)
  # metadata$nCount_integratedRNA_norm <- colSums(exprs_srt)
  # metadata <- as_tibble(metadata) %>% group_by(CellType) %>% mutate(nCount_integratedRNA_normmean = mean(nCount_integratedRNA_norm))
  # exprs_srt <- apply(exprs_srt, 1,function(x){
  #   return((x/metadata$nCount_integratedRNA_normmean)*100)
  # })
  # scdata <- ExpressionSet(assayData=t(exprs_srt))
  # pData(scdata) <- data.frame(metadata)
  # rownames(scdata) <- colnames(exprs_srt)
  # colnames(scdata) <- rownames(exprs_srt)

  # deconvolute
  NLandL.prop = music_prop(bulk.eset = Vit_BulkRNAseq, 
                           sc.eset = scdata, 
                           clusters = celltype_label,
                           samples = samples, verbose = T)
  
  return(NLandL.prop$Est.prop.weighted)
}


#' getMarkerGenes
#'
#' @param scdata single cell data
#' @param IterDEgenes DE genes from bulk data
#' @param input input 
#'
#' @examples
#'     x <- getMarkerGenes()
#'     
getMarkerGenes <- function(scdata = NULL, IterDEgenes = NULL, input = NULL){
  
  if (is.null(scdata) | is.null(input$select_celltype)) 
    return(NULL)
  
  if(!is.null(IterDEgenes)){
    featuresData <- fData(scdata)[rownames(fData(scdata)) %in% IterDEgenes,]
  } else {
    featuresData <- fData(scdata)
  }
  gene_scores <- featuresData[,paste0(input$select_celltype,"_marker_score_CellType")]
  top_n_markers <- as.numeric(input$select_top_markers)
  top_n_markers <- ifelse(is.na(top_n_markers), length(gene_scores),
                          ifelse(top_n_markers > length(gene_scores), length(gene_scores), top_n_markers))
  marker_genes <- rownames(featuresData)[order(gene_scores, decreasing = FALSE)[1:top_n_markers]]
 
  return(marker_genes)
}

#' getAllMarkerGenes
#'
#' @param scdata single cell data
#' @param data reference bulk data
#' @param input input 
#'
#' @examples
#'     x <- getAllMarkerGenes()
#'     
getAllMarkerGenes <- function(scdata = NULL, data = NULL, input = NULL){
  if (is.null(scdata)) return(NULL)
  
  # pull those genes that are in bulk data
  featuresData <- fData(scdata)
  rownames_fdata <- rownames(featuresData)
  rownames_fdata <- intersect(rownames_fdata, rownames(data))
  featuresData <- featuresData[rownames_fdata, ]
  
  # select cell types and their ranks
  gene_scores <- featuresData[,paste0(isolate(input$condition),"_marker_score_CellType"), drop = FALSE]
  top_n_markers <- as.numeric(isolate(input$top_genes))
  top_n_markers <- ifelse(is.na(top_n_markers), nrow(gene_scores),
                          ifelse(top_n_markers > nrow(gene_scores), nrow(gene_scores), top_n_markers))
  
  # pull genes
  gene_list <- c()
  marker_genes <- list()
  for(i in 1:ncol(gene_scores)){
    marker_genes[[i]] <- rownames(featuresData)[order(gene_scores[,i], decreasing = FALSE)[1:top_n_markers]]
  }
  all_markers <- unlist(marker_genes)
  duplicates <- names(table(all_markers)[table(all_markers) > 1])
  gene_list <- unique(all_markers)
  gene_list <- gene_list[!gene_list %in% duplicates]
  return(gene_list)
}
  