#' dprofilerdeconvolute
#'
#' Module to perform and visualize deconvolution results.
#'
#' @param input input variables
#' @param output output objects
#' @param session session
#' @param dc de results
#' @param scdata single cell data
#' @param sc_marker_table single cell marker table
#' @param parent_session parent session
#' 
#' @return DE panel
#'
#' @examples
#'     x <- dprofilerdeconvolute() 
#'     
#' @export
#' 
dprofilerdeconvolute <- function(input = NULL, output = NULL, session = NULL, parent_session = NULL, data = NULL) {
  
  startDec <- data$startDec
  data <- data$selectedData
  
  ###
  # Reactive Values ####
  ###
  
  DecResults <- reactiveValues(prop = NULL)
  ScoreTable <- reactiveVal()
  Summarise <- reactive(input$summarisescores)
  
  ###
  # Start Deconvolute ####
  ###
  
  observeEvent(startDec(), {
    waiter_show(html = spin_ring(), color = transparent(.5))
    withProgress(message = 'Running RNA Deconvolution', value = 0, {
      mixtures()
    })
    waiter_hide()
    showTab("menutabs","discover", session = parent_session)
    showTab("DeconvoluteBox","deconvoluteresults", session = parent_session)
    updateTabsetPanel(session = parent_session, "DeconvoluteBox", "deconvoluteresults")
  })
  
  # counts, normalization if asked
  counts <- reactive({
    if(grepl("CPM",isolate(input$deconvolute_norm))){
      countdata <- edgeR::cpm(data()$count)
      param <- strsplit(isolate(input$deconvolute_norm), split = "\\+")[[1]][2]
    } else {
      countdata <- data()$count
      param <- isolate(input$deconvolute_norm)
    }
    countdata <- getNormalizedMatrix(countdata, method = param)
    countdata
  })

  # Deconvolution
  mixtures <- reactive({
    
      # select genes
      degenes <- getAllMarkerGenes(data()$markers, data()$scdata, data()$count, input)
      
      # deconvolute
      prop <- deconvolute(counts(), degenes, data()$scdata, input)
      DecResults$prop <- prop 
  })
  
  ScoreTable <- reactive({
    if(is.null(DecResults$prop)) return(NULL)
    # if(!is.null(data()$score)){
    #   if(!is.null(DecResults$prop)){
    #     scoretable <- cbind(data()$score, DecResults$prop)
    #   } 
    # } else {
    #   scoretable <- DecResults$prop
    # }
    scoretable <- DecResults$prop
    data.frame(Sample = rownames(scoretable), scoretable[,!colnames(scoretable) %in% "Sample"], row.names = NULL)
  })
  
  ###
  # Heatmap Reactive Events ####
  ###
  
  # prepare heat data
  data_de_tmm <- reactive({
    marker_genes <- getHeatmapMarkerGenes(data()$markers, data()$count, input)
    if(is.null(marker_genes)) return(NULL)

    # prep data for heatmap
    heatdata <- prepHeatData(counts(), input)
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
             collapsible = TRUE, title = session$ns("Markers"), status = "info",
             solidHeader = TRUE, width = NULL,
             draggable = TRUE, height = 600,
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

  ###
  # Window Settings ####
  ###
  
  # hide initial tabs
  hideTab("DeconvoluteBox","deconvoluteresults", session = parent_session)

  ###
  # Main Observable ####
  ###
  
  # Observe for Tables and Plots
  observe({

    # Score and deconvolution paper
    getDeconvoluteTableDetails(output, session, "MembershipScoresIterDE", ScoreTable(),
                               modal = FALSE, highlight = TRUE)

    # Summary of the scRNA data
    getSCRNAEmbedding(output, "uploadSummaryTSNE", data()$scdata, input$conditions_from_meta0)

  })

  return(list(DecResults = DecResults, Summarise = Summarise, ScoreTable = ScoreTable))
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
#' @export
#' 
getDeconvoluteUI<- function (id) {
  ns <- NS(id)
  list(
    tabBox(id = "DeconvoluteBox",
           width = NULL,
           tabPanel(title = "Conditions",
                    fluidRow(
                      shinydashboard::box(title = "Select Annotations",
                                          solidHeader = T, status = "info",  width = 12, collapsible = TRUE,
                                          uiOutput(ns("conditionSelector")),
                                          column(4,actionButtonDE(ns("deconvolute"), "Start", 
                                                                  styleclass = "primary", style = 'margin-top:21px'))
                                          
                      )
                    ),
                    value = "deconvoluteconditions"
           ),
           tabPanel(title = "Cellular Compositions",
                    fluidRow(
                      # column(2,actionButtonDE("recordcompositionalscores", "Record Scores", styleclass = "primary", 
                      #                         style = "width: 100%; margin-top: 25px; margin-left: 0px"))
                    ),
                    fluidRow(
                      shinydashboard::box(title = "RNA Deconvolution",
                                          solidHeader = T, status = "info",  width = 12, collapsible = TRUE,
                                          #p(strong("Note:"), " P65_NL has a low membership score and ", strong("estimated melanocyte proportion of P65_NL is lower than other non-lesional samples"), 
                                          #  " suggesting that profile of P65_NL might be similar to lesional samples due to ", strong("its melanocytes being low in number.")),
                                          DT::dataTableOutput(ns("MembershipScoresIterDE")),
                                          column(2,actionButtonDE(ns("summarisescores"), "Summarise Scores", styleclass = "primary", 
                                                                  style = "width: 100%; margin-top: 25px; margin-left: 0px"))
                      )
                    ),
                    fluidRow(
                      uiOutput(ns("heatmapUI")),
                      shinydashboard::box(title = "RNA Deconvolution",
                                          solidHeader = T, status = "info",  width = 6, collapsible = TRUE, height = 600, 
                                          fluidRow(
                                            column(12,div(style = 'overflow: scroll', 
                                                          plotOutput(ns("uploadSummaryTSNE")))
                                            ),
                                          )
                      )
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
#' @export
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

#' deconvolute
#' 
#' the deconvolution function based on MuSiC algorithm
#'
#' @param data Bulk expression data set
#' @param DEgenes DE genes for limiting the genes of scRNA and Bulk RNA data sets
#' @param scdata single cell ExpressionSet Object
#' @param input input
#'
#' @examples
#'     x <- deconvolute()
#'     
#' @export
#'  
deconvolute <- function(data = NULL, DEgenes = NULLL, scdata = NULL, input = NULL){
  if (is.null(data)) return(NULL)
  
  # parameters
  celltypes <- isolate(input$condition)
  celltype_label <- isolate(input$conditions_from_meta0)
  samples <- isolate(input$deconvolute_samples)
  method <- isolate(input$deconvolute_methods)
  
  # Bulk Data
  genes_data <- rownames(data)
  DEgenes <- intersect(genes_data, DEgenes)
  data_de <- data[DEgenes,]
  Vit_BulkRNAseq <- ExpressionSet(assayData=as.matrix(data_de))
  pData(Vit_BulkRNAseq) <- data.frame(row.names = colnames(data), columns =  colnames(data))
  
  # Single cell data
  metadata <- pData(scdata)
  scdata <- scdata[,metadata[,celltype_label] %in% celltypes] 
  metadata <- metadata[metadata[,celltype_label] %in% celltypes,]
  scdata <- scdata[rownames(scdata) %in% DEgenes,]
  
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
  if(method == "MuSIC"){
    NLandL.prop = music_prop(bulk.eset = Vit_BulkRNAseq, 
                             sc.eset = scdata, 
                             clusters = celltype_label,
                             samples = samples, 
                             verbose = T)
    res <- NLandL.prop$Est.prop.weighted
  } else if(method =="BisqueRNA"){
    res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset = Vit_BulkRNAseq, 
                                                  sc.eset = scdata, 
                                                  cell.types = celltype_label,
                                                  subject.names = samples,
                                                  use.overlap = FALSE, 
                                                  markers = DEgenes)
    res <- t(res$bulk.props)
  } else if(method =="SCDC"){
    res <- SCDC_prop(bulk.eset = Vit_BulkRNAseq, 
                     sc.eset = scdata, 
                     ct.varname = celltype_label,
                     sample = samples,
                     ct.sub = celltypes)
    res <- res$prop.est.mvw
  }
  
  return(res)
}


#' getHeatmapMarkerGenes
#'
#' @param sc_marker_table single cell marker table
#' @param data reference bulk data
#' @param input input 
#'
#' @examples
#'     x <- getHeatmapMarkerGenes()
#'     
#' @export
#'     
getHeatmapMarkerGenes <- function(sc_marker_table = NULL, data = NULL, input = NULL){
  if (is.null(sc_marker_table) || is.null(input$select_celltype)) return(NULL)
  
  # pull those genes that are in bulk data
  marker_table_genes <- unique(sc_marker_table$gene)
  common_genes <- intersect(marker_table_genes, rownames(data))
  sc_marker_table <- sc_marker_table[sc_marker_table$gene %in% common_genes,]
  
  # take out additional genes
  sc_marker_table <- sc_marker_table[!grepl("^AC[0-9]|^AL[0-9]|^AL[0-9]|^MT-", sc_marker_table$gene),]
  
  # grep cell type specific markers
  sc_marker_table <- sc_marker_table[sc_marker_table$cluster %in% input$select_celltype,]
  top_n_markers <- as.numeric(input$select_top_markers)
  top_n_markers <- ifelse(is.na(top_n_markers), nrow(sc_marker_table),
                          ifelse(top_n_markers > nrow(sc_marker_table), nrow(sc_marker_table), top_n_markers))
  marker_genes <- sc_marker_table$gene[order(sc_marker_table$avg_log2FC, decreasing = TRUE)[1:top_n_markers]]
  
  return(marker_genes)
}

#' getAllMarkerGenes
#'
#' @param sc_marker_table single cell marker table
#' @param scdata single cell data object
#' @param data bulk data
#' @param columns columns of the bulk data
#' @param input input 
#'
#' @examples
#'     x <- getAllMarkerGenes()
#'     
#' @export
#'  
getAllMarkerGenes <- function(sc_marker_table = NULL, scdata = NULL, data = NULL, input = NULL){
  if (is.null(data)) return(NULL)
  # data <- data[,columns]
  
  # if all genes are requested, return the intersecting genes
  if(isolate(input$allgenes) == "Yes"){
    gene_list <- intersect(rownames(data),rownames(scdata))
    return(gene_list)
  }
  
  # select cell types, and other conditions
  sc_marker_table <- sc_marker_table[sc_marker_table$Level %in% isolate(input$conditions_from_meta0),]
  sc_marker_table <- sc_marker_table[sc_marker_table$cluster %in% isolate(input$condition),]
  sc_marker_table <- sc_marker_table[sc_marker_table$pct.1 > as.numeric(isolate(input$pct1)) & 
                                       sc_marker_table$pct.2 < as.numeric(isolate(input$pct2)) & 
                                       sc_marker_table$avg_log2FC > as.numeric(isolate(input$logFC)) & 
                                       sc_marker_table$p_val_adj < as.numeric(isolate(input$padj)),]
  
  # pull those genes that are in bulk data
  marker_table_genes <- unique(sc_marker_table$gene)
  common_genes <- intersect(marker_table_genes, rownames(data))
  sc_marker_table <- sc_marker_table[sc_marker_table$gene %in% common_genes,]
  
  # delete duplicate genes
  num_genes <- table(sc_marker_table$gene)
  duplicates <- names(num_genes[num_genes>1])
  sc_marker_table <- sc_marker_table[!sc_marker_table$gene %in% duplicates,]
  
  # take out additional genes
  sc_marker_table <- sc_marker_table[!grepl("^AC[0-9]|^AL[0-9]|^MT-", sc_marker_table$gene),]
  
  # Find mean gene abundance in bulk data
  initial_markers <- unique(sc_marker_table$gene)
  data_de <- data[rownames(data) %in% initial_markers,]
  data_de <- as.data.frame(data_de)
  data_de$gene <- rownames(data_de)
  data_de$celltype_markers <- sapply(data_de$gene, function(x) sc_marker_table$cluster[which(x == initial_markers)[1]])
  datax_de_mean_count <- data.frame(Count = rowMeans(data_de[,1:(ncol(data_de)-2)]), gene = data_de$gene, celltype = data_de$celltype_markers)
  
  # eliminate outlier gene
  datax_de_mean_count <- datax_de_mean_count %>%
    group_by(celltype) %>% 
    mutate(Outlier = remove_outliers(Count))
  outlier_genes <- datax_de_mean_count$gene[datax_de_mean_count$Outlier]
  sc_marker_table <- sc_marker_table[!sc_marker_table$gene %in% outlier_genes,]
  
  # pick top genes
  sc_marker_table %>%
    group_by(cluster) %>%
    top_n(n = as.numeric(isolate(input$top_genes)), wt = avg_log2FC) -> marker_table_topgenes 
  gene_list <- unique(marker_table_topgenes$gene)
  
  return(gene_list)
}

#' getSCRNAEmbedding
#' 
#' get single cell RNA Embedding
#'
#' @param output output
#' @param summary summary output name
#' @param data single cell ExpressionSet Object
#' @param ident column in scRNA metadata to visualize 
#'
#' @examples
#'       x <- getSCRNAEmbedding()
#'     
#' @export
#'    
getSCRNAEmbedding <- function (output = NULL, summary = NULL, data = NULL, ident = NULL) 
{
  if (is.null(output)) 
    return(NULL)
  
  if(is.null(ident))
    return(NULL)
  
  if(!is.null(data)){
    
    if(any(!ident %in% colnames(pData(data))))
      return(NULL)
    
    tsne_data <- pData(data)[,c(ident,"x","y")]
    tsne_data <- as_tibble(tsne_data) %>% group_by_at(1) %>% summarize(x = mean(x), y = mean(y))
    output[[summary]] <- renderPlot({
      SignallingSingleCell::plot_tsne_metadata(data, color_by = ident,
                                               legend_dot_size = 4, text_sizes = c(20, 10, 5, 8, 8, 15)) + 
        geom_text(tsne_data, mapping = aes(x = x, y = y, label = tsne_data[[1]]), size=5)
    })
  }
}

remove_outliers <- function(x, na.rm = TRUE, probs=c(.25, .75), ...) {
  qnt <- quantile(x, probs=probs, na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- ifelse((x < (qnt[1] - H)) | (x > (qnt[2] + H)), TRUE, FALSE)
  y
}

###
# old functions ####
###

#' prepDeconvolute
#' 
#' Prepares the container for deconvolution methods
#'
#' @param dc reactive object of DE analysis
#' @param scdata single cell ExpressionSet Object
#' @param sc_marker_table single cell markers table
#' @param parent_session parent session
#'
#' @examples
#'        x <- prepDeconvolute()
#'     
#' @export
#' 
prepDeconvolute_old <- function(dc = NULL, scdata = NULL, sc_marker_table = NULL, parent_session = NULL){
  if (is.null(dc)) return(NULL)
  waiter_show(html = spin_ring(), color = transparent(.5))
  withProgress(message = 'Running RNA Deconvolution', value = 0, {
    mixtures <- callModule(dprofilerdeconvolute, "deconvolute", dc, scdata, sc_marker_table, parent_session)
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
#' @param sc_marker_table single cell marker table
#' @param parent_session parent session
#' 
#' @return DE panel
#'
#' @examples
#'     x <- dprofilerdeconvolute() 
#'     
#' @export
#' 
dprofilerdeconvolute_old <- function(input = NULL, output = NULL, session = NULL, dc = NULL, 
                                 scdata = NULL, sc_marker_table = NULL, parent_session = NULL) {
  if(is.null(dc)) return(NULL)
  
  # get columns
  if(!is.null(dc()$cols)) columns <- dc()$cols
  else columns <- colnames(dc()$count)
  
  # counts, normalization if asked
  counts <- reactive({
    if(grepl("CPM",isolate(input$deconvolute_norm))){
      results <- edgeR::cpm(dc()$count[,columns])
      param <- strsplit(isolate(input$deconvolute_norm), split = "\\+")[[1]][2]
    } else {
      results <- dc()$count[,columns]
      param <- isolate(input$deconvolute_norm)
    }
    results <- getNormalizedMatrix(results, method = param)
    results
  })
  
  # Deconvolution
  mixtures <- reactive({
    # select genes
    degenes <- getAllMarkerGenes(sc_marker_table, scdata, dc()$count, columns, input)
    
    # deconvolute
    mixture <- deconvolute(counts(), degenes, scdata, input)
    updateTabsetPanel(session = parent_session, "DeconvoluteBox", "deconvoluteresults")
    mixture
  })
  
  # prepare heat data
  data_de_tmm <- reactive({
    marker_genes <- getHeatmapMarkerGenes(sc_marker_table, dc()$count, input)
    if(is.null(marker_genes)) return(NULL)
    
    # prep data for heatmap
    heatdata <- prepHeatData(counts(), input)
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
             collapsible = TRUE, title = session$ns("Markers"), status = "info", 
             solidHeader = TRUE, width = NULL,
             draggable = TRUE, height = 600,
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
    
    # If the scores are not available, should be NULL
    if(!is.null(dc()$CrossScore)){
      ScoreTable <- MergeScoreTables(dc()$CrossScore()$IterDEscore, dc()$IntraScore()$Score)
    } else {
      ScoreTable <- NULL
    } 
    
    # Score and deconvolution paper
    Reads <- data.frame(Reads = colSums(dc()$count[,columns]))
    ScoreTable <- cbind(ScoreTable, mixtures())
    ScoreTable <- cbind(Reads, ScoreTable)
    getDeconvoluteTableDetails(output, session, "MembershipScoresIterDE", ScoreTable, 
                               modal = FALSE, highlight = TRUE)
    
    # Summary of the scRNA data
    getSCRNAEmbedding(output, "uploadSummaryTSNE", scdata, input$conditions_from_meta0)
    
  })
  
  return(mixtures = mixtures)
}
