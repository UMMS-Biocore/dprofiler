#' dconnectServer
#'
#' Sets up shinyServer to be able to run DEBrowser interactively.
#'
#' @note \code{dconnectServer}
#' @param input, input params from UI
#' @param output, output params to UI
#' @param session, session variable
#' @return the panel for main plots;
#'
#' @export
#' @importFrom shiny actionButton actionLink addResourcePath column 
#'             conditionalPanel downloadButton downloadHandler 
#'             eventReactive fileInput fluidPage helpText isolate 
#'             mainPanel need numericInput observe observeEvent 
#'             outputOptions parseQueryString plotOutput radioButtons 
#'             reactive reactiveValues renderPlot renderUI runApp 
#'             selectInput shinyApp  shinyServer  shinyUI sidebarLayout 
#'             sidebarPanel sliderInput  stopApp  tabPanel tabsetPanel 
#'             textInput textOutput titlePanel uiOutput tags HTML
#'             h4 img icon updateTabsetPanel updateTextInput  validate 
#'             wellPanel checkboxInput br p checkboxGroupInput onRestore
#'             reactiveValuesToList renderText onBookmark onBookmarked 
#'             updateQueryString callModule enableBookmarking htmlOutput
#'             onRestored NS reactiveVal withProgress tableOutput
#'             selectizeInput fluidRow div renderPrint renderImage
#'             verbatimTextOutput imageOutput renderTable incProgress
#'             a h3 strong h2 withMathJax updateCheckboxInput
#'             showNotification updateSelectInput
#' @importFrom shinyjs show hide enable disable useShinyjs extendShinyjs
#'             js inlineCSS onclick
#' @importFrom DT datatable dataTableOutput renderDataTable formatStyle
#'             styleInterval formatRound
#' @importFrom ggplot2 aes aes_string geom_bar geom_point ggplot
#'             labs scale_x_discrete scale_y_discrete ylab
#'             autoplot theme_minimal theme geom_density
#'             geom_text element_blank margin
#' @importFrom plotly renderPlotly plotlyOutput plot_ly add_bars event_data
#'             hide_legend %>% group_by ggplotly
#' @importFrom gplots heatmap.2 redblue bluered
#' @importFrom igraph layout.kamada.kawai  
#' @importFrom grDevices dev.off pdf colorRampPalette 
#' @importFrom graphics barplot hist pairs par rect text plot
#' @importFrom stats aggregate as.dist cor cor.test dist
#'             hclust kmeans na.omit prcomp var sd model.matrix
#'             p.adjust runif cov mahalanobis quantile as.dendrogram
#'             density
#' @importFrom utils read.csv read.table write.table update.packages
#'             download.file read.delim data install.packages
#'             packageDescription installed.packages
#' @importFrom DOSE enrichDO
#' @importFrom enrichplot gseaplot dotplot
#' @importMethodsFrom DOSE summary
#' @importMethodsFrom AnnotationDbi as.data.frame as.list colnames
#'             exists sample subset head mappedkeys ncol nrow subset 
#'             keys mapIds select
#' @importMethodsFrom GenomicRanges as.factor setdiff
#' @importMethodsFrom IRanges as.matrix "colnames<-" mean
#'             nchar paste rownames toupper unique which
#'             as.matrix lapply "rownames<-" gsub
#' @importMethodsFrom S4Vectors eval grep grepl levels sapply t 
#' @importMethodsFrom SummarizedExperiment cbind order rbind
#' @importFrom jsonlite fromJSON
#' @importFrom methods new
#' @importFrom stringi stri_rand_strings
#' @importFrom annotate geneSymbols
#' @importFrom reshape2 melt
#' @importFrom Harman harman reconstructData
#' @importFrom clusterProfiler compareCluster enrichKEGG enrichGO gseGO bitr
#' @importFrom DESeq2 DESeq DESeqDataSetFromMatrix results estimateSizeFactors
#'             counts lfcShrink
#' @importFrom edgeR calcNormFactors equalizeLibSizes DGEList glmLRT
#'             exactTest estimateCommonDisp glmFit topTags
#' @importFrom shinydashboard dashboardHeader dropdownMenu messageItem
#'             dashboardPage dashboardSidebar sidebarMenu dashboardBody
#'             updateTabItems menuItem tabItems tabItem menuSubItem
#' @importFrom limma lmFit voom eBayes topTable
#' @importFrom sva ComBat
#' @importFrom RCurl getURL
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import shinyBS
#' @import colourpicker
#' @import RColorBrewer
#' @import heatmaply
#' @import apeglm
#' @import ashr

dconnectServer <- function(input, output, session) {
  options(warn = -1)
  tryCatch(
    {
      if (!interactive()) {
        options( shiny.maxRequestSize = 30 * 1024 ^ 2,
                 shiny.fullstacktrace = FALSE, shiny.trace=FALSE, 
                 shiny.autoreload=TRUE, warn =-1)
      }
      
      choicecounter <- reactiveValues(nc = 0)
      
      output$programtitle <- renderUI({
        # togglePanels(0, c(0), session)
        getProgramTitle(session)
      })
      
      uploadeddata <- reactiveVal()
      filtd <- reactiveVal()
      batch <- reactiveVal()
      sel <- reactiveVal()
      dc <- reactiveVal()
      compsel <- reactive({
        cp <- 1
        if (!is.null(input$compselect_dataprep))
          cp <- input$compselect_dataprep
        cp
      })
      
      observe({
        
        # Data Loading Event
        uploadeddata(callModule(dataLoadServer, "load", "Filter", session))
        updateTabItems(session, "MenuItems", "Upload")
        
        # Data Filtering Event
        observeEvent (input$Filter, {
          if(!is.null(uploadeddata()$load())){
            updateTabItems(session, "MenuItems", "Filter")
            filtd(callModule(debrowserlowcountfilter, "lcf", uploadeddata()$load()))
          }
        })
        
        # Batch Correction Event
        observeEvent (input$Batch, {
          if(!is.null(filtd()$filter())){
            updateTabItems(session, "MenuItems", "BatchEffect")
            batch(callModule(debrowserbatcheffect, "batcheffect", filtd()$filter()))
          }
        })
        
        
        # Regular DE analysis Event
        observeEvent (input$goDE, {
          if(is.null(batch())) batch(setBatch(filtd()))
          updateTabItems(session, "MenuItems", "CondSelect")
          sel(debrowsercondselect(input, output, session,
                                  batch()$BatchEffect()$count, batch()$BatchEffect()$meta))
          choicecounter$nc <- sel()$cc()
        })
        observeEvent (input$goDEFromFilter, {
          if(is.null(batch())) batch(setBatch(filtd()))
          updateTabItems(session, "MenuItems", "CondSelect")
          sel(debrowsercondselect(input, output, session,
                                  batch()$BatchEffect()$count, batch()$BatchEffect()$meta))
          choicecounter$nc <- sel()$cc()
        })
        observeEvent (input$startDE, {
          if(!is.null(batch()$BatchEffect()$count)){
            togglePanels(0, c(0), session)
            res <- prepDataContainer(batch()$BatchEffect()$count, sel()$cc(), input)
            if(is.null(res$de)) return(NULL)
            dc(res$de)
            updateTabItems(session, "MenuItems", "DEAnalysis")
            buttonValues$startDE <- TRUE
          }
        })

        output$compselectUI <- renderUI({
          if (!is.null(sel()) && !is.null(sel()$cc()))
            getCompSelection("compselect_dataprep",sel()$cc())
        })

        output$cutOffUI <- renderUI({
          cutOffSelectionUI(paste0("DEResults", compsel()))
        })
        output$deresUI <- renderUI({
          getDEResultsUI(paste0("DEResults",compsel()))
        })
        output$iterderesUI <- renderUI({
          getIterDEResultsUI(paste0("DEResults",compsel()))
        })
      })
      
      # check if conditions are ready for DE Analysis 
      output$condReady <- reactive({
        if (!is.null(sel()))
          choicecounter$nc <- sel()$cc()
        choicecounter$nc
      })
      
      # auxiliary values of buttons and reactive values
      buttonValues <- reactiveValues(goQCplots = FALSE, goDE = FALSE,
                                     startDE = FALSE)

    },
    err=function(errorCondition) {
      cat("in err handler")
      message(errorCondition)
    },
    warn=function(warningCondition) {
      cat("in warn handler")
      message(warningCondition)
    })
}