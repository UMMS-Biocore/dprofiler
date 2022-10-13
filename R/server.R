#' dprofilerServer
#'
#' Sets up shinyServer to be able to run DEBrowser interactively.
#'
#' @note \code{dprofilerServer}
#' @param input, input params from UI
#' @param output, output params to UI
#' @param session, session variable
#' @return the panel for main plots;
#'
#' @export
#' @importFrom debrowser actionButtonDE applyFiltersNew batchMethod changeClusterOrder 
#'             checkMetaData debrowserdataload fileTypes getAfterLoadMsg
#'             getHistogramUI getIQRPlotUI getLegendRadio getLevelOrder getLogo
#'             getMean getMostVariedList getNormalizedMatrix getPCAPlotUI
#'             getPCAcontolUpdatesJS getSampleDetails getSampleNames
#'             getSelectInputBox getStartupMsg getTableDetails
#'             get_conditions_given_selection heatmapControlsUI histogramControlsUI
#'             lcfMetRadio niceKmeans normalizationMethods pcaPlotControlsUI
#'             plotData plotSizeMarginsUI prepGroup prepHeatData removeCols
#'             runDESeq2 runEdgeR runHeatmap runHeatmap2 selectedInput setBatch
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
#'             showNotification updateSelectInput tagList setProgress
#' @importFrom ggplot2 aes aes_string geom_bar geom_point ggplot
#'             labs scale_x_discrete scale_y_discrete ylab
#'             autoplot theme_minimal theme geom_density
#'             geom_text element_blank margin scale_y_continuous xlab facet_grid
#'             element_text
#' @importFrom shinyjs show hide enable disable useShinyjs extendShinyjs
#'             js inlineCSS onclick disabled
#' @importFrom plotly renderPlotly plotlyOutput plot_ly add_bars event_data
#'             hide_legend %>% group_by ggplotly config
#' @importFrom grDevices dev.off pdf colorRampPalette rainbow
#' @importFrom graphics barplot hist pairs par rect text plot 
#' @importFrom stats aggregate as.dist cor cor.test dist
#'             hclust kmeans na.omit prcomp var sd model.matrix
#'             p.adjust runif cov mahalanobis quantile as.dendrogram
#'             density reorder
#' @importFrom grid grid.newpage grid.draw
#' @importFrom utils read.csv read.table write.table update.packages
#'             download.file read.delim data install.packages
#'             packageDescription installed.packages combn
#' @importFrom jsonlite fromJSON toJSON
#' @importFrom reshape2 melt
#' @importFrom shinydashboard dashboardHeader dropdownMenu messageItem
#'             dashboardPage dashboardSidebar sidebarMenu dashboardBody
#'             updateTabItems menuItem tabItems tabItem menuSubItem tabBox
#' @importFrom SignallingSingleCell plot_tsne_metadata plot_density_ridge
#' @importFrom cluster silhouette 
#' @importFrom Biobase ExpressionSet pData exprs fData
#' @importFrom rlang .data
#' @importFrom MuSiC music_prop 
#' @importFrom VennDiagram draw.pairwise.venn
#' @importFrom waiter spin_ring transparent use_waiter waiter_hide waiter_show
#' @importFrom dplyr summarize group_by_at as_tibble group_by mutate top_n left_join
#' @importFrom shinybusy add_busy_spinner
#' @importFrom limma lmFit voom eBayes topTable
#' @importFrom DT datatable dataTableOutput renderDataTable formatStyle
#'             styleInterval formatRound
#' @importFrom gplots heatmap.2 redblue bluered
#' @importFrom sva ComBat
#' @importFrom RCurl getURL
#' @importFrom DESeq2 DESeq DESeqDataSetFromMatrix results estimateSizeFactors
#'             counts lfcShrink
#' @importFrom Harman harman reconstructData
#' @importFrom edgeR calcNormFactors equalizeLibSizes DGEList glmLRT
#'             exactTest estimateCommonDisp glmFit topTags
#' @importFrom shinyFiles shinyDirButton shinyDirChoose parseDirPath
#' @importMethodsFrom IRanges as.matrix "colnames<-" mean
#'             nchar paste rownames toupper unique which
#'             as.matrix lapply "rownames<-" gsub
#' @importMethodsFrom S4Vectors eval grep grepl levels sapply t 
#' @importMethodsFrom SummarizedExperiment cbind order rbind
#' @importFrom methods new
#' @import apeglm
#' @import ashr
#' @import BisqueRNA
#' @import SCDC
#' @import Seurat
#' @import irlba
#' @import flexclust
#' @import mongolite
#' 
dprofilerServer <- function(input = NULL, output = NULL, session = NULL) {
  options(warn = -1)
  tryCatch(
    {
      if (!interactive()) {
        options( shiny.maxRequestSize = 30 * 1024 ^ 2,
                 shiny.fullstacktrace = FALSE, shiny.trace=FALSE, 
                 shiny.autoreload=TRUE, warn =-1)
      }

      ####
      ## Data Upload ####
      ####
      
      uploadeddata <- callModule(dataLoadServer, "load", session, mongo)
      
      ####
      ## Data Filtering Events ####
      ####
      
      filtd <- callModule(dprofilerFilter, "lcf", session, uploadeddata)
      
      ####
      ## Batch Correction Events ####
      ####
      
      batch <- callModule(dprofilerBatch, "batcheffect", session, filtd)
      
      ####
      ## DE Analysis Events ####
      ####
      
      selectedData <- list(filtd = filtd, batch = batch)
      cond_dc <- dprofilercondselect(input, output, session, session, data = selectedData)
      dc <- callModule(dprofilerdeanalysis, "deresults", session, cond_dc)
      
      # aggregate data for deconvolution and profiling
      selectedData <- list(load = uploadeddata$load, filtd = filtd, batch = batch, cond_dc = cond_dc, dc = dc) 
      
      ####
      ## DE Analysis Main Plot Events ####
      ####
      
      ### Main plots UI ####
      selectedMain <- callModule(dprofilermainplot, "deresults", dc)
      
      ### Heatmaps UI #### 
      selectedHeat <- reactiveVal()
      observe({
        if (!is.null(selectedMain) && !is.null(selectedMain$selGenes())) {
          withProgress(message = 'Creating plot', style = "notification", value = 0.1, {
            selectedHeat(callModule(dprofilerheatmap, "deresults", dc$count()[selectedMain$selGenes(), cond_dc$DEparameters()$cols]))
          })
        }
      })
      
      selgenename <- reactiveVal()
      observe({
        if (!is.null(selectedMain) && !is.null(selectedMain$shgClicked())
            && selectedMain$shgClicked()!=""){
          selgenename(selectedMain$shgClicked())
          if (!is.null(selectedHeat()) && !is.null(selectedHeat()$shgClicked()) &&
              selectedHeat()$shgClicked() != ""){
            js$resetInputParam("deresults-hoveredgenenameclick")
          }
        }
      })
      observe({
        if (!is.null(selectedHeat()) && !is.null(selectedHeat()$shgClicked()) && 
            selectedHeat()$shgClicked() != ""){
          selgenename(selectedHeat()$shgClicked())
        }
      })
      
      ### Bar Main and BoxMain Plots UI ####
      observe({
        if (!is.null(selgenename()) && selgenename()!=""){
          withProgress(message = 'Creating Bar/Box plots', style = "notification", value = 0.1, {
            callModule(dprofilerbarmainplot, "deresults", dc$count(),
                       cond_dc$DEparameters()$cols, cond_dc$DEparameters()$conds, selgenename())
            callModule(dprofilerboxmainplot, "deresults", dc$count(),
                       cond_dc$DEparameters()$cols, cond_dc$DEparameters()$conds, selgenename())
          })
        }
      })
      
      ####
      ## Compositional Profiling Events ####
      ####
      
      # deconvolute after iterative DE analysis
      cond_dec <- callModule(dprofilercondselect, "deconvolute", session, data = selectedData)
      dec <- callModule(dprofilerdeconvolute, "deconvolute", session, data = cond_dec)
      
      ####
      ## Comparative Profiling Events ####
      ####
      
      # profiling with bulk reference 
      cond_prof <- callModule(dprofilercondselect, "profiling", session, data = selectedData)
      prof <- callModule(dprofilerprofiling, "profiling", session, data = cond_prof)
      
      ####
      ## Summarizer Events ####
      ####
      
      summariser <- callModule(dprofilersummariser, "summariser", session, dc, dec, prof)
      
      ####
      ## Mongo Server Events ####
      ####
      
      mongo <- callModule(mongoServer, "mongoserver", session, uploadeddata, summariser)
      
      ####
      ## Auxiliary Shiny commands ####
      ####
      
      # hide initial tabs
      hideTab("DataProcessingBox","batcheffect")
      hideTab("menutabs","discover")
      
      # count data and metadata files for help section
      output$metaFile <-  renderTable({
        read.delim(system.file("extdata", "www", "metaFile.txt",
                               package = "dprofiler"), header=TRUE, sep = " ", skipNul = TRUE)
      })
      output$countFile <-  renderTable({
        read.delim(system.file("extdata", "www", "countFile.txt",
                               package = "dprofiler"), header=TRUE, sep = " ", skipNul = TRUE)
      })

      output$logo <- renderUI({
        getLogo()
      })
      output$startup <- renderUI({
        getStartupMsg()
      })
      output$afterload <- renderUI({
        getAfterLoadMsg()
      })

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