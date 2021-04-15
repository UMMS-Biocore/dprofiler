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

dprofilerServer <- function(input, output, session) {
  options(warn = -1)
  tryCatch(
    {
      if (!interactive()) {
        options( shiny.maxRequestSize = 30 * 1024 ^ 2,
                 shiny.fullstacktrace = FALSE, shiny.trace=FALSE, 
                 shiny.autoreload=TRUE, warn =-1)
      }
      
      ## Reactive Values #### 
      BigImage <- reactiveVal(FALSE)
      uploadeddata <- reactiveVal()
      filtd <- reactiveVal()
      batch <- reactiveVal()
      sel <- reactiveVal()
      sel_prof <- reactiveVal()
      dc <- reactiveVal()
      dec <- reactiveVal()
      
      observe({
        
        ## Start from the Upload Tab ####
        updateTabItems(session, "MenuItems", "Upload")
        
        ## Data Loading Event ####
        uploadeddata(callModule(dataLoadServer, "load", "Filter", session))

        ## Data Filtering Event ####
        observeEvent (input$Filter, {
          if(!is.null(uploadeddata()$load())){
            shinyjs::toggle(id = "dataprocessing")
            updateTabItems(session, "MenuItems", "DataProcessing")
            filtd(callModule(debrowserlowcountfilter, "lcf", uploadeddata()$load()))
          }
        })
        
        ## Batch Correction Event ####
        observeEvent (input$Batch, {
          if(!is.null(filtd()$filter())){
            updateTabsetPanel(session, "DataProcessingBox","batcheffect")
            batch(callModule(debrowserbatcheffect, "batcheffect", filtd()$filter()))
          }
        })
        
        ## Heterogeneity Analysis Events ####
        observeEvent (input$goDE, {
          if(is.null(batch())) batch(setBatch(filtd()))
          sel(dprofilercondselect(input, output, session,
                                  batch()$BatchEffect()$count, batch()$BatchEffect()$meta))
          updateTabItems(session, "MenuItems", "DEAnalysis")
        })
        
        observeEvent (input$goDEFromFilter, {
          if(is.null(batch())) batch(setBatch(filtd()))
          sel(dprofilercondselect(input, output, session,
                                  batch()$BatchEffect()$count, batch()$BatchEffect()$meta))
          updateTabItems(session, "MenuItems", "DEAnalysis")
          
        })
        
        
        observeEvent(input$startDE, {
          if(!is.null(batch()$BatchEffect()$count)){
            waiter_show(html = waiting_screen, color = transparent(.5))
            res <- prepDataContainer(batch()$BatchEffect()$count, input, session)
            waiter_hide()
            if(is.null(res)) return(NULL)
            dc(res)
            # updateTabItems(session, "MenuItems", "DEAnalysis")
            buttonValues$startDE <- TRUE
          }
        })
        
        output$cutOffUI <- renderUI({
          cutOffSelectionUI("deresults")
        })
        
        output$ScoreCutOffUI <- renderUI({
          ScoreCutOffSelectionUI("deresults")
        })
        
        output$deresUI <- renderUI({
          getDEResultsUI("deresults")
        })

        ## Main plots UI ####
        selectedMain <- reactiveVal()
        observe({
          selectedMain(callModule(dprofilermainplot, "deresults", dc()))
        })
        
        ## Heatmaps UI #### 
        selectedHeat <- reactiveVal()
        observe({
          if (!is.null(selectedMain()) && !is.null(selectedMain()$selGenes())) {
            withProgress(message = 'Creating plot', style = "notification", value = 0.1, {
              selectedHeat(callModule(dprofilerheatmap, "deresults", dc()$init_dedata[selectedMain()$selGenes(), dc()$cols]))
            })
          }
        })
        
        selgenename <- reactiveVal()
        observe({
          if (!is.null(selectedMain()) && !is.null(selectedMain()$shgClicked())
              && selectedMain()$shgClicked()!=""){
            selgenename(selectedMain()$shgClicked())
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
        
        ## Bar Main and BoxMain Plots UI ####
        observe({
          if (!is.null(selgenename()) && selgenename()!=""){
            withProgress(message = 'Creating Bar/Box plots', style = "notification", value = 0.1, {
              callModule(dprofilerbarmainplot, "deresults", dc()$init_dedata, 
                         dc()$cols, dc()$conds, selgenename())
              callModule(dprofilerboxmainplot, "deresults", dc()$init_dedata, 
                         dc()$cols, dc()$conds, selgenename())
            })
          }
        })
        
        ## Cellular Heterogeneity Events ####
        observeEvent (input$gotodeconvolute, {
          if(!is.null(dc())){
            sel(callModule(dprofilercondselect, "deconvolute", uploadeddata()$load()$sc_count))
            updateTabItems(session, "MenuItems", "CellComp")
          }
        })
        
        observeEvent(input$deconvolute, {
          if(!is.null(dc())){
            deconvolute <- prepDeconvolute(dc, uploadeddata()$load()$sc_count, session)
          } 
        })
        
        output$cellcompUI <- renderUI({
          getDeconvoluteUI("deconvolute")
        })
        
        ## Profiling Events ####
        observeEvent (input$gotoprofile, {
          if(!is.null(dc())){
            sel(callModule(dprofilercondselect, "profiling", uploadeddata()$load()$prof_count, uploadeddata()$load()$prof_meta))
            updateTabItems(session, "MenuItems", "Profile")
          }
        })
        
        observeEvent(input$startprofiling, {
          if(!is.null(dc())){
            waiter_show(html = waiting_screen, color = transparent(.5))
            profiling <- callModule(dprofilerprofiling, "profiling", dc(), uploadeddata()$load()$prof_count,
                                    uploadeddata()$load()$prof_meta)
            waiter_hide()
          }
        })
        
        output$ProfilingUI <- renderUI({
          getProfilingUI("profiling")
        })
        
        
      })
      
      
      # auxiliary values of buttons and reactive values
      buttonValues <- reactiveValues(goQCplots = FALSE, goDE = FALSE,
                                     startDE = FALSE)
      
      # count data and metadata files for help section
      output$metaFile <-  renderTable({
        read.delim(system.file("extdata", "www", "metaFile.txt",
                               package = "debrowser"), header=TRUE, skipNul = TRUE)
      })
      output$countFile <-  renderTable({
        read.delim(system.file("extdata", "www", "countFile.txt",
                               package = "debrowser"), header=TRUE, skipNul = TRUE)
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