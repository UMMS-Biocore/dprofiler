#' dprofilerFilter
#'
#' Module to filter low count genes/regions, adapted from debrowser::debrowserlowcountfilter
#' 
#' @param input, input variables
#' @param output, output objects
#' @param session, session 
#' @param parent_session, parent_session
#' @param uploadeddata, loaded data
#' 
#' @return main plot
#'
#' @examples
#'     x <- dprofilerFilter()
#'
dprofilerFilter <- function (input = NULL, output = NULL, session = NULL, parent_session = NULL, uploadeddata = NULL) 
{
  ldata <- uploadeddata$load
  
  ###
  # reactive Values ####
  ###
  
  fdata <- reactiveValues(count = NULL, meta = NULL)
  
  filtereddata <- reactive({
    if (!is.null(fdata$count)) {
      ret <- fdata
    } else {
      ret <- NULL
    }
    return(ret)
  })
  
  BatchButton <- reactive(input$Batch)
  DEButton <- reactive(input$goDE)
  DeconvoluteButton <- reactive(input$gotodeconvolute)
  Profile <- reactive(input$gotoprofile)

  ###
  # Filter ####
  ###
  
  observeEvent(uploadeddata$FilterButton(), {
    updateTabItems(parent_session, "MenuItems", "DataProcessing")
  })
  
  observeEvent(input$submitLCF, {
    if (is.null(ldata()$count))
      return(NULL)
    filtd <- ldata()$count
    filtd[, colnames(filtd)] <- apply(filtd[, colnames(filtd)],
                                      2, function(x) as.integer(x))
    if (input$lcfmethod == "Max") {
      filtd <- subset(filtd, apply(filtd, 1, max, na.rm = TRUE) >=
                        as.numeric(input$maxCutoff))
    }
    else if (input$lcfmethod == "Mean") {
      filtd <- subset(filtd, rowMeans(filtd, na.rm = TRUE) >=
                        as.numeric(input$meanCutoff))
    }
    else if (input$lcfmethod == "CPM") {
      cpmcount <- edgeR::cpm(filtd)
      filtd <- subset(filtd, rowSums(cpmcount > as.numeric(input$CPMCutoff),
                                     na.rm = TRUE) >= as.numeric(input$numSample))
    }
    fdata$count <- filtd
    fdata$meta <- ldata()$meta
  })
  
  output$cutoffLCFMet <- renderUI({
    ret <- textInput(session$ns("maxCutoff"), "Filter features where Max Value <",
                     value = "10")
    if (input$lcfmethod == "Mean") {
      ret <- textInput(session$ns("meanCutoff"), "Filter features where Row Means <",
                       value = "10")
    }
    else if (input$lcfmethod == "CPM") {
      ret <- list(textInput(session$ns("CPMCutoff"), "Filter features where CPM <",
                            value = "1"), textInput(session$ns("numSample"),
                                                    "at least # of samples", value = toString(ncol(ldata()$count) -
                                                                                                1)))
    }
    ret
  })

  ###
  # Main Observable ####
  ###
  
  # initially hide additional buttons
  shinyjs::hide("goDE")
  shinyjs::hide("gotodeconvolute")
  shinyjs::hide("gotoprofile")
  shinyjs::hide("Batch")
  
  # if batch correction is requested, dont let other kind of analysis
  observeEvent(input$Batch, {
    shinyjs::hide("goDE")
    shinyjs::hide("gotodeconvolute")
    shinyjs::hide("gotoprofile")
  })
  
  # if filtering is complete, show buttons
  observeEvent(input$submitLCF, {
    shinyjs::show("Batch")
    shinyjs::show("goDE")
    shinyjs::show("gotodeconvolute")
    shinyjs::show("gotoprofile")
  })
  
  # if DE analysis is requested, dont allow other methods anymore
  observeEvent(input$goDE, {
    shinyjs::hide("gotodeconvolute")
    shinyjs::hide("gotoprofile")
  })
  
  observe({
    getSampleDetails(output, "uploadSummary", "sampleDetails",
                     ldata())
    getSampleDetails(output, "filteredSummary", "filteredDetails",
                     filtereddata())
    getTableDetails(output, session, "loadedtable", data = ldata()$count,
                    modal = TRUE)
    callModule(debrowser::debrowserhistogram, "beforeFiltering", ldata()$count)
    if (!is.null(filtereddata()$count) && nrow(filtereddata()$count) >
        2) {
      getTableDetails(output, session, "filteredtable",
                      data = filtereddata()$count, modal = TRUE)
      callModule(debrowser::debrowserhistogram, "afterFiltering",
                 filtereddata()$count)
    }
  })
  list(filter = filtereddata, BatchButton = BatchButton, DEButtonFilter = DEButton, DeconvoluteFilter = DeconvoluteButton, Profile = Profile)
}

#' dprofilerBatch
#'
#' Module to correct batch effect, adapted from debrowser::debrowserbatcheffect().
#' 
#' @param input, input variables
#' @param output, output objects
#' @param session, session 
#' @param parent_session, parent session 
#' @param filtd, filtered data
#' 
#' @return panel
#' @export
#'
#' @examples
#'     x <- dprofilerBatch()
#'
dprofilerBatch <- function (input, output, session, parent_session = NULL, filtd = NULL) 
{
  ldata <- filtd$filter
  
  ###
  # reactive Values ####
  ###
  
  batchdata <- reactiveValues(count = NULL, meta = NULL)

  batcheffectdata <- reactive({
    ret <- NULL
    if (!is.null(batchdata$count)) {
      ret <- batchdata
    }
    return(ret)
  })
  
  DEButton <- reactive(input$goDE)
  DeconvoluteButton <- reactive(input$gotodeconvolute)
  Profile <- reactive(input$gotoprofile)
  
  ###
  # Batch ####
  ###
  
  observeEvent(filtd$BatchButton(), {
    showTab("DataProcessingBox","batcheffect", select = TRUE, session = parent_session)
  })
  
  observeEvent(input$submitBatchEffect, {
    if (is.null(ldata()$count))
      return(NULL)
    countData <- ldata()$count
    withProgress(message = "Normalization", detail = "Normalization",
                 value = NULL, {
                   if (input$norm_method != "none") {
                     countData <- debrowser::getNormalizedMatrix(ldata()$count,
                                                      method = input$norm_method)
                   }
                 })
    withProgress(message = "Batch Effect Correction", detail = "Adjusting the Data",
                 value = NULL, {
                   if (input$batchmethod == "CombatSeq" | input$batchmethod ==
                       "Combat") {
                     batchdata$count <- debrowser::correctCombat(input, countData,
                                                      ldata()$meta, method = input$batchmethod)
                   }
                   else if (input$batchmethod == "Harman") {
                     batchdata$count <- debrowser::correctHarman(input, countData,
                                                      ldata()$meta)
                   }
                   else {
                     batchdata$count <- countData
                   }
                 })
    if (is.null(batchdata$count))
      return(NULL)
    batchdata$meta <- ldata()$meta
  })
  output$batchfields <- renderUI({
    if (!is.null(ldata()$meta))
      list(conditionalPanel(condition = paste0("input['",
                                               session$ns("batchmethod"), "']!='none'"), debrowser::selectGroupInfo(ldata()$meta,
                                                                                                         input, session$ns("treatment"), "Treatment"),
                            debrowser::selectGroupInfo(ldata()$meta, input, session$ns("batch"),
                                            "Batch")))
  })
  
  ###
  # Main Observable ####
  ###
  
  shinyjs::hide("goDE")
  shinyjs::hide("gotodeconvolute")
  shinyjs::hide("gotoprofile")
  
  observeEvent(input$submitBatchEffect, {
    shinyjs::show("goDE")
    shinyjs::show("gotodeconvolute")
    shinyjs::show("gotoprofile")
  })
  
  # if DE analysis is requested, dont allow other methods anymore
  observeEvent(input$goDE, {
    shinyjs::hide("gotodeconvolute")
    shinyjs::hide("gotoprofile")
  })
  
  observe({
    getSampleDetails(output, "uploadSummary", "sampleDetails",
                     ldata())
    getSampleDetails(output, "filteredSummary", "filteredDetails",
                     batcheffectdata())
    getTableDetails(output, session, "beforebatchtable",
                    ldata()$count, modal = TRUE)
    callModule(debrowser::debrowserpcaplot, "beforeCorrectionPCA", ldata()$count,
               ldata()$meta)
    callModule(debrowser::debrowserIQRplot, "beforeCorrectionIQR", ldata()$count)
    callModule(debrowser::debrowserdensityplot, "beforeCorrectionDensity",
               ldata()$count)
    if (!is.null(batcheffectdata()$count) && nrow(batcheffectdata()$count) >
        2) {
      withProgress(message = "Drawing the plot", detail = "Preparing!",
                   value = NULL, {
                     getTableDetails(output, session, "afterbatchtable",
                                     batcheffectdata()$count, modal = TRUE)
                     callModule(debrowser::debrowserpcaplot, "afterCorrectionPCA",
                                batcheffectdata()$count, batcheffectdata()$meta)
                     callModule(debrowser::debrowserIQRplot, "afterCorrectionIQR",
                                batcheffectdata()$count)
                     callModule(debrowser::debrowserdensityplot, "afterCorrectionDensity",
                                batcheffectdata()$count)
                   })
    }
  })
  list(BatchEffect = batcheffectdata, DEButtonBatch = DEButton, DeconvoluteBatch = DeconvoluteButton, Profile = Profile)
}