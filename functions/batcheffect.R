#' debrowserbatcheffect
#'
#' Module to correct batch effect
#' 
#' @param input, input variables
#' @param output, output objects
#' @param session, session 
#' @param ldata, loaded data
#' @return main plot
#'
#' @return panel
#' @export
#'
#' @examples
#'     x <- debrowserbatcheffect()
#'
debrowserbatcheffect <- function(input, output, session, ldata = NULL) {
    if(is.null(ldata)) return(NULL)
    batchdata <- reactiveValues(count=NULL, meta = NULL)
    observeEvent(input$submitBatchEffect, {
    if (is.null(ldata$count)) return (NULL)

    countData <- ldata$count
    withProgress(message = 'Normalization', detail = "Normalization", value = NULL, {
        if (input$norm_method != "none"){
            countData <- getNormalizedMatrix(ldata$count, method=input$norm_method)
        }
    })
    withProgress(message = 'Batch Effect Correction', detail = "Adjusting the Data", value = NULL, {
    if (input$batchmethod == "Combat"){
        batchdata$count <- correctCombat(input, countData, ldata$meta)
    }
    else if (input$batchmethod == "Harman"){
        batchdata$count <- correctHarman(input, countData, ldata$meta)
    }
    else{
        batchdata$count <-  countData
    }
    })
    if (is.null(batchdata$count)) return(NULL)
    batchdata$meta <- ldata$meta
  })
  
  output$batchfields <- renderUI({
    if (!is.null(ldata$meta))
        list( conditionalPanel(condition = paste0("input['", session$ns("batchmethod"),"']!='none'"),
             selectGroupInfo( ldata$meta, input, session$ns("treatment"), "Treatment"),
             selectGroupInfo( ldata$meta, input, session$ns("batch"), "Batch")))
  })
  
  batcheffectdata <- reactive({
    ret <- NULL
    if(!is.null(batchdata$count)){
      ret <- batchdata
    }
    return(ret)
  })
  
  observe({
    getSampleDetails(output, "uploadSummary", "sampleDetails", ldata)
    getSampleDetails(output, "filteredSummary", "filteredDetails", batcheffectdata())
    getTableDetails(output, session, "beforebatchtable", ldata$count, modal = TRUE)
    callModule(debrowserpcaplot, "beforeCorrectionPCA", ldata$count, ldata$meta)
    callModule(debrowserIQRplot, "beforeCorrectionIQR",  ldata$count)
    callModule(debrowserdensityplot, "beforeCorrectionDensity", ldata$count)
    if ( !is.null(batcheffectdata()$count ) && nrow(batcheffectdata()$count)>2 ){
      withProgress(message = 'Drawing the plot', detail = "Preparing!", value = NULL, {
       getTableDetails(output, session, "afterbatchtable", batcheffectdata()$count, modal = TRUE)
       callModule(debrowserpcaplot, "afterCorrectionPCA",  batcheffectdata()$count, batcheffectdata()$meta)
       callModule(debrowserIQRplot, "afterCorrectionIQR",  batcheffectdata()$count)
       callModule(debrowserdensityplot, "afterCorrectionDensity", batcheffectdata()$count)
      })
    }
  })
  
  list(BatchEffect=batcheffectdata)
}


#' batchEffectUI
#' Creates a panel to coorect batch effect
#'
#' @param id, namespace id
#' @return panel
#' @examples
#'     x <- batchEffectUI("batcheffect")
#'
#' @export
#'
batchEffectUI <- function (id) {
  ns <- NS(id)

  list(
    fluidRow(
        shinydashboard::box(title = "Batch Effect Correction and Normalization",
        solidHeader = TRUE, status = "info",  width = 12, 
        fluidRow(
            column(5,div(style = 'overflow: scroll',
                tableOutput(ns("uploadSummary")),
                DT::dataTableOutput(ns("sampleDetails"))),
                uiOutput(ns("beforebatchtable"))
            ),
            column(2,
            shinydashboard::box(title = "Options",
                solidHeader = TRUE, status = "info",
                width = 12, 
                normalizationMethods(id),
                batchMethod(id),
                uiOutput(ns("batchfields")),
                actionButtonDE(ns("submitBatchEffect"), label = "Submit", styleclass = "primary")
           )
          ),
          column(5,div(style = 'overflow: scroll', 
                tableOutput(ns("filteredSummary")),
                DT::dataTableOutput(ns("filteredDetails"))),
                uiOutput(ns("afterbatchtable"))
          )
        ),
        conditionalPanel(condition = paste0("input['", ns("submitBatchEffect"),"']"),
        actionButtonDE("goDE", "DE Analysis", styleclass = "primary"))),
    shinydashboard::box(title = "Plots",
        solidHeader = TRUE, status = "info",  width = 12, 
        fluidRow(column(1, div()),
            tabsetPanel( id = ns("batchTabs"),
                tabPanel(id = ns("PCA"), "PCA",
                    column(5,
                        getPCAPlotUI(ns("beforeCorrectionPCA"))),
                    column(2,  
                        shinydashboard::box(title = "PCA Controls",
                        solidHeader = T, status = "info",  width = 12, 
                        tabsetPanel( id = ns("pcacontrols"),
                        tabPanel ("Before",
                        pcaPlotControlsUI(ns("beforeCorrectionPCA"))),
                        tabPanel ( "After",
                        pcaPlotControlsUI(ns("afterCorrectionPCA")))))),
                    column(5,
                        getPCAPlotUI(ns("afterCorrectionPCA")))
                ),
                tabPanel(id = ns("IQR"), "IQR",
                    column(5,
                        getIQRPlotUI(ns("beforeCorrectionIQR"))),
                    column(2, div()),
                    column(5,
                        getIQRPlotUI(ns("afterCorrectionIQR")))
                ),
                tabPanel(id = ns("Density"), "Density",
                    column(5,
                        getDensityPlotUI(ns("beforeCorrectionDensity"))),
                    column(2, div()),
                    column(5,
                        getDensityPlotUI(ns("afterCorrectionDensity")))
                )
            )
        )
      )
    ), getPCAcontolUpdatesJS())
}