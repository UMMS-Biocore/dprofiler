#' debrowserdataload
#'
#' Module to load count data and metadata
#' 
#' @param input, input variables
#' @param output, output objects
#' @param session, session 
#' @param nextpagebutton, the name of the next page button after loading the data
#' @param parent_session a parameter to pass the session
#' @return main plot
#'
#' @return panel
#'
dataLoadServer <- function(input = NULL, output = NULL, session = NULL, nextpagebutton = NULL, parent_session = NULL) {
  if (is.null(input)) return(NULL)
  ldata <- reactiveValues(count=NULL, meta=NULL)
  loadeddata <- reactive({
    ret <- NULL
    if(!is.null(ldata$count)){
      ldata$count <- ldata$count[,sapply(ldata$count, is.numeric)]
      ret <- list(count = ldata$count, meta = ldata$meta)
    }
    return(ret)
  })
  output$dataloaded <- reactive({
    return(!is.null(loadeddata()))
  })
  
  observe({
    query <- parseQueryString(session$clientData$url_search)
    jsonobj<-query$jsonobject

    if (!is.null(jsonobj))
    {
      raw <- RCurl::getURL(jsonobj, .opts = list(ssl.verifypeer = FALSE),
                           crlf = TRUE)
      jsondata<-data.frame(fromJSON(raw, simplifyDataFrame = TRUE),
                           stringsAsFactors = TRUE)
      rownames(jsondata)<-jsondata[, 1]
      jsondata<-jsondata[,c(3:ncol(jsondata))]
      jsondata[,c(1:ncol(jsondata))] <- sapply(
        jsondata[,c(1:ncol(jsondata))], as.numeric)
      jsondata <- jsondata[,sapply(jsondata, is.numeric)]
      ldata$count <- jsondata
      
      metadatatable <- NULL
      jsonmet <-query$meta
      if(!is.null(jsonmet)){
        raw <- RCurl::getURL(jsonmet, .opts = list(ssl.verifypeer = FALSE),
                             crlf = TRUE)
        metadatatable<-data.frame(fromJSON(raw, simplifyDataFrame = TRUE),
                                  stringsAsFactors = TRUE)
        
      }else{
        metadatatable <- cbind(colnames(ldata$count), 1)
        colnames(metadatatable) <- c("Sample", "Batch")
      }
      ldata$meta <- metadatatable
      input$Filter
    }
  })
  
  # Event for uploading the demo file
  observeEvent(input$demo3, {
    load("demo/demodata3.Rda")
    ldata$count <- demodata
    ldata$meta <- metadatatable
    updateTabsetPanel(session = parent_session, "UploadBox", "uploadsummary")
  })
  
  # Event for uploading any file
  observeEvent(input$uploadFile, {
    if (is.null(input$countdata)) return (NULL)
    checkRes <- checkCountData(input)
    
    if (checkRes != "success"){
      showNotification(checkRes, type = "error")
      return(NULL)
    }
    counttable <-as.data.frame(
      try(
        read.delim(input$countdata$datapath, 
                   header=T, sep=input$countdataSep, 
                   row.names=1, strip.white=TRUE ), TRUE))
    counttable <- counttable[,sapply(counttable, is.numeric)]
    metadatatable <- c()
    if (!is.null(input$metadata$datapath)){
      metadatatable <- as.data.frame(
        try(
          read.delim(input$metadata$datapath, 
                     header=TRUE, sep=input$metadataSep, strip.white=TRUE), TRUE))
      
      checkRes <- checkMetaData(input, counttable)
      if (checkRes != "success"){
        showNotification(checkRes, type = "error")
        return(NULL)
      }
    }
    else{
      metadatatable <- cbind(colnames(counttable), 1)
      colnames(metadatatable) <- c("Sample", "Batch")
    }
    if (is.null(counttable)) 
    {stop("Please upload the count file")}
    ldata$count <- counttable
    ldata$meta <- metadatatable
  })
  
  # generate the button the next page
  output$nextButton <- renderUI({
    actionButtonDE(nextpagebutton, label = nextpagebutton, styleclass = "primary")
  })
  
  observe({
    getSampleDetails(output, "uploadSummary", "sampleDetails", loadeddata())
  })
  list(load=loadeddata)
}

#' dataLoadUI
#' 
#' Creates a panel to upload the data
#'
#' @param id, namespace id
#' @return panel
#' @examples
#'     x <- dataLoadUI("load")
#'
#' @export
#'
dataLoadUI<- function (id) {
  ns <- NS(id)
  list(conditionalPanel(condition =  paste0("!output['", ns("dataloaded"),"']"), fluidRow(
    fileUploadBox(id, "countdata", "Count Data"),
    fileUploadBox(id, "metadata", "Metadata")
  ),
  fluidRow(column(12,
                  actionButtonDE(ns("uploadFile"), label = "Upload", styleclass = "primary"), 
                  actionButtonDE(ns("demo3"),  label = "Load Demo PRJNA554241", styleclass = "primary"))))
  )
}

#' dataSummaryUI
#' 
#' Creates a panel to view the summary of uploaded data
#'
#' @param id, namespace id
#' @return panel
#' @examples
#'     x <- dataLoadUI("load")
#'
#' @export
#'
dataSummaryUI<- function(id) {
  ns <- NS(id)
  list(
  fluidRow(
    shinydashboard::box(title = "Upload Summary", solidHeader = TRUE, status = "info",
                        width = 12, 
                        fluidRow(
                          column(12, 
                                 tableOutput(ns("uploadSummary"))
                          )),
                        fluidRow(
                          column(12,div(style = 'overflow: scroll', 
                                        DT::dataTableOutput(ns("sampleDetails")))
                          )
                        ),
                        fluidRow(
                          column(12,uiOutput(ns("nextButton"))) 
                        )
    )
  ))
}