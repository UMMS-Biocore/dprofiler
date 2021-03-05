#' dataLoadServer
#'
#' Module to load Bulk
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
  
  # reactive values of count and metadata
  ldata <- reactiveValues(count=NULL, meta=NULL, 
                          sc_count=NULL, sc_count_UMILabel = NULL, sc_count_Cluster_Label = NULL,
                          prof_count=NULL, prof_meta=NULL)
  
  # If data sets are loaded, react
  loadeddata <- reactive({
    ret <- NULL
    
    # check if count or meta data
    if(!is.null(ldata$count)){
      ldata$count <- ldata$count[,sapply(ldata$count, is.numeric)]
      
      if(!is.null(ldata$prof_count)){
        ldata$prof_count <- ldata$prof_count[,sapply(ldata$count, is.numeric)]
      }
      
      # collect data and metadata
      ret <- list(count = ldata$count, meta = ldata$meta, 
                  sc_count = ldata$sc_count,
                  prof_count = ldata$prof_count, prof_meta = ldata$prof_meta)
      
      # switch to upload summary page
      updateTabsetPanel(session = parent_session, "UploadBox", "uploadsummary")
    }
    return(ret)
  })
  
  # Event for uploading the demo file
  observeEvent(input$demo3, {
    load("demo/demodata3.Rda")
    ldata$count <- demodata
    ldata$meta <- metadatatable
    ldata$sc_count <- demoscdata
    rm(demoscdata)
    ldata$prof_count <- demoprofdata
    ldata$prof_meta <- profmetadatatable
  })
  
  # Event for uploading any file
  observeEvent(input$uploadFile, {
    if (is.null(input$countdata)) return (NULL)
    
    ###
    # check count data and import 
    ###
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
    metadatatable <- c()
    if (!is.null(input$metadata$datapath)){
      metadatatable <- as.data.frame(
        try(
          read.delim(input$metadata$datapath, 
                     header=TRUE, sep=input$metadataSep, strip.white=TRUE), TRUE))
    }
    else{
      metadatatable <- cbind(colnames(counttable), 1)
      colnames(metadatatable) <- c("Sample", "Batch")
    }
    ldata$count <- counttable
    ldata$meta <- metadatatable
    
    ###
    # check sc count data and import 
    ###
    ldata$sc_count <- readRDS(input$sccountdata)
    colnames_sc_count <- colnames(pData(ldata$sc_count))
    colnames_sc_count[colnames_sc_count == input$sccountdataumilabel] <- "UMI_sum_raw"
    colnames_sc_count[colnames_sc_count == input$sccountdataclusterlabel] <- "CellType"
    colnames_sc_count[colnames_sc_count == input$sccountdatasamplelabel] <- "Patient"
    colnames(pData(ldata$sc_count)) <- colnames_sc_count
    
    ###
    # check prof data and import 
    ###
    profcounttable <-as.data.frame(
      try(
        read.delim(input$profilecountdata$datapath,
                   header=T, sep=input$profilecountdataSep,
                   row.names=1, strip.white=TRUE ), TRUE))
    profmetadatatable <- c()
    if (!is.null(input$profilemetadata$datapath)){
      profmetadatatable <- as.data.frame(
        try(
          read.delim(input$profilemetadata$datapath,
                     header=TRUE, sep=input$profilemetadataSep, strip.white=TRUE), TRUE))
    }
    else{
      profmetadatatable <- cbind(colnames(profcounttable), 1)
      colnames(profmetadatatable) <- c("Sample", "Batch")
    }
    ldata$prof_count <- profcounttable
    ldata$prof_meta <- profmetadatatable
    
    # check if count table is null
    if(is.null(counttable)) 
        top("Please upload the count file")
    
  })
  
  # generate the button the next page
  output$nextButton <- renderUI({
    actionButtonDE(nextpagebutton, label = nextpagebutton, styleclass = "primary")
  })
  
  observe({
    getSampleDetails(output, "uploadSummary", "sampleDetails", loadeddata())
    getProfileSampleDetails(output, "uploadSummaryprof", "sampleDetailsprof", loadeddata())
    getSCRNASampleDetails(output, "uploadSummarysc", "sampleDetailssc", loadeddata())
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
  list(
    fluidRow(
      column(12,
             actionButtonDE(ns("uploadFile"), label = "Upload", styleclass = "primary"), 
             actionButtonDE(ns("demo3"),  label = "Load Demo PRJNA554241", styleclass = "primary"))
    ),
    fluidRow(
      fileUploadBox(id, "countdata", "Bulk Count Data"),
      fileUploadBox(id, "metadata", "Bulk Metadata")
    ),
    fluidRow(
      fileUploadBox(id, "profilecountdata", "Profiling Count Data"),
      fileUploadBox(id, "profilemetadata", "Profiling Metadata")
    ),
    fluidRow(
      scfileUploadBox(id, "sccountdata", "scRNA Count Data Object"),
    )
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
    column(12,uiOutput(ns("nextButton"))),
    shinydashboard::box(title = "Bulk Upload Summary", solidHeader = TRUE, status = "info",
                        width = 6, 
                        fluidRow(
                          column(12, 
                                 tableOutput(ns("uploadSummary"))
                          )),
                        fluidRow(
                          column(12,div(style = 'overflow: scroll', 
                                        DT::dataTableOutput(ns("sampleDetails")))
                          )
                        )
    ),
    shinydashboard::box(title = "Profile Upload Summary", solidHeader = TRUE, status = "info",
                        width = 6, 
                        fluidRow(
                          column(12, 
                                 tableOutput(ns("uploadSummaryprof"))
                          )),
                        fluidRow(
                          column(12,div(style = 'overflow: scroll', 
                                        DT::dataTableOutput(ns("sampleDetailsprof")))
                          )
                        )
    )
  ),
  fluidRow(
    shinydashboard::box(title = "scRNA Upload Summary", solidHeader = TRUE, status = "info",
                        width = 12, 
                        fluidRow(
                          column(12, 
                                 tableOutput(ns("uploadSummarysc"))
                          )
                        ),
                        fluidRow(
                          column(6,div(style = 'overflow: scroll', 
                                       plotOutput(ns("sampleDetailsscdensity")))
                          ),
                          column(6,div(style = 'overflow: scroll', 
                                       plotOutput(ns("sampleDetailssctsne")))
                          )
                        )
    )
  )
  )
}

getSCRNASampleDetails <- function (output = NULL, summary = NULL, details = NULL, data = NULL) 
{
  if (is.null(data)) 
    return(NULL)
  output[[summary]] <- renderTable({
    countdata <- data$sc_count
    samplenums <- dim(countdata)[2]
    rownums <- dim(countdata)[1]
    patientnums <- length(unique(pData(countdata)$Patient))
    result <- rbind(samplenums, rownums, patientnums)
    rownames(result) <- c("# of barcodes", "# of rows (genes/regions)", "# of samples")
    colnames(result) <- "Value"
    result
  }, digits = 0, rownames = TRUE, align = "lc")
  output[[paste0(details,"density")]] <- renderPlot({
    plot_density_ridge(data$sc_count, color_by = "CellType", title = "UMIs per Cell Type", val = "UMI_sum_raw") + 
      xlim(0,7500)
  })
  output[[paste0(details,"tsne")]] <- renderPlot({
    plot_tsne_metadata(data$sc_count, color_by = "CellType", title = "Cell Types",
                       legend_dot_size = 4, text_sizes = c(20, 10, 5, 10, 15, 15))
  })
}

getProfileSampleDetails <- function (output = NULL, summary = NULL, details = NULL, data = NULL) 
{
  if (is.null(data)) 
    return(NULL)
  output[[summary]] <- renderTable({
    countdata <- data$prof_count
    samplenums <- dim(countdata)[2]
    rownums <- dim(countdata)[1]
    result <- rbind(samplenums, rownums)
    rownames(result) <- c("# of samples", "# of rows (genes/regions)")
    colnames(result) <- "Value"
    result
  }, digits = 0, rownames = TRUE, align = "lc")
  output[[details]] <- DT::renderDataTable({
    dat <- colMeans(data$prof_count)
    dat <- cbind(names(dat), dat)
    dat[, c("dat")] <- format(round(as.numeric(dat[,c("dat")], digits = 2), digits = 2), big.mark = ",", 
                              scientific = FALSE)
    if (!is.null(data$prof_meta)) {
      met <- data$prof_meta
      dat <- cbind(met, dat[, "dat"])
      rownames(dat) <- NULL
      colnames(dat)[ncol(dat)] <- "read Means"
    }
    else {
      rownames(dat) <- NULL
      colnames(dat) <- c("samples", "read Means")
    }
    dat
  })
}

sepRadio <- function (id, name) 
{
  ns <- NS(id)
  radioButtons(inputId = ns(name), label = "Separator", 
               choices = c(Comma = ",", Semicolon = ";", Tab = "\t"), selected = "\t",
               inline=T)
}

fileUploadBox <- function (id = NULL, inputId = NULL, label = NULL) 
{
  ns <- NS(id)
  shinydashboard::box(title = paste0(label, " File"), 
                      solidHeader = TRUE, status = "info", width = 6, 
                      helpText(paste0("Upload your '", label, " File'")), 
                      fileInput(inputId = ns(inputId), label = NULL, accept = fileTypes()), 
                      sepRadio(id, paste0(inputId, "Sep")))
}

scfileUploadBox <- function (id = NULL, inputId = NULL, label = NULL) 
{
  ns <- NS(id)
  shinydashboard::box(title = paste0(label, " File"), 
                      solidHeader = TRUE, status = "info", width = 6, 
                      helpText(paste0("Upload your '", label, " File'")), 
                      fileInput(inputId = ns(inputId), label = NULL, accept = fileTypes()),
                      fluidRow(
                        column(6,textInput(inputId = ns(paste0(inputId,"umilabel")), label = "UMI Column Name", value = "UMI_sum_raw")),
                        column(6,textInput(inputId = ns(paste0(inputId,"clusterlabel")), label = "Cell Type Column Name", value = "CellType"))
                      ),
                      fluidRow(
                        column(6,textInput(inputId = ns(paste0(inputId,"samplelabel")), label = "Samples Column Name", value = "Patient")),
                      )
                      
  )
}