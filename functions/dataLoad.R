#' dataLoadServer
#'
#' Module to load Bulk
#' 
#' @param input, input variables
#' @param output, output objects
#' @param session, session 
#' @param nextpagebutton, the name of the next page button after loading the data
#' @param parent_session a parameter to pass the session for sending update variables to main server
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
  observeEvent(input$demo, {
    load("demo/demodata_trimsc.Rda")
    ldata$count <- demodata
    ldata$meta <- metadatatable
    ldata$sc_count <- demoscdata
    rm(demoscdata)
    ldata$prof_count <- demoprofdata
    ldata$prof_meta <- profmetadatatable
  })
  
  # Event for uploading the demo file with no single cell data
  observeEvent(input$demo_nosc, {
    load("demo/demodata_nosc.Rda")
    ldata$count <- demodata
    ldata$meta <- metadatatable
    ldata$prof_count <- demoprofdata
    ldata$prof_meta <- profmetadatatable
  })
  
  # Cell Type Selector Menu for single cell data
  output$identSelector <- renderUI({
    numeric_columns <- colnames(pData(ldata$sc_count))[sapply(as.data.frame(pData(ldata$sc_count)),is.numeric)]
    character_columns <- colnames(pData(ldata$sc_count))[!sapply(as.data.frame(pData(ldata$sc_count)),is.numeric)]
    list(
      column(2,
             selectInput(session$ns("selectident"), label = "Select Ident", 
                         choices = character_columns, selected = "CellType")),
      column(2,
             selectInput(session$ns("selectumi"), label = "Select UMI column", 
                         choices = numeric_columns, selected = "UMI_sum_raw"))
    )
  })
  
  # Event for uploading any file
  observeEvent(input$uploadFile, {
    if (is.null(input$countdata)) return (NULL)
    
    ###
    # check count data and import 
    ###
    counttable <-as.data.frame(
      try(
        read.delim(input$countdata$datapath, 
                   header=T, sep=input$countdataSep, 
                   strip.white=TRUE ), TRUE))
    counttable <- counttable[,sapply(counttable, is.numeric)]
    metadatatable <- c()
    if (!is.null(input$metadata)){
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
    else {
      metadatatable <- cbind(colnames(counttable), 1)
      colnames(metadatatable) <- c("Sample", "Batch")
    }
    ldata$count <- counttable
    ldata$meta <- metadatatable
    
    ###
    # check sc count data and import 
    ###
    if(!is.null(input$sccountdata)){
      ldata$sc_count <- readRDS(input$sccountdata$datapath)
    }
    
    ###
    # check prof data and import 
    ### 
    if(!is.null(input$profilecountdata)){
      profcounttable <-as.data.frame(
        try(
          read.delim(input$profilecountdata$datapath,
                     header=T, sep=input$profilecountdataSep,
                     strip.white=TRUE ), TRUE))
      profcounttable <- profcounttable[,sapply(profcounttable, is.numeric)]
      profmetadatatable <- c()
      if (!is.null(input$profilemetadata)){
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
    }
    
    # check if count table is null
    if(is.null(counttable)) 
        stop("Please upload the count file")
    
  })
  
  # generate the button the next page
  output$nextButton <- renderUI({
    actionButtonDE(nextpagebutton, label = nextpagebutton, styleclass = "primary")
  })
  
  observe({
    getSampleDetails(output, "uploadSummary", "sampleDetails", loadeddata())
    getProfileSampleDetails(output, "uploadSummaryprof", "sampleDetailsprof", loadeddata())
    getSCRNASampleDetails(output, "uploadSummarysc", "sampleDetailssc", loadeddata()$sc_count, 
                          input$selectident, input$selectumi)
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
             actionButtonDE(ns("demo"),  label = "Load Demo PRJNA554241", styleclass = "primary"),
             actionButtonDE(ns("demo_nosc"),  label = "Load Demo PRJNA554241 (no scRNA)", styleclass = "primary"))
    ),
    fluidRow(
      fileUploadBox(id, "countdata", "metadata", "Reference Bulk Expression Data"),
      fileUploadBox(id, "profilecountdata", "profilemetadata", "Profiling Bulk Expression Data (Optional)")
    ),
    fluidRow(
      scfileUploadBox(id, "sccountdata", "scRNA Expression Data Object (Optional)"),
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
    shinydashboard::box(title = "Reference Bulk Data Summary", solidHeader = TRUE, status = "info",
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
    shinydashboard::box(title = "Profiling Bulk Data Summary", solidHeader = TRUE, status = "info",
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
    shinydashboard::box(title = "scRNA Data Summary", solidHeader = TRUE, status = "info",
                        width = 12, 
                        fluidRow(
                          column(3, 
                                 tableOutput(ns("uploadSummarysc"))
                          ),
                          uiOutput(ns('identSelector'))
                        ),
                        fluidRow(
                          column(6,div(style = 'overflow: scroll', 
                                       plotOutput(ns("sampleDetailsscdensity")))
                          ),
                          column(6,div(style = 'overflow: scroll', 
                                       plotOutput(ns("sampleDetailssctsne")))
                          ),
                        )
    )
  )
  )
}

#' getProfileSampleDetails
#' 
#' get single cell RNA samples details
#'
#' @param output output
#' @param summary summary output name
#' @param details details output name
#' @param data data
#'
#' @return
#' @export
#'
#' @examples
#' 
getProfileSampleDetails <- function (output = NULL, summary = NULL, details = NULL, data = NULL) 
{
  if (is.null(data$prof_count)) 
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

#' getSCRNASampleDetails
#' 
#' get single cell RNA samples details
#'
#' @param output output
#' @param summary summary output name
#' @param details details output name
#' @param data single cell ExpressionSet Object
#'
#' @return
#' @export
#'
#' @examples
getSCRNASampleDetails <- function (output = NULL, summary = NULL, details = NULL, data = NULL, 
                                   ident = NULL, UMI_column = NULL) 
{
  if (is.null(data)) 
    return(NULL)
  
  if(is.null(ident) | is.null(UMI_column))
    return(NULL)
  
  output[[summary]] <- renderTable({
    countdata <- data
    samplenums <- dim(countdata)[2]
    rownums <- dim(countdata)[1]
    result <- rbind(samplenums, rownums)
    rownames(result) <- c("# of barcodes", "# of rows (genes/regions)")
    colnames(result) <- "Value"
    result
  }, digits = 0, rownames = TRUE, align = "lc")
  output[[paste0(details,"density")]] <- renderPlot({
    plot_density_ridge(data, color_by = ident, title = UMI_column, val = UMI_column)
  })
  output[[paste0(details,"tsne")]] <- renderPlot({
    plot_tsne_metadata(data, color_by = ident, title = ident,
                       legend_dot_size = 4, text_sizes = c(20, 10, 5, 10, 15, 15))
  })
}


#' fileUploadBox
#' 
#' file upload module for bulk data and metadata 
#'
#' @param id namespace id
#' @param inputId_count input data file ID
#' @param inputId_meta input metadata file ID
#' @param label label
#'
#' @return
#' @export
#'
#' @examples
fileUploadBox <- function (id = NULL, inputId_count = NULL, inputId_meta = NULL, label = NULL) 
{
  ns <- NS(id)
  shinydashboard::box(title = label, 
                      solidHeader = TRUE, status = "info", width = 6, 
                      
                      helpText(paste0("Upload Data")), 
                      tags$div(style = "margin-bottom:-30px",
                        fileInput(inputId = ns(inputId_count), label = NULL, accept = fileTypes())
                      ),
                      tags$div(style = "margin-top:-30px;margin-bottom:30px",
                               sepRadio(id, paste0(inputId_count, "Sep"))
                      ),
                      
                      helpText(paste0("Upload MetaData (Optional)")), 
                      tags$div(style = "margin-bottom:-30px",
                               fileInput(inputId = ns(inputId_meta), label = NULL, accept = fileTypes())
                      ),
                      tags$div(style = "margin-top:-30px",
                               sepRadio(id, paste0(inputId_meta, "Sep"))
                      )
                      
                      )
}

#' scfileUploadBox
#' 
#' file upload module for single cell ExpressionSet object
#'
#' @param id namespace id
#' @param inputId_count input data file ID
#' @param inputId_meta input metadata file ID
#' @param label label
#'
#' @return
#' @export
#'
#' @examples
scfileUploadBox <- function (id = NULL, inputId = NULL, label = NULL) 
{
  ns <- NS(id)
  options(shiny.maxRequestSize=30*1024^2)
  shinydashboard::box(title = label, 
                      solidHeader = TRUE, status = "info", width = 6, 
                      helpText(paste0("Upload Data")),
                      fileInput(inputId = ns(inputId), label = NULL, accept = fileTypes()),
                      # fluidRow(
                      #   column(6,textInput(inputId = ns(paste0(inputId,"umilabel")), label = "UMI Column Name", value = "UMI_sum_raw")),
                      #   column(6,textInput(inputId = ns(paste0(inputId,"clusterlabel")), label = "Cell Type Column Name", value = "CellType"))
                      # ),
                      # fluidRow(
                      #   column(6,textInput(inputId = ns(paste0(inputId,"samplelabel")), label = "Samples Column Name", value = "Patient")),
                      # )
                      
  )
}

#' sepRadio
#'
#' Radio button for seperators: an inline version of debrowser's sepRadio function
#'
#' @param id module id 
#' @param name 
#'
#' @return
#' @export
#'
#' @examples
sepRadio <- function (id, name) 
{
  ns <- NS(id)
  radioButtons(inputId = ns(name), label = "Separator", 
               choices = c(Comma = ",", Semicolon = ";", Tab = "\t"), selected = "\t",
               inline=T)
}
