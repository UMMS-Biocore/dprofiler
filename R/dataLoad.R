#' dataLoadServer
#'
#' Data loading module for bulk and single cell data sets.
#' 
#' @param input input variables
#' @param output output objects
#' @param session session 
#' @param nextpagebutton the name of the next page button after loading the data
#' @param parent_session the parent session
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
    load(system.file("extdata", "demo", "demovitiligo.Rda",
                     package = "dprofiler"))
    ldata$sc_count <- demovitiligoscdata
    ldata$count <- demodata
    ldata$meta <- metadatatable
    ldata$prof_count <- demoprofdata
    ldata$prof_meta <- profmetadatatable
  })
  
  # Cell Type Selector Menu for single cell data
  output$identSelector <- renderUI({
    
    if(!is.null(ldata$sc_count)){
      numeric_columns <- colnames(pData(ldata$sc_count))[sapply(as.data.frame(pData(ldata$sc_count)),is.numeric)]
      character_columns <- colnames(pData(ldata$sc_count))[!sapply(as.data.frame(pData(ldata$sc_count)),is.numeric)]
      list(
        column(2,
               selectInput(session$ns("selectident"), label = "Select Identification", 
                           choices = character_columns, selected = "CellType")),
        column(2,
               selectInput(session$ns("selectumi"), label = "Select UMI column", 
                           choices = numeric_columns, selected = numeric_columns[1]))
      )
    }
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
    if(!is.null(input$profiledmeta)){
      
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
    getSCRNASampleDetails(output, "uploadSummarysc", "sampleDetailssc", loadeddata()$sc_count, 
                          input$selectident, input$selectumi)
    getProfileSampleDetails(output, "uploadSummaryprof", "sampleDetailsprof", loadeddata())
  })
  
  list(load=loadeddata)
}

#' dataLoadUI
#' 
#' Creates panels to upload the data. Adapted from debrowser::dataLoadUI(). 
#'
#' @param id namespace id
#' 
#' @examples
#'     x <- dataLoadUI("load")
#'
dataLoadUI<- function (id) {
  ns <- NS(id)
  list(
    fluidRow(
      column(12,
             column(12,
                    actionButtonDE(ns("uploadFile"), label = "Upload", styleclass = "primary"), 
                    actionButtonDE(ns("demo"),  label = "Load Demo Vitiligo", styleclass = "primary")))
    ),
    fluidRow(
      column(6,
             fileUploadBox(id, "countdata", "metadata", "Bulk Expression Data"),
             scfileUploadBox(id, "sccountdata", "scRNA Expression Data Object (Optional)"),
      ),
      column(6,
             profilefileUploadBox(id, "profilecountdata", "profilemetadata", 
                                  "profiledmeta", "profiledmetakey", "Reference Bulk Expression Data (Optional)")  
      )
    )
  )
}

#' dataSummaryUI
#' 
#' Creates data table and figures to summarize uploaded data
#'
#' @param id, namespace id
#' 
#' @examples
#'     x <- dataSummaryUI("load")
#'
dataSummaryUI<- function(id) {
  ns <- NS(id)
  list(
    fluidRow(
      column(12,uiOutput(ns("nextButton")),
             p(strong("Note:")," We analyze ", strong("lesional (L) and non-lesional (NL) samples of 5 Vitiligo samples."), " We will analyze and score each sample of this dataset ", strong("to reveal critical phenotypic information"), 
               " for each individual sample. For more information: ", a("PRJNA554241",href="https://www.ncbi.nlm.nih.gov/bioproject/PRJNA554241."), " Vitiligo is an autoimmune skin disease defined by T cell–mediated destruction of melanocytes.")
      ),
      shinydashboard::box(title = "Bulk Data Summary", solidHeader = TRUE, status = "info", height = 700,
                          width = 4, 
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
      # shinydashboard::box(title = "Reference Bulk Data Summary", solidHeader = TRUE, status = "info",
      #                     width = 6, 
      #                     fluidRow(
      #                       column(12, 
      #                              tableOutput(ns("uploadSummaryprof"))
      #                       )),
      #                     fluidRow(
      #                       column(12,div(style = 'overflow: scroll', 
      #                                     DT::dataTableOutput(ns("sampleDetailsprof")))
      #                       )
      #                     )
      # )
      shinydashboard::box(title = "scRNA Data Summary", solidHeader = TRUE, status = "info", height = 700,
                          width = 8, 
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

#' dataProfileSummaryUI
#' 
#' Creates data table and figures to summarize uploaded data
#'
#' @param id, namespace id
#' 
#' @examples
#'     x <- dataProfileSummaryUI("load")
#'
dataProfileSummaryUI<- function(id) {
  ns <- NS(id)
  list(
    fluidRow(
      shinydashboard::box(title = "Reference Bulk Data Summary", solidHeader = TRUE, status = "info",
                          width = 12, 
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
    )
  )
}

#' getProfileSampleDetails
#' 
#' get Profile Data Sample Details
#'
#' @param output output
#' @param summary summary output name
#' @param details details output name
#' @param data data
#'
#' @examples
#'      x <- getProfileSampleDetails()
#'      
getProfileSampleDetails <- function (output = NULL, summary = NULL, details = NULL, data = NULL) 
{
  if (is.null(output)) 
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
      dat <- met
      # dat <- cbind(met, dat[, "dat"])
      # rownames(dat) <- NULL
      # colnames(dat)[ncol(dat)] <- "read Means"
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
#' @param ident column in scRNA metadata to visualize 
#' @param UMI_column column in scRNA metadata with total UMI counts
#'
#' @examples
#'       x <- getSCRNASampleDetails()
#'       
getSCRNASampleDetails <- function (output = NULL, summary = NULL, details = NULL, data = NULL, 
                                   ident = NULL, UMI_column = NULL) 
{
  if (is.null(output)) 
    return(NULL)
  
  if(is.null(ident) | is.null(UMI_column))
    return(NULL)
  
  if(!is.null(data)){
    
    if(any(!c(UMI_column,ident) %in% colnames(pData(data))))
      return(NULL)
    
    tsne_data <- pData(data)[,c(ident,"x","y")]
    tsne_data <- as_tibble(tsne_data) %>% group_by_at(1) %>% summarize(x = mean(x), y = mean(y))
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
      SignallingSingleCell::plot_density_ridge(data, color_by = ident, title = UMI_column, val = UMI_column) + 
        theme(legend.position = "none")
    })
    output[[paste0(details,"tsne")]] <- renderPlot({
      SignallingSingleCell::plot_tsne_metadata(data, color_by = ident,
                                               legend_dot_size = 4, text_sizes = c(20, 10, 5, 8, 8, 15)) + 
        geom_text(tsne_data, mapping = aes(x = x, y = y, label = tsne_data[[1]]), size=5)
    })
  } else {
    output[[summary]] <- renderTable({
      # countdata <- data
      # samplenums <- dim(countdata)[2]
      # rownums <- dim(countdata)[1]
      # result <- rbind(samplenums, rownums)
      # rownames(result) <- c("# of barcodes", "# of rows (genes/regions)")
      # colnames(result) <- "Value"
      # result
      
      result <- matrix("no Single Cell data was provided", nrow = 1, ncol = 1)
      colnames(result) <- "Message"
      rownames(result) <- NULL
      result 
    }, digits = 0, rownames = TRUE, align = "lc")
  }
}


#' fileUploadBox
#' 
#' file upload box for bulk data and metadata. Adapted from debrowser::fileUploadBox()
#'
#' @param id namespace id
#' @param inputId_count input data file ID
#' @param inputId_meta input metadata file ID
#' @param label label
#'
#' @examples
#'       x <- fileUploadBox("meta", "count", "metadata", "Metadata")
#' 
fileUploadBox <- function (id = NULL, inputId_count = NULL, inputId_meta = NULL, label = NULL) 
{
  ns <- NS(id)
  shinydashboard::box(title = label, 
                      solidHeader = TRUE, status = "info", width = 12, 
                      
                      helpText(paste0("Upload Data Table")), 
                      tags$div(style = "margin-bottom:-30px",
                               fileInput(inputId = ns(inputId_count), label = NULL, accept = fileTypes())
                      ),
                      tags$div(style = "margin-top:-30px;margin-bottom:30px",
                               sepRadio(id, paste0(inputId_count, "Sep"))
                      ),
                      
                      helpText(paste0("Upload MetaData Table (Optional)")), 
                      tags$div(style = "margin-bottom:-30px",
                               fileInput(inputId = ns(inputId_meta), label = NULL, accept = fileTypes())
                      ),
                      tags$div(style = "margin-top:-30px",
                               sepRadio(id, paste0(inputId_meta, "Sep"))
                      ),
  )
}

#' profilefileUploadBox
#' 
#' file upload box for profile bulk data and metadata. Adapted from debrowser::fileUploadBox()
#'
#' @param id namespace id
#' @param inputId_count input data file ID
#' @param inputId_meta input metadata file ID
#' @param inputId_dmeta dmeta API URL
#' @param inputId_dmeta_key dmeta token
#' @param label label
#'
#' @examples
#'       x <- profilefileUploadBox("meta", "count", "metadata", "Metadata")
#' 
profilefileUploadBox <- function (id = NULL, inputId_count = NULL, inputId_meta = NULL, inputId_dmeta = NULL, inputId_dmeta_key = NULL, label = NULL) 
{
  ns <- NS(id)
  shinydashboard::box(title = label, 
                      solidHeader = TRUE, status = "info", width = 12, 
                      
                      helpText(paste0("Upload Data Table")), 
                      tags$div(style = "margin-bottom:-30px",
                               fileInput(inputId = ns(inputId_count), label = NULL, accept = fileTypes())
                      ),
                      tags$div(style = "margin-top:-30px;margin-bottom:30px",
                               sepRadio(id, paste0(inputId_count, "Sep"))
                      ),
                      helpText(paste0("Upload MetaData Table (Optional)")), 
                      tags$div(style = "margin-bottom:-30px",
                               fileInput(inputId = ns(inputId_meta), label = NULL, accept = fileTypes())
                      ),
                      tags$div(style = "margin-top:-30px;margin-bottom:30px",
                               sepRadio(id, paste0(inputId_meta, "Sep"))
                      ),
                      helpText(paste0("Dmeta API of Project and Series")),
                      tags$div(style = "margin-bottom:30px",
                               textInput(inputId = ns(inputId_dmeta), label = NULL)
                      ),
                      helpText(paste0("Dmeta API key")),
                      tags$div(style = "margin-bottom:30px",
                               textInput(inputId = ns(inputId_dmeta_key), label = NULL)
                      )
  )
}

#' scfileUploadBox
#' 
#' file upload box for single cell ExpressionSet object.
#'
#' @param id namespace id
#' @param inputId input data file ID
#' @param label label
#'
#' @examples
#'       x <- fileUploadBox("meta", "metadata", "Metadata")
#'       
scfileUploadBox <- function (id = NULL, inputId = NULL, label = NULL) 
{
  ns <- NS(id)
  options(shiny.maxRequestSize=30*1024^2)
  shinydashboard::box(title = label, 
                      solidHeader = TRUE, status = "info", width = 12, 
                      helpText(paste0("Upload ExpressionSet Object (.rds)")),
                      fileInput(inputId = ns(inputId), label = NULL, accept = fileTypes())
  )
}

#' sepRadio
#'
#' Radio button for seperators: an inline version adapted from debrowser::sepRadio().
#'
#' @param id module id 
#' @param name name
#'
#' @examples
#'       x <- sepRadio("meta", "metadata")
#'       
sepRadio <- function (id, name) 
{
  ns <- NS(id)
  radioButtons(inputId = ns(name), label = "Separator", 
               choices = c(Comma = ",", Semicolon = ";", Tab = "\t"), selected = "\t",
               inline=T)
}
