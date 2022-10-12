###
# Computational Pheno. Profiling Data ####
###

#' dataLoadServer
#'
#' Data loading module for bulk and single cell data sets.
#' 
#' @param input input variables
#' @param output output objects
#' @param session session 
#' @param parent_session the parent session
#' @param mongo the data provided by mongo server module
#'
#' @return panel
#'     
#' @export
#' 
dataLoadServer <- function(input = NULL, output = NULL, session = NULL, parent_session = NULL, mongo = NULL) {
  if (is.null(input)) return(NULL)
  
  ###
  ## Reactive Values ####
  ###
  
  # reactive values of count and metadata
  ldata <- reactiveValues(count = NULL, meta = NULL, reference_path = NULL,
                          screference_path = NULL, sc_count = NULL, sc_marker_table = NULL,
                          bulkreference_path = NULL, prof_count = NULL, prof_meta = NULL, prof_marker_table = NULL, 
                          mongo = list(serverloc = NULL, mongousername = NULL, mongopassword = NULL,
                                       mongohost = NULL, mongodbname = NULL, referencerepotext = NULL))
  
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
      ret <- list(count = ldata$count, meta = ldata$meta, reference_path = ldata$reference_path, 
                  screference_path = ldata$screference_path, sc_count = ldata$sc_count, sc_marker_table = ldata$sc_marker_table,
                  bulkreference_path = ldata$bulkreference_path, prof_count = ldata$prof_count, 
                  prof_meta = ldata$prof_meta, prof_marker_table = ldata$prof_marker_table, 
                  mongo = list(serverloc = ldata$mongo$serverloc, mongousername = ldata$mongo$mongousername, 
                               mongopassword = ldata$mongo$mongopassword, mongohost = ldata$mongo$mongohost, 
                               mongodbname = ldata$mongo$mongodbname, referencerepotext = ldata$mongo$referencerepotext))
    }
    return(ret)
  })
  
  # get reactive buttons
  demoButton <- reactive(input$startdemo)
  FilterButton <- reactive(input$Filter)
  shinyjs::hide("Filter")
  
  ###
  ## Reactive Selectors ####
  ###
  
  # Condition Selector for Reference Profile Data
  output$identPCASelector <- renderUI({

    if(!is.null(ldata$prof_meta)){
      character_columns <- colnames(ldata$prof_meta)[!sapply(as.data.frame(ldata$prof_meta),is.numeric)]
      character_columns <- character_columns[sapply(as.data.frame(ldata$prof_meta[,character_columns]),function(x) length(unique(x)) < round(length(x)/4))]
      selected_columns <- character_columns[sapply(as.data.frame(ldata$prof_meta[,character_columns]),function(x) length(unique(x)) > 1 && length(unique(x)) < 5)]
      if(is.null(selected_columns)) selected_columns <- character_columns
      list(
        column(5,
               selectInput(session$ns("selectcondpca"), label = "Select Identification",
                           choices = character_columns, selected = selected_columns[1]))
      )
    }
  })
  
  # Condition Selector for Reference Profile Data
  output$identPCASelectornew <- renderUI({

    if(!is.null(ldata$prof_meta)){
      character_columns <- colnames(ldata$prof_meta)[!sapply(as.data.frame(ldata$prof_meta),is.numeric)]
      few_option_columns <- character_columns[sapply(as.data.frame(ldata$prof_meta[,character_columns]),function(x) length(unique(x)) < round(length(x)/4))]
      selected_columns <- character_columns[sapply(as.data.frame(ldata$prof_meta[,character_columns]),function(x) length(unique(x)) > 1 && length(unique(x)) < 7)]
      if(is.null(selected_columns)) selected_columns <- character_columns
      list(
        column(5,
               selectInput(session$ns("selectcondpcanew"), label = "Select Identification",
                           choices = few_option_columns, selected = selected_columns[2]))
      )
    }
  })
  
  # Cell Type Selector Menu for single cell data
  output$identSelector <- renderUI({
    
    if(!is.null(ldata$sc_count)){
      numeric_columns <- colnames(pData(ldata$sc_count))[sapply(as.data.frame(pData(ldata$sc_count)),is.numeric)]
      character_columns <- colnames(pData(ldata$sc_count))[!sapply(as.data.frame(pData(ldata$sc_count)),is.numeric)]
      list(
        column(2,
               selectInput(session$ns("selectumi"), label = "Select Ridge Variable", 
                           choices = numeric_columns, selected = numeric_columns[1])),
        column(2,
               selectInput(session$ns("selectident"), label = "Select Identification", 
                           choices = character_columns, selected = "CellType"))
      )
    }
  })
  
  ###
  ## Event for uploading the demo file ####
  ###
  
  observeEvent(input$startdemo, {
    
    if(input$demoselect == "Psoriasis"){
      
      # load demo data directory
      ldata$reference_path <- paste0(path.package("dprofiler"), "/inst/extdata/demo/Psoriasis/")
      
      # Load MongoDB Server 
      ldata$mongo$serverloc <- "Local"
      ldata$mongo$mongohost <- "0.0.0.0"
      ldata$mongo$mongodbname <- "Psoriasis"
      
      # Load Count Data 
      ldata$count <- readRDS(paste0(ldata$reference_path, "TargetData/countdata.rds"))
      ldata$meta <- readRDS(paste0(ldata$reference_path, "TargetData/metadata.rds"))
      
      # load Single Cell Reference
      ldata$screference_path <- paste0(ldata$reference_path, "SingleCellData/")
      ListOfscReferences <- list.dirs(ldata$screference_path, full.names = FALSE, recursive = FALSE)
      updateSelectInput(session, "selectscReference", choices = c("No Reference", ListOfscReferences), 
                        selected = ListOfscReferences[1])
      
      # load Bulk Reference
      ldata$bulkreference_path <- paste0(ldata$reference_path, "ReferenceBulkData/")
      ListOfbulkReferences <- list.dirs(ldata$bulkreference_path , full.names = FALSE, recursive = FALSE)
      updateSelectInput(session, "selectbulkReference", choices = c("No Reference", ListOfbulkReferences), 
                        selected = ListOfbulkReferences[1])
      
      # show start analysis button
      shinyjs::show("Filter")
      
    } else {
      
      # load demo data directory
      ldata$reference_path <- paste0(path.package("dprofiler"), "/inst/extdata/demo/KidsFirstDRC/")
      
      # Load MongoDB Server 
      ldata$mongo$serverloc <- "Local"
      ldata$mongo$mongohost <- "0.0.0.0"
      ldata$mongo$mongodbname <- "KidsFirstDRC"
      
      # Load Count Data 
      ldata$count <- readRDS(paste0(ldata$reference_path, "TargetData/countdata.rds"))
      ldata$meta <- readRDS(paste0(ldata$reference_path, "TargetData/metadata.rds"))
      
      # load Single Cell Reference
      ldata$screference_path <- paste0(ldata$reference_path, "SingleCellData/")
      ListOfscReferences <- list.dirs(ldata$screference_path, full.names = FALSE, recursive = FALSE)
      updateSelectInput(session, "selectscReference", choices = c("No Reference", ListOfscReferences), 
                        selected = ListOfscReferences[1])
      
      # load Bulk Reference
      updateSelectInput(session, "selectbulkReference", choices = c("No Reference"), 
                        selected = c("No Reference"))
      
      # show start analysis button
      shinyjs::show("Filter")
      
    }
    
    # switch to Main Sidebar Menu
    updateTabsetPanel(session = parent_session, "menutabs", "dataprep")
  })
  
  ###
  ## Event for updating reference folders ####
  ###
  
  observeEvent(mongo(), {
    
    ldata$reference_path <- mongo()
    
    # load Single Cell Reference
    ldata$screference_path <- paste0(ldata$reference_path, "/SingleCellData/")
    ListOfscReferences <- list.dirs(ldata$screference_path, full.names = FALSE, recursive = FALSE)
    updateSelectInput(session, "selectscReference", choices = c("No Reference", ListOfscReferences), 
                      selected = ListOfscReferences[1])
    
    # load Bulk Reference
    ldata$bulkreference_path <- paste0(ldata$reference_path, "/ReferenceBulkData/")
    ListOfbulkReferences <- list.dirs(ldata$bulkreference_path , full.names = FALSE, recursive = FALSE)
    updateSelectInput(session, "selectbulkReference", choices = c("No Reference", ListOfbulkReferences), 
                      selected = ListOfbulkReferences[1])
    
  })
  
  ###
  ## Event for Target Data Upload ####
  ###
  
  observeEvent(input$uploadFile, {
    if (is.null(input$countdata)) return (NULL)
    
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
    } else {
      metadatatable <- cbind(colnames(counttable), 1)
      colnames(metadatatable) <- c("Sample", "Batch")
    }
    ldata$count <- counttable
    ldata$meta <- metadatatable
    
    # show start analysis button
    shinyjs::show("Filter")
  })
  
  ###
  ### Event for Single cell data upload ####
  ###
  
  observeEvent(input$uploadscReference, {

    # if(ldata$reference_path == paste0(path.package("dprofiler"), "/inst/extdata/demo/KidsFirstDRC/") || 
    #    ldata$reference_path == paste0(path.package("dprofiler"), "/inst/extdata/demo/Psoriasis/")) {
    #   import_path <- paste0(ldata$screference_path, input$selectscReference)
    #   load(url("https://www.dropbox.com/s/5e0c5jla8za69u2/countdata.Rdata?dl=1"))
    #   ldata$sc_count <- countdata
    #   ldata$sc_marker_table <- readRDS(paste0(import_path, "/markerdata.rds"))
    # } else {
      # if scdata imported from external source, take it first!
      if(!input$externalscReference==""){
        import_path <- input$externalscReference$datapath
      } else if(!is.null(ldata$screference_path)){
        import_path <- paste0(ldata$screference_path, input$selectscReference)
      } else {
        import_path <- NULL
      }
      
      # validate the scRNA Reference by checking files
      valid <- validatescReference(import_path, session)
      
      # import sc data
      if(valid$file_check_list$count){
        if(!is.null(import_path)){
          ldata$sc_count <- readRDS(paste0(import_path, "/countdata.rds"))
          marker_file <- paste0(import_path, "/markerdata.rds")
          if(file.exists(marker_file)){
            ldata$sc_marker_table <- readRDS(marker_file)
          }
        } else {
          ldata$sc_count <- NULL
          ldata$sc_marker_table <- NULL
        }
      }
    # }
  })
  
  ###
  ### Event for Reference bulk data upload ####
  ###
    
  observeEvent(input$uploadbulkReference, {
    
    # if reference bulk imported from external source, take it first!
    if(!input$externalbulkReference==""){
      import_path <- input$externalbulkReference$datapath
    } else if(!is.null(ldata$bulkreference_path)){
      import_path <- paste0(ldata$bulkreference_path, input$selectbulkReference)
    } else {
      import_path <- NULL
    }
    
    # validate the Bulk Reference by checking files
    valid <- validatebulkReference(import_path, session)
    
    # import reference bulk data
    if(valid$file_check_list$count && valid$file_check_list$meta){
      if(!is.null(import_path)){
        ldata$prof_count <- readRDS(paste0(import_path, "/countdata.rds"))
        ldata$prof_meta <- readRDS(paste0(import_path, "/metadata.rds"))
        marker_file <- paste0(import_path, "/markerdata.rds")
        if(file.exists(marker_file)){
          ldata$prof_marker_table <- readRDS(marker_file)
        }
      } else {
        ldata$prof_count <- NULL
        ldata$prof_meta <- NULL
        ldata$prof_marker_table <- NULL
      }
    }

  })
  
  ###
  ## Event for External reference selection ####
  ###
  
  # folder selector for the reference repository
  volumes <- c(Home = fs::path_home())
  shinyDirChoose(input=input, id="screferencerepo", roots = volumes,
                 session = session)
  shinyDirChoose(input=input, id="bulkreferencerepo", roots = volumes,
                 session = session)
  
  # update Reference Repository text
  observeEvent(input$screferencerepo, {
    updateTextInput(session, "externalscReference", value = parseDirPath(roots=volumes, selection=input$screferencerepo))
  })
  observeEvent(input$bulkreferencerepo, {
    updateTextInput(session, "externalbulkReference", value = parseDirPath(roots=volumes, selection=input$bulkreferencerepo))
  })
  
  ###
  ## Main Observable ####
  ###
  
  # generate the button the next page
  # output$nextButton <- renderUI({
  #   actionButtonDE(nextpagebutton, label = nextpagebutton, styleclass = "primary")
  # })
  
  observe({

    # Sample Detail Windows
    getSampleDetails(output, "uploadSummary", "sampleDetails", loadeddata())
    getSCRNASampleDetails(output, "uploadSummarysc", "sampleDetailssc", loadeddata()$sc_count,
                          input$selectident, input$selectumi)
    getProfileSampleDetails(output, "uploadSummaryprof", "sampleDetailsprof", loadeddata(), 
                            input$selectcondpca, input$selectcondpcanew)
    
  })
  
  list(load=loadeddata, demoButton = demoButton, FilterButton = FilterButton)
}

#' MainLoadPage
#' 
#' Create tab items for the Main Page which includes the data loading modules (bulk, single cell, bulk reference), 
#' and mongo server panes 
#'     
#' @export
#' 
MainLoadPage <- function() {
  ns <- NS("load")
  list(
    tabBox(id = "AnalysisBox",
           width = NULL,
           tabPanel(title = "Dprofiler Start Up",
                    fluidRow(
                      column(2,selectInput(inputId = ns("demoselect"), label = "Select Demo",
                                           choices = c("Psoriasis", "KidsFirstDRC"))),
                      # actionButtonDE(ns("demopsoriasis"),  label = "Load Demo Psoriasis", styleclass = "primary",
                      #                style = "margin-top: 25px; margin-left: 0px"),
                      # actionButtonDE(ns("demokidsfirstdrc"),  label = "Load Demo KidsFirstDRC", styleclass = "primary",
                      #                style = "margin-top: 25px; margin-left: 0px")
                      column(1,actionButtonDE(ns("startdemo"),  label = "Load Demo", styleclass = "primary",
                                              style = "margin-top: 25px; margin-left: 0px")),
                      column(1,actionButtonDE(ns("Filter"), label = "Start Analysis", styleclass = "primary",
                                              style = "margin-top: 25px; margin-left: 0px"))
                    )
           )
    ),
    tabBox(id = "UploadBox",
           width = NULL,
           tabPanel(title = "Import Bulk RNA-Seq",
                    dataUI("load")
           ),
           tabPanel(title = "Reference scRNA-Seq",
                    dataSCUI("load"),
                    value = "uploadsummary"
           ),
           tabPanel(title = "Reference Bulk RNA-Seq",
                    dataProfileUI("load"),
                    value = "uploadprofilesummary"
           ),
           tabPanel(title = "Dprofiler Database",
                    DBUI("mongoserver"),
                    value = "dprofilerdatabase"
           )
    )
  )
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
#' @export
#' 
dataUI<- function (id) {
  ns <- NS(id)
  list(
    fluidRow(
      column(12,
                    column(1,actionButtonDE(ns("uploadFile"), label = "Import Data", styleclass = "primary",
                                   style = "margin-top: 25px; margin-left: 0px")),
                    # column(1,actionButtonDE(ns("Filter"), label = "Start Analysis", styleclass = "primary",
                    #                style = "margin-top: 25px; margin-left: 0px"))
                    # column(2,selectInput(inputId = ns("demoselect"), label = "Select Demo",
                    #                      choices = c("Psoriasis", "KidsFirstDRC"))),
                    # # actionButtonDE(ns("demopsoriasis"),  label = "Load Demo Psoriasis", styleclass = "primary",
                    # #                style = "margin-top: 25px; margin-left: 0px"),
                    # # actionButtonDE(ns("demokidsfirstdrc"),  label = "Load Demo KidsFirstDRC", styleclass = "primary",
                    # #                style = "margin-top: 25px; margin-left: 0px")
                    # column(1,actionButtonDE(ns("startdemo"),  label = "Load Demo", styleclass = "primary",
                    #                style = "margin-top: 25px; margin-left: 0px"))
      )
    ),
    fluidRow(
      column(6,
             fileUploadBox(id, "countdata", "metadata", "Bulk Expression Data"),
      ),
      column(6,
             fileSummaryBox(id, "uploadSummary", "sampleDetails", "Summary", "countdatabox"),
      )
    )
  )
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
#' @export
#' 
fileUploadBox <- function (id = NULL, inputId_count = NULL, inputId_meta = NULL, label = NULL) 
{
  ns <- NS(id)
  shinydashboard::box(title = label, 
                      solidHeader = TRUE, status = "info", width = 12, height = 420,
                      
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
                      )
  )
}

#' fileSummaryBox
#' 
#' file summary box for detailing the gene expression data and metadata
#'
#' @param id namespace id
#' @param upload upload id
#' @param sample sample id
#' @param label label
#' @param idbox id of the shiny box
#'
#' @examples
#'       x <- fileSummaryBox("meta", "count", "metadata", "Metadata")
#'     
#' @export
#' 
fileSummaryBox <- function (id = NULL, upload = NULL, sample = NULL, label = NULL, idbox = NULL) 
{
  ns <- NS(id)
  shinydashboard::box(id = ns(idbox), title = label, solidHeader = TRUE, status = "info",
                      width = 12,
                      fluidRow(
                        column(12,
                               tableOutput(ns(upload))
                        )),
                      fluidRow(
                        column(12,div(style = 'overflow: scroll',
                                      DT::dataTableOutput(ns(sample)))
                        )
                      )
  )
}

###
# Compositional Profiling Data ####
###

#' dataSCUI
#' 
#' Creates data table and figures to summarize uploaded data
#'
#' @param id, namespace id
#' 
#' @examples
#'     x <- dataSCUI("load")
#'     
#' @export
#' 
dataSCUI<- function(id) {
  ns <- NS(id)
  list(
    fluidRow(
      scfileUploadBox(id, "selectscReference", "externalscReference", "uploadscReference", "screferencerepo", "scRNA Expression Data Object (Optional)")
    ),
    fluidRow(
      shinydashboard::box(title = "Summary", solidHeader = TRUE, status = "info", # height = 700,
                          width = 12, 
                          fluidRow(
                            column(2, 
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

#' scfileUploadBox
#' 
#' file upload box for single cell ExpressionSet object.
#'
#' @param id namespace id
#' @param inputId_select input ID for selecting references
#' @param inputId_external input ID for custom external reference
#' @param inputId_button input ID for import button
#' @param inputId_dirbutton input ID for directory select button
#' @param label label
#'
#' @examples
#'       x <- scfileUploadBox("meta", "metadata", "Metadata")
#'     
#' @export
#'     
scfileUploadBox <- function (id = NULL, inputId_select = NULL, inputId_external = NULL, inputId_button = NULL, inputId_dirbutton = NULL, label = NULL) 
{
  ns <- NS(id)
  options(shiny.maxRequestSize=30*1024^2)
  list(
    column(2,disabled(textInput(inputId = ns(inputId_external), label = getIconLabel("scRNA Reference Repository", 
                                message = "This reference will be prioritized over existing references if loaded")))),
    column(1, shinyDirButton(id = ns(inputId_dirbutton), label = "Browse", title = "Select Single Cell Reference Repository",
                          class = "btn action-button btn-primary",
                          style = "width: 100%; margin-top: 25px; margin-left: 0px", buttonType = "blue")),
    column(2,selectInput(inputId = ns(inputId_select), label = getIconLabel("Select scRNA Reference", 
                         message = "This reference will be prioritized if no external repository is provided"),
                         choices = c("No Reference"))),
    column(1,actionButtonDE(inputId = ns(inputId_button),  label = "Load", styleclass = "primary", 
                            style = "width: 100%; margin-top: 25px; margin-left: 0px"))
  )
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
#' @export
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
      result <- matrix("no Single Cell data was provided", nrow = 1, ncol = 1)
      colnames(result) <- "Message"
      rownames(result) <- NULL
      result 
    }, digits = 0, rownames = TRUE, align = "lc")
  }
}


#' validatescReference
#' 
#' Validates if the single cell refernece folder has all necessary files
#'
#' @param import_path path to the single cell reference directory
#' @param session session
#'
validatescReference <- function(import_path, session){
  
  # get file names
  file_check_list <- list()
  file_check_list$count <- file.exists(paste0(import_path,"/countdata.rds"))
  file_check_list$marker <- file.exists(paste0(import_path,"/markerdata.rds"))
  
  # check if the files exists
  validity_message <- ""
  if(!file_check_list$count){
    validity_message <- "No count table is provided, the single cell reference is invalid"
  } else {
    if(!file_check_list$marker){
      validity_message <- "No marker table is provided, hence Compositional Profiling can only be conducted with all genes"
    } 
  }
  
  # show message box
  if(any(!unlist(file_check_list))){
    showModal(modalDialog(
      title = "Single Cell Reference",
      validity_message,
      easyClose = TRUE,
      footer = NULL
    ))
  }
  
  # revert back to no reference
  if(!file_check_list$count){
    updateSelectInput(session, "selectscReference", selected = "No Reference")
  }
  
  return(list(file_check_list = file_check_list, message = validity_message))
}

###
# Comparative Profiling Data ####
###

#' dataProfileUI
#' 
#' Creates data table and figures to summarize uploaded data
#'
#' @param id, namespace id
#' 
#' @examples
#'     x <- dataProfileUI("load")
#'     
#' @export
#' 
dataProfileUI<- function(id) {
  ns <- NS(id)
  list(
    fluidRow(
      profilefileUploadBox(id, "selectbulkReference", "externalbulkReference", "uploadbulkReference", "bulkreferencerepo", "Reference Bulk Expression Data (Optional)")
    ),
    fluidRow(
      shinydashboard::box(title = "PCA 1", solidHeader = TRUE, status = "info",
                          width = 6,
                          fluidRow(
                            uiOutput(ns('identPCASelector')),
                            column(12,div(style = 'overflow: scroll',
                                          plotOutput(ns("sampleDetailsprofpca1")))
                            )
                          ),
      ),
      shinydashboard::box(title = "PCA 2", solidHeader = TRUE, status = "info",
                          width = 6,
                          fluidRow(
                            uiOutput(ns('identPCASelectornew')),
                            column(12,div(style = 'overflow: scroll',
                                          plotOutput(ns("sampleDetailsprofpca2")))
                            )
                          ),
      ),
      shinydashboard::box(title = "Summary", solidHeader = TRUE, status = "info",
                          width = 12, 
                          fluidRow(
                            column(12, 
                                   tableOutput(ns("uploadSummaryprof"))
                            )),
                          fluidRow(
                            column(12,div(style = 'overflow: scroll', 
                                          DT::dataTableOutput(ns("sampleDetailsprof")))
                            ),
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
#' @param ident1 column in profile metadata for PCA plot 1
#' @param ident2 column in profile metadata for PCA plot 2
#'
#' @examples
#'      x <- getProfileSampleDetails()
#'     
#' @export
#'  
getProfileSampleDetails <- function (output = NULL, summary = NULL, details = NULL, data = NULL,
                                     ident1 = NULL, ident2 = NULL) 
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
      dat <- met
    }
    else {
      rownames(dat) <- NULL
      colnames(dat) <- c("samples", "read Means")
    }
    dat
  })
  meta.data <- data$prof_meta
  if(!is.null(ident1) && !is.null(ident2)){
    PCs <- colnames(meta.data)[grepl("^PC[1-9]",colnames(meta.data))]
    var.perc <- sapply(PCs,function(x) as.numeric(strsplit(x, split = "_")[[1]][2]))
    colnames(meta.data)[grepl("^PC[1-9]",colnames(meta.data))] <- c("PC1","PC2")
    meta.data$Group.1 <- meta.data[,ident1]
    output[[paste0(details,"pca1")]] <- renderPlot({
      ggplot(meta.data, ggplot2::aes(x = PC1, y = PC2, colour = Group.1), diag = "blank") + geom_point() + xlab(paste0("PC1 (%", var.perc[1],")")) + ylab(paste0("PC2 (%", var.perc[2],")"))
    })
    meta.data$Group.2 <- meta.data[,ident2]
    output[[paste0(details,"pca2")]] <- renderPlot({
      ggplot(meta.data, ggplot2::aes(x = PC1, y = PC2, colour = Group.2), diag = "blank") + geom_point() + xlab(paste0("PC1 (%", var.perc[1],")")) + ylab(paste0("PC2 (%", var.perc[2],")"))
    })
  }
}

#' profilefileUploadBox
#' 
#' file upload box for profile bulk data and metadata. Adapted from debrowser::fileUploadBox()
#'
#' @param id namespace id
#' @param inputId_select input ID for selecting references
#' @param inputId_external input ID for custom external reference
#' @param inputId_button input ID for import button
#' @param inputId_dirbutton input ID for directory select button
#' @param label label
#'
#' @examples
#'       x <- profilefileUploadBox("meta", "count", "metadata", "Metadata")
#'     
#' @export
#' 
profilefileUploadBox <- function (id = NULL, inputId_select = NULL, inputId_external = NULL, inputId_button = NULL, inputId_dirbutton = NULL, label = NULL) 
{
  ns <- NS(id)
  list(
    column(2,disabled(textInput(inputId = ns(inputId_external), label = getIconLabel("Bulk Reference Repository", 
                                                                                     message = "This reference will be prioritized over existing references if loaded")))),
    column(1, shinyDirButton(id = ns(inputId_dirbutton), label = "Browse", title = "Select Single Cell Reference Repository",
                             class = "btn action-button btn-primary",
                             style = "width: 100%; margin-top: 25px; margin-left: 0px", buttonType = "blue")),
    column(2,selectInput(inputId = ns(inputId_select), label = getIconLabel("Select Bulk Reference", 
                                                                            message = "This reference will be prioritized if no external repository is provided"),
                         choices = c("No Reference"))),
    column(1,actionButtonDE(inputId = ns(inputId_button),  label = "Load", styleclass = "primary", 
                            style = "width: 100%; margin-top: 25px; margin-left: 0px"))
  )
}

#' validatebulkReference
#' 
#' Validates if the bulk reference folder has all necessary files
#'
#' @param import_path path to the bulk reference directory
#' @param session session
#'
validatebulkReference <- function(import_path, session){
  
  # get file names
  file_check_list <- list()
  file_check_list$count <- file.exists(paste0(import_path,"/countdata.rds"))
  file_check_list$marker <- file.exists(paste0(import_path,"/markerdata.rds"))
  file_check_list$meta <- file.exists(paste0(import_path,"/metadata.rds"))
  
  # check if the files exists
  validity_message <- ""
  if(!file_check_list$count){
    validity_message <- "No count table is provided, the bulk reference is invalid"
  } else if(!file_check_list$meta){
    validity_message <- "No metadata table is provided, the bulk reference is invalid"
  } else {
    if(!file_check_list$marker){
      validity_message <- paste0(validity_message, "No marker table is provided, hence Compositional Profiling can only be conducted with all genes \n ")
    } 
  }
  
  # show message box
  if(any(!unlist(file_check_list))){
    showModal(modalDialog(
      title = "Bulk Reference",
      validity_message,
      easyClose = TRUE,
      footer = NULL
    ))
  }
  
  # revert back to no reference option
  if(!file_check_list$count || !file_check_list$meta){
    updateSelectInput(session, "selectbulkReference", selected = "No Reference")
  }
  
  return(list(file_check_list = file_check_list, message = validity_message))
}

###
# Auxiliary Support ####
###

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
#' @export
#'      
sepRadio <- function (id, name) 
{
  ns <- NS(id)
  radioButtons(inputId = ns(name), label = "Separator", 
               choices = c(Comma = ",", Semicolon = ";", Tab = "\t"), selected = "\t",
               inline=T)
}
