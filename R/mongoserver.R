#' mongoServer
#' 
#' a module server to manage mongo server
#'
#' @param input input
#' @param output output
#' @param session session
#' @param uploadeddata uploaded data from initial dprofiler module
#' @param comppheno results from computational phenotypic profiling
#' @param composprof results from cmopositional profiling
#' @param comparprof results from comparative profiling
#' 
#' @examples
#'     x <- mongoServer()
#'     
#' @export
#' 
# mongoServer <- function (input = NULL, output = NULL, session = NULL, parent_session = NULL,
#                          uploadeddata = NULL, comppheno = NULL, composprof = NULL, comparprof = NULL) {
mongoServer <- function (input = NULL, output = NULL, session = NULL, parent_session = NULL,
                         uploadeddata = NULL, summarisedData = NULL) {
  
  ###
  # Reference Repository Settings ####
  ###
  
  # folder selector for the reference repository
  volumes <- c(Home = fs::path_home())
  shinyDirChoose(input=input, id="referencerepo", roots = volumes,
                 session = session)
  observeEvent(input$referencerepo, {
    updateTextInput(session, "referencerepotext", value = parseDirPath(roots=volumes, selection=input$referencerepo))
  })
  BrowsedRepo <- reactive(input$referencerepotext)
  
  ###
  # Main Mongo Server Settings ####
  ###
  
  # choose between local or remote server
  observeEvent(input$serverloc, {
    if(input$serverloc=="Local"){
      shinyjs::hide("mongousername")
      shinyjs::hide("mongopassword")
    } else {
      shinyjs::show("mongousername")
      shinyjs::show("mongopassword")
    }
  })
  
  # connect to server and get samples
  observeEvent(input$mongoconnect, {
    dprofilerMongo <- connectMongo(input, uploadeddata$load()$mongo$mongodbname)
    if(!is.null(dprofilerMongo)){
      UpdateMongoSearchBox(dprofilerMongo, input, session)
      getDataMongo(dprofilerMongo, input, output)
    } else {
      showModal(modalDialog(
        title = "MongoDB Error",
        "No server has found, please start mongodb first or enter a valid hostname",
        easyClose = TRUE,
        footer = NULL
      ))
    }
  })
  
  ###
  # Demo Settings for Mongo Server ####
  ###
  
  # connect to server and get samples for the demo from data load module
  observeEvent(uploadeddata$demoButton(), {
    updateSelectInput(session, "serverloc", selected = uploadeddata$load()$mongo$serverloc)
    updateTextInput(session, "mongohost", value = uploadeddata$load()$mongo$mongohost)
    updateTextInput(session, "mongodbname", value = uploadeddata$load()$mongo$mongodbname)
    updateTextInput(session, "referencerepotext", value = uploadeddata$load()$reference_path)
    dprofilerMongo <- connectMongo(input, uploadeddata$load()$mongo$mongodbname, uploadeddata$load()$mongo$mongohost)
    if(!is.null(dprofilerMongo)){
      if("MongoServer" %in% list.dirs(uploadeddata$load()$reference_path, full.names = FALSE)){
        mongodb_table <- readRDS(paste0(uploadeddata$load()$reference_path, "MongoServer/profiles.rds"))
        insertDataMongo(dprofilerMongo, output, mongodb_table)
        UpdateMongoSearchBox(dprofilerMongo, input, session)
      }
      getDataMongo(dprofilerMongo, input, output)
    } else {
      showModal(modalDialog(
        title = "MongoDB Error",
        "No server has found, please start mongodb first or enter a valid hostname",
        easyClose = TRUE,
        footer = NULL
      ))
    }
  })
  
  ###
  # Reactive Values for Mongo Server ####
  ###
  
  Mongo <- reactiveValues(ScoreTable = NULL)
  
  # get modal button for the Shiny Module 
  SubmitModal <- function() {
    ns <- session$ns
    modalDialog(
      tags$h2('Enter Series Name'),
      textInput(ns('newdataset'), ''),
      footer=tagList(
        actionButton(ns('submit'), 'Submit'),
        modalButton('cancel')
      )
    )
  }
  
  # only store the information if the user clicks submit
  observeEvent(input$submit, {
    removeModal()
    dprofilerMongo <- connectMongo(input)
    insertDataMongo(dprofilerMongo, output, Mongo$ScoreTable, input$newdataset)
    updateTabItems(parent_session, "MenuItems", "Upload")
    updateTabsetPanel(parent_session, "UploadBox", "dprofilerdatabase")
    UpdateMongoSearchBox(dprofilerMongo, input, session)
    getDataMongo(dprofilerMongo, input, output)
    selectDataMongo(dprofilerMongo, input, session, input$newdataset)
  })
  
  ###
  # Profiling Events for Mongo Server ####
  ###

  # get scores from the Summary module
  observeEvent(summarisedData$InsertMongo(), {
    Mongo$ScoreTable <- summarisedData$ScoreTable()
    showModal(SubmitModal())
  })
  
  return(BrowsedRepo) 
}

#' mongoServerUI
#' 
#' Creates a side panel for server set up!
#'
#' @param id namespace id
#' 
#' @examples
#'     x <- mongoServerUI()
#'     
#' @export
#' 
mongoServerUI<- function (id) {
  ns <- NS(id)
  list(
    selectInput(inputId = ns("serverloc"), label = "Server Location", choices = c("Local","Remote","Dmeta"), selected = "Local"),
    textInput(inputId = ns("mongousername"), label = "Username"),
    textInput(inputId = ns("mongopassword"), label = "Password"),
    textInput(inputId = ns("mongohost"), label = "Host (Database URL)"),
    textInput(inputId = ns("mongodbname"), label = "Database Name", value = ""),
    disabled(textInput(inputId = ns("referencerepotext"), 
                       label = getIconLabel("Reference Repository", 
                                            message = "Main directory for scRNA and Bulk RNA references"))),
    column(6,shinyDirButton(id = ns("referencerepo"), label = "Browse", title = "Select Reference Repository", 
                            class = "btn action-button btn-primary",
                            buttonType = "blue", style = "width: 100%; margin-top: 5px; margin-left: 0px; margin-right: 0px")),
    column(6,actionButtonDE(inputId = ns("mongoconnect"),  label = "Load", styleclass = "primary",
                            style = "width: 100%; margin-top: 5px; margin-left: 0px; margin-right: 0px"))
  )
}

#' DBUI
#' 
#' Creates the Mongo Database Pane
#'
#' @param id namespace id
#' 
#' @examples
#'     x <- DBUI("load")
#'     
#' @export
#' 
DBUI<- function (id) {
  ns <- NS(id)
  list(
    # database entries
    fluidRow(
      shinydashboard::box(title = "Samples and Membership Scores",
                          solidHeader = T, status = "info",  width = 12, collapsible = TRUE,
                          # uiOutput(ns("SeriesPane")),
                          column(2,selectInput(ns("SeriesPaneSelector"), "Choose Series", choices = "All Series", selected = "All Series")),
                          DT::dataTableOutput(ns("DataPane"))
      )
    )
  )
  
}

#' connectMongo
#' 
#' connect to a Mongo Server
#'
#' @param input input
#' @param database database
#' @param host host  
#'
#' @examples
#'       mongo_connection <- connectMongo(input)
#'     
connectMongo <- function (input = NULL, database = NULL, host = NULL) 
{
  if(is.null(input)) return(NULL)
  
  # get URL, choose between local and remote servers
  if(input$serverloc == "Local"){
    url <- sprintf(
      "mongodb://%s/%s",
      ifelse(is.null(host), input$mongohost, host),
      ifelse(is.null(database), input$mongodbname, database)
    )
  } else {
    
    # if there is no username or password, return NULL
    if(is.null(input$mongousername) || is.null(input$mongopassword)){
      return(NULL)
    }
    
    url = sprintf(
      "mongodb+srv://%s:%s@%s/%s",
      input$mongousername,
      input$mongopassword,
      ifelse(is.null(host), input$mongohost, host),
      ifelse(is.null(database), input$mongodbname, database)
    )
  }
  
  # connect to MongoDB server
  dprofilerMongo <- NULL
  try({
    dprofilerMongo <- mongolite::mongo(
      collection = "profiles",
      url = url,
      options = ssl_options(weak_cert_validation = TRUE)
    )
    
    # make field unique
    dprofilerMongo$run('{"createIndexes": "profiles", "indexes":[{"key":{"Sample":1}, "name": "Sample_Index", "unique": true}]}')
  })
  
  return(dprofilerMongo)
}

#' getDataMongo
#' 
#' get Phenotypic data from the MongoDB server
#'
#' @param dprofilerMongo Dprofiler Mongo connection
#' @param input input 
#' @param output output
#' @param highlight choose if you wanna highlight scores
#'
#' @examples
#'       x <- getDataMongo()
#'     
getDataMongo <- function (dprofilerMongo = NULL, input = NULL, output = NULL, highlight = TRUE) 
{
  if(is.null(dprofilerMongo)) return(NULL)
  
  # get data
  data <- dprofilerMongo$find('{}')
  
  # order columns
  if(nrow(data) > 0){
    
    # order columns
    colnames_table <- colnames(data)
    if("Series" %in% colnames_table){
      data <- data[,c("Sample","Series", colnames_table[!grepl("Sample|Series", colnames_table)])]
    } else {
      data <- data[,c("Sample", colnames_table[!grepl("Sample|Series", colnames_table)])]
    }
    
    # load data table
    output[["DataPane"]] <- DT::renderDataTable({
      if (!is.null(data)){
        if(input$SeriesPaneSelector != "All Series"){
          data <- data[data$Series == input$SeriesPaneSelector,]
        }
        dttable <- DT::datatable(data, extensions = 'Buttons',
                                 rownames = FALSE, filter = list(
                                   position = 'top', clear = FALSE
                                 ),
                                 options = list(server = TRUE,
                                                dom = "Blfrtip",
                                                buttons = list("copy", 
                                                               list(extend = "collection", 
                                                                    buttons = c("csv", "excel", "pdf"), 
                                                                    text = "Download")
                                                ), # end of buttons customization
                                                # customize the length menu
                                                lengthMenu = list(c(10, 20,  50, -1), # declare values
                                                                  c(10, 20, 50, "All") # declare titles
                                                ), # end of length Menu customization
                                                pageLength = 10))
        numeric_names <- colnames(data)[!colnames(data) %in% c("Sample","Series")]
        dttable <- dttable %>% DT::formatRound(numeric_names, digits=3)
        if(highlight){
          colours <- rainbow(length(numeric_names))
          for(i in 1:length(numeric_names)){
            dttable <-  dttable %>% DT::formatStyle(numeric_names[i],
                                                    background = DT::styleColorBar(c(0,1), colours[i]))
          }
        } 
        dttable
      }
    })
  }
}

#' UpdateMongoSearchBox
#' 
#' updates the search box of Series in the Mongo Server
#'
#' @param dprofilerMongo Dprofiler Mongo connection
#' @param output output
#' @param table MongoDB database sample table
#'
#' @examples
#'       x <- UpdateMongoSearchBox()
#'     
UpdateMongoSearchBox <- function (dprofilerMongo = NULL, input = NULL, session = NULL) 
{
  if(is.null(dprofilerMongo)) return(NULL)
  
  # get data
  data <- dprofilerMongo$find('{}')
  
  # get choices
  choices <-  c("All Series", unique(data$Series))
  selected <- "All Series"
  
  # load data table
  updateSelectInput(session, inputId = "SeriesPaneSelector", choices = choices, selected = selected)
}

#' selectDataMongo
#' 
#' select a series from the search box of Series in the Mongo Server
#'
#' @param dprofilerMongo Dprofiler Mongo connection
#' @param output output
#' @param table MongoDB database sample table
#'
#' @examples
#'       x <- selectDataMongo()
#'     
selectDataMongo <- function (dprofilerMongo = NULL, input = NULL, session = NULL, newSeries = NULL) 
{
  if(is.null(dprofilerMongo)) return(NULL)
  
  # get data
  data <- dprofilerMongo$find('{}')
  
  # get choices
  choices <-  c("All Series", unique(data$Series))
  selected <- newSeries
  
  # load data table
  updateSelectInput(session, inputId = "SeriesPaneSelector", choices = choices, selected = selected)
}

#' insertDataMongo
#' 
#' insert Phenotypic data to the MongoDB server
#'
#' @param dprofilerMongo Dprofiler Mongo connection
#' @param output output
#' @param tablename Table name in shiny 
#' @param table MongoDB database sample table
#' @param dataname Optional, ask for the name (Series) of the newly inserted dataset
#'
#' @examples
#'       x <- insertDataMongo()
#'     
insertDataMongo <- function (dprofilerMongo = NULL, output = NULL, table = NULL, dataname = NULL) 
{
  if(is.null(dprofilerMongo)) return(NULL)
  
  if(!is.data.frame(table))
    table <- as.data.frame(table)
  
  if(!is.null(dataname)){
    table$Series <- dataname
  }
  
  # order columns
  colnames_table <- colnames(table)
  if("Series" %in% colnames_table){
    table <- table[,c("Sample","Series", colnames_table[!grepl("Sample|Series", colnames_table)])]
  } else {
    table <- table[,c("Sample", colnames_table[!grepl("Sample|Series", colnames_table)])]
  }
  
  # send updated data to mongo
  samples <- table[,"Sample", drop = FALSE] 
  sample_data <- table
  for(i in 1:nrow(samples)){
    samples_json <- gsub("\\]|\\[","",jsonlite::toJSON(samples[i,,drop=FALSE]))
    sample_data_json <- paste0('{"$set":', gsub("\\]|\\[","",jsonlite::toJSON(sample_data[i,,drop=FALSE], na = "string")), '}')
    dprofilerMongo$update(samples_json, sample_data_json, upsert = TRUE)
  }
}