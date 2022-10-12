###
# Main Condition Selector Module ####
###

#' dprofilercondselect
#'
#' Condition selection for DE analysis and reference single cell data. Adapted from debrowsercondselect().
#' 
#' @param input input variables
#' @param output output objects
#' @param session session 
#' @param datax bulk or single cell data
#' @param marker_table marker table of either reference scRNA or bulk RNA data
#' @param profiling if TRUE, condition selection is conducted for external data sets (default: FALSE) 
#'
#' @examples
#'     x <- dprofilercondselect()
#'     
#' @export
#' 
# dprofilercondselect <- function(input = NULL, output = NULL, session = NULL, data = NULL, metadata = NULL, marker_table = NULL, profiling = FALSE) {
dprofilercondselect <- function(input = NULL, output = NULL, session = NULL, parent_session = NULL, data = NULL) {
  
  selectedData <- reactiveVal()
  DEparameters <- reactiveVal()
  startDE <- reactiveVal()
  startDec <- reactiveVal()
  startProf <- reactiveVal()
  
  ###
  # Bulk Data ####
  ###
  
  ###
  ## Reactive data setting for bulk data ####
  ###

  if(session$ns("")==""){
    
    PreselectedData <- reactiveValues(count = NULL, meta = NULL)
      
    observeEvent(selectedData()$DEButtonBatch(),{
      PreselectedData$count <- data$batch$BatchEffect()$count
      PreselectedData$meta <- data$batch$BatchEffect()$meta
      updateTabItems(parent_session, "MenuItems", "DEAnalysis")
    })
    observeEvent(selectedData()$DEButtonFilter(),{
      PreselectedData$count <- data$filtd$filter()$count
      PreselectedData$meta <- data$filtd$filter()$meta
      updateTabItems(parent_session, "MenuItems", "DEAnalysis")  
    })
    
    selectedData <- reactive({
      return(list(count = PreselectedData$count,
                  meta = PreselectedData$meta,
                  DEButtonBatch = data$batch$DEButtonBatch,
                  DEButtonFilter = data$filtd$DEButtonFilter))
    })
  
  
    ###
    ## Reactive Values for DE analysis ####
    ###
    
    startDE <- reactive(input$startDE)
    
    ###
    ## Reactive Values for parameters ####
    ###
    
    DEparameters <- reactiveVal()
    observeEvent(startDE(),{
      
      # get the DE parameters and pass it to the DE analysis module
      DEparameters(prepDEparameters(input, session))
      
    })
    
    ###
    ## Reactive events for bulk data ####
    ###
    
    observeEvent(input[["demethod"]],{
      if(isolate(input[["demethod"]])=="DESeq2"){
        showMulti("demethod_deseq2")
        hideMulti(c("demethod_edger","demethod_limma"))
      } else if(isolate(input[["demethod"]])=="EdgeR"){
        showMulti("demethod_edger")
        hideMulti(c("demethod_deseq2","demethod_limma"))
      } else {
        showMulti("demethod_limma")
        hideMulti(c("demethod_deseq2","demethod_edger"))
      }
    })

    observeEvent(input[["iterde"]],{
      if(isolate(input[["iterde"]])=="Log2FC+Padj"){
        showMulti("iterde_logfoldchange")
        hideMulti("iterde_topstat")
      } else {
        showMulti("iterde_topstat")
        hideMulti("iterde_logfoldchange")
      }
    })
  }
  
  ###
  # Single Cell Data ####
  ###
  
  ###
  ## Reactive data setting for single cell data ####
  ###
  
  if(grepl("deconvolute",session$ns(""))){
    
    PreselectedData <- reactiveValues(count = NULL, score = NULL)
    PreselectedscData <- reactiveValues(count = NULL, markers = NULL)
    
    observeEvent(selectedData()$DeconvoluteBatch(),{
      PreselectedData$count <- data$batch$BatchEffect()$count
      PreselectedData$score <- NULL
      PreselectedscData$count <- data$load()$sc_count
      PreselectedscData$markers <- data$load()$sc_marker_table
      updateTabItems(parent_session, "MenuItems", "CellComp")
    })
    observeEvent(selectedData()$DeconvoluteFilter(),{
      PreselectedData$count <- data$filtd$filter()$count
      PreselectedData$score <- NULL
      PreselectedscData$count <- data$load()$sc_count
      PreselectedscData$markers <- data$load()$sc_marker_table
      updateTabItems(parent_session, "MenuItems", "CellComp")
    })
    observeEvent(selectedData()$DeconvoluteDC(),{
      PreselectedData$count <- data$dc$CountData()$count[,data$cond_dc$DEparameters()$cols]
      PreselectedData$score <- data$dc$ScoreTable()
      PreselectedscData$count <- data$load()$sc_count
      PreselectedscData$markers <- data$load()$sc_marker_table
      updateTabItems(parent_session, "MenuItems", "CellComp")  
    }) 
    
    selectedData <- reactive({
      return(list(
        count =  PreselectedData$count,
        score = PreselectedData$score,
        scdata = PreselectedscData$count,
        markers = PreselectedscData$markers,
        DeconvoluteBatch = data$batch$DeconvoluteBatch,
        DeconvoluteFilter = data$filt$DeconvoluteFilter,
        DeconvoluteDC = data$dc$Deconvolute))
    })
    
    ###
    ## Reactive Values for DE analysis ####
    ###
    
    startDec <- reactive(input$deconvolute)
    
    ###
    ## Reactive events for sc data ####
    ###
    
    # Update Cell Types based on selected Metadata column
    observeEvent(input[["conditions_from_meta0"]],{
      metadata <- pData(selectedData()$scdata)
      if(!is.null(input[["conditions_from_meta0"]])){
        if(input[["conditions_from_meta0"]]!="No Selection"){
          updateSelectInput(session, "condition", choices = unique(metadata[,input[["conditions_from_meta0"]]]),
                            selected = unique(metadata[,input[["conditions_from_meta0"]]]))
        }
      }
    })
    
    # Show or Hide, Gene selection options in the scRNA condition selector 
    observeEvent(input[["allgenes"]],{
      metadata <- pData(selectedData()$scdata)
      character_columns <- colnames(metadata)[!sapply(as.data.frame(metadata),is.numeric)]
      if(isolate(input[["allgenes"]]) == "Yes"){
        if(input[["conditions_from_meta0"]]!="No Selection"){
          updateSelectInput(session, "conditions_from_meta0", choices = character_columns,
                            selected = selectedInput("conditions_from_meta", 0, NULL, input))
        }
        hideMulti(c("top_genes","logFC","padj","pct1","pct2"))
      } else {
        if(input[["conditions_from_meta0"]]!="No Selection"){
          selected <- ifelse(selectedInput("conditions_from_meta", 0, NULL, input) %in% unique(selectedData()$markers$Level),
                             selectedInput("conditions_from_meta", 0, NULL, input), 
                             unique(selectedData()$markers$Level))
          updateSelectInput(session, "conditions_from_meta0", choices = unique(selectedData()$markers$Level),
                            selected = selected)
        }
        showMulti(c("top_genes","logFC","padj","pct1","pct2"))
      }
    }) 
  }
  
  ###
  # Reference Bulk Data ####
  ###
  
  ###
  ## Reactive data setting for profiling data ####
  ###
  
  if(grepl("profiling",session$ns(""))){
    
    PreselectedData <- reactiveValues(count = NULL, score = NULL)
    PreselectedbulkData <- reactiveValues(count = NULL, meta = NULL, markers = NULL)
    
    observeEvent(selectedData()$ProfileBatch(),{
      PreselectedData$count <- data$batch$BatchEffect()$count
      PreselectedData$meta <- data$batch$BatchEffect()$meta
      PreselectedData$score <- NULL
      PreselectedbulkData$count <- data$load()$prof_count
      PreselectedbulkData$meta <- data$load()$prof_meta
      PreselectedbulkData$markers <- data$load()$prof_marker_table
      updateTabItems(parent_session, "MenuItems", "Profile")
    })
    observeEvent(selectedData()$ProfileFilter(),{
      PreselectedData$count <- data$filtd$filter()$count
      PreselectedData$meta <- data$filtd$filter()$meta
      PreselectedData$score <- NULL
      PreselectedbulkData$count <- data$load()$prof_count
      PreselectedbulkData$meta <- data$load()$prof_meta
      PreselectedbulkData$markers <- data$load()$prof_marker_table
      updateTabItems(parent_session, "MenuItems", "Profile")
    })
    observeEvent(selectedData()$ProfileDC(),{
      # metadata
      metadatatable <- data$dc$CountData()$meta
      columns <- data$cond_dc$DEparameters()$cols
      sample_column_ind <- which(apply(metadatatable, 2, function(x) sum(x %in% columns) == length(columns)))
      sample_id <- colnames(metadatatable)[sample_column_ind]
      PreselectedData$meta <- metadatatable[match(columns, metadatatable[,sample_id]),]
      
      # other samples
      PreselectedData$count <- data$dc$CountData()$count[,data$cond_dc$DEparameters()$cols]
      PreselectedData$score <- data$dc$ScoreTable()
      PreselectedbulkData$count <- data$load()$prof_count
      PreselectedbulkData$meta <- data$load()$prof_meta
      PreselectedbulkData$markers <- data$load()$prof_marker_table
      
      updateTabItems(parent_session, "MenuItems", "Profile")  
    }) 
    
    selectedData <- reactive({
      return(list(
        count =  PreselectedData$count,
        meta =  PreselectedData$meta,
        score = PreselectedData$score,
        prof_count = PreselectedbulkData$count,
        prof_meta = PreselectedbulkData$meta,
        markers = PreselectedbulkData$markers,
        ProfileBatch = data$batch$Profile,
        ProfileFilter = data$filt$Profile,
        ProfileDC = data$dc$Profile))
    })
    
    ###
    ## Reactive Values for Profiling ####
    ###
    
    startProf <- reactive(input$startprofiling)
    
  }
  
  ###
  # Main Condition Selector ####
  ###
  
  # Select Condition Selector UI based on data type
  output$conditionSelector <- renderUI({
    if(grepl("deconvolute",session$ns(""))){
      selected <- selectScRNAConditions(selectedData()$scdata, session, input)
    } 
    if(grepl("profiling",session$ns(""))){
      selected <- selectProfilingConditions(selectedData()$prof_count, selectedData()$prof_meta, selectedData()$markers, session, input)
    }
    if(session$ns("")==""){
      selected <- selectConditions(selectedData()$count, selectedData()$meta, session, input)
    }
    selected
  })
  
  # Select Covariates
  output$covariateSelector <- renderUI({
    if(session$ns("")==""){
      list(
        column(12,
               debrowser::getCovariateDetails(1, input, metadata = selectedData()$meta))
      )
    }
  })
  
  return(list(selectedData = selectedData, DEparameters = DEparameters, startDE = startDE, startDec = startDec, startProf = startProf))
}

###
# Condition Selector for Comp Pheno Profiling ####
###

#' condSelectUI
#' 
#' Creates a panel to select samples for each condition
#'
#' @examples
#'     x <- condSelectUI()
#'     
#' @export
#' 
condSelectUI<- function(){
  list(
    fluidRow(
      shinydashboard::box(title = "Comparison Selection",
                          solidHeader = TRUE, status = "info",  width = 10, height = NULL, collapsible = TRUE,
                          fluidRow(
                            uiOutput("conditionSelector"),
                            column(12, style='padding-bottom:25px;',
                                   column(2,  
                                          strong(style ='border-bottom: 1px solid #000000; padding-bottom:5px', 
                                                 "Scoring Parameters:"))
                            ),
                            column(12, 
                                   getIterMethodDetails(1, input),
                            ),
                            column(12, style='padding-bottom:25px',
                                   column(2,  
                                          strong(style ='border-bottom: 1px solid #000000; padding-bottom:5px', 
                                                 "DE Analysis Parameters:"))
                            ),
                            column(12,
                                   column(2,selectInput("demethod",
                                                  getIconLabel("DE Method", message = "Method for DE Analysis"),
                                                  c("DESeq2", "EdgeR", "Limma"), "DESeq2")),
                                   getMethodDetails(1, input),
                            ),
                            uiOutput("covariateSelector"),
                            column(12,
                                   actionButtonDE("startDE", "Start", styleclass = "primary")
                            )
                          ))
    ))
}

#' prepDataContainer 
#'
#' Prepares the data container that stores values used within DE analysis. Adapted from debrowser::prepDataContainer.
#' OLD FUNCTION DELETE LATER
#'
#' @param input, input parameters
#' @param session session
#' @param data, loaded dataset
#' @param meta, loaded metadata
#' 
#' @examples
#'     x <- prepDataContainer()
#'     
#' @export
#' 
prepDataContainer_old <- function(input = NULL, session = NULL, data = NULL, meta = NULL) {
  if (is.null(data)) return(NULL)
  
  inputconds <- reactiveValues(demethod_params = list(), conds = list(), dclist = list())
  cnt <- 1
  inputconds$conds <- list()
  inputconds$conds[[1]] <- isolate(input[[paste0("condition",1)]])
  inputconds$conds[[2]] <- isolate(input[[paste0("condition",2)]])
  
  #Get parameters for each method
  inputconds$demethod_params <- NULL
  covariate <- isolate(paste(input[[paste0("covariate",cnt)]], collapse = "|"))
  covariate <- ifelse(covariate == "", "NoCovariate", covariate)
  if (isolate(input[[paste0("demethod",cnt)]]) == "DESeq2"){
    inputconds$demethod_params[cnt] <- paste(
      isolate(input[[paste0("demethod",cnt)]]),
      covariate,
      isolate(input[[paste0("fitType",cnt)]]),
      isolate(input[[paste0("betaPrior",cnt)]]),
      isolate(input[[paste0("testType",cnt)]]),
      isolate(input[[paste0("shrinkage",cnt)]]), sep=",")
  }
  else if (isolate(input[[paste0("demethod",cnt)]]) == "EdgeR"){
    inputconds$demethod_params[cnt]<- paste(
      isolate(input[[paste0("demethod",cnt)]]),
      covariate,
      isolate(input[[paste0("edgeR_normfact",cnt)]]),
      isolate(input[[paste0("dispersion",cnt)]]),
      isolate(input[[paste0("edgeR_testType",cnt)]]), 
      "", sep=",")
  }
  else if (isolate(input[[paste0("demethod",cnt)]]) == "Limma"){
    inputconds$demethod_params[cnt] <- paste(
      isolate(input[[paste0("demethod",cnt)]]),
      covariate,
      isolate(input[[paste0("limma_normfact",cnt)]]),
      isolate(input[[paste0("limma_fitType",cnt)]]),
      isolate(input[[paste0("normBetween",cnt)]]), 
      isolate(input[[paste0("datatype",cnt)]]), sep=",")
  }
  
  # condition inputs for iterative de analysis
  inputconds$demethod_params[cnt] <- paste(inputconds$demethod_params[cnt],
                                           isolate(input[[paste0("scoremethod",cnt)]]),
                                           isolate(input[[paste0("minscore",cnt)]]),
                                           isolate(input[[paste0("iterde",cnt)]]),
                                           isolate(input[[paste0("logfoldchange",cnt)]]),
                                           isolate(input[[paste0("padj",cnt)]]),
                                           isolate(input[[paste0("topstat",cnt)]]), sep = ","
  )
  
  conds <- c(rep(paste0("Cond", 1), length(inputconds$conds[[1]])), 
             rep(paste0("Cond", 2), length(inputconds$conds[[2]])))
  cols <- c(paste(inputconds$conds[[1]]), 
            paste(inputconds$conds[[2]]))
  params <- unlist(strsplit(inputconds$demethod_params[1], ","))
  withProgress(message = 'Running Computational Profiling', 
               detail = paste0("DEmethod: ", params[1], " ScoringMethod: ", params[6]),
               value = 0, {
                 initd <- callModule(dprofilerdeanalysis, "deresults", data = data, metadata = meta, 
                                     columns = cols, conds = conds, params = params, parent_session = session)
                 if (!is.null(initd$dat()) && nrow(initd$dat()) > 1){
                   inputconds$dclist[[1]] <- list(conds = conds, 
                                                  cols = cols, 
                                                  count=initd$dat(),
                                                  DEgenes = initd$DEgenes,
                                                  count_iter = initd$iterdat(),
                                                  IterDEgenes = initd$IterDEgenes,
                                                  CrossScore = initd$CrossScore, 
                                                  IntraScore = initd$IntraScore,
                                                  ScoreTable = initd$ScoreTable,
                                                  InsertMongo = initd$InsertMongo,
                                                  cleaned_columns = initd$cleaned_columns,
                                                  demethod_params = inputconds$demethod_params[1])
                 } else {
                   return(NULL)
                 }
                 incProgress(1)
               })
  
  if(length(inputconds$dclist) <1) return(NULL)
  
  return(inputconds$dclist[[1]])
}

#' prepDEparameters
#'
#' Prepares the parameters for DE analysis. Adapted from debrowser::prepDataContainer.
#'
#' @param input, input parameters
#' @param session session
#' 
#' @examples
#'     x <- prepDEparameters()
#'     
#' @export
#' 
prepDEparameters <- function(input = NULL, session = NULL) {
  if (is.null(isolate(input$conditions_from_meta1))) return(NULL)
  
  inputconds <- list()
  cnt <- 1
  inputconds$conds <- list()
  inputconds$conds[[1]] <- isolate(input[[paste0("condition",1)]])
  inputconds$conds[[2]] <- isolate(input[[paste0("condition",2)]])
  
  #Get parameters for each method
  inputconds$demethod_params <- NULL
  covariate <- isolate(paste(input[[paste0("covariate",cnt)]], collapse = "|"))
  covariate <- ifelse(covariate == "", "NoCovariate", covariate)
  if (isolate(input[["demethod"]]) == "DESeq2"){
    inputconds$demethod_params[cnt] <- paste(
      isolate(input[["demethod"]]),
      covariate,
      isolate(input[["fitType"]]),
      isolate(input[["betaPrior"]]),
      isolate(input[["testType"]]),
      isolate(input[["shrinkage"]]), sep=",")
  }
  else if (isolate(input[["demethod"]]) == "EdgeR"){
    inputconds$demethod_params[cnt]<- paste(
      isolate(input[["demethod"]]),
      covariate,
      isolate(input[["edgeR_normfact"]]),
      isolate(input[["dispersion"]]),
      isolate(input[["edgeR_testType"]]), 
      "", sep=",")
  }
  else if (isolate(input[["demethod"]]) == "Limma"){
    inputconds$demethod_params[cnt] <- paste(
      isolate(input[["demethod"]]),
      covariate,
      isolate(input[["limma_normfact"]]),
      isolate(input[["limma_fitType"]]),
      isolate(input[["normBetween"]]), 
      isolate(input[["datatype"]]), sep=",")
  }
  
  # condition inputs for iterative de analysis
  inputconds$demethod_params[cnt] <- paste(inputconds$demethod_params[cnt],
                                           isolate(input[["scoremethod"]]),
                                           isolate(input[["minscore"]]),
                                           isolate(input[["iterde"]]),
                                           isolate(input[["logfoldchange"]]),
                                           isolate(input[["padj"]]),
                                           isolate(input[["topstat"]]), sep = ","
  )
  
  conds <- c(rep(paste0("Cond", 1), length(inputconds$conds[[1]])), 
             rep(paste0("Cond", 2), length(inputconds$conds[[2]])))
  cols <- c(paste(inputconds$conds[[1]]), 
            paste(inputconds$conds[[2]]))
  params <- unlist(strsplit(inputconds$demethod_params, ","))
  
  return(list(conds = conds, cols = cols, params = params))
}

#' selectConditions
#'
#' Selects user input conditions, multiple if present, to be used in DE analysis. 
#'
#' @param Dataset, used dataset 
#' @param metadata, metadatatable to select from metadata
#' @param session, session
#' @param input, input params
#'
#' @examples
#'     x<- selectConditions()
#'     
#' @export
#' 
selectConditions<-function(Dataset = NULL,
                           metadata = NULL,
                           session = NULL,
                           input = NULL) {
  if (is.null(Dataset)) return(NULL)
  
  selectedSamples <- function(num){
    if (is.null(input[[paste0(session$ns("condition"), num)]]))
      getSampleNames(colnames(Dataset), num %% 2 )
    else
      input[[paste0(session$ns("condition"), num)]]
  }
  
  allsamples <- getSampleNames( colnames(Dataset), "all" )
  selected1 <- selectedSamples(1)
  selected2 <- selectedSamples(2)
  to_return <- list(
    column(12, getMetaSelector(metadata = metadata, session = session, input=input, n = 1),
           debrowser::getGroupSelector(metadata, input, 1, 1),
           debrowser::getGroupSelector(metadata, input, 1, 2),
           getConditionSelectorFromMeta(metadata, session = session, input, 1,
                                        1, allsamples, selected1),
           getConditionSelectorFromMeta(metadata, session = session, input, 1,
                                        2, allsamples, selected2)
    )
    # column(12,
    #        debrowser::getCovariateDetails(1, input, metadata = metadata))
  )
  
  
  # check DE conditions
  if (!is.null(selectedInput("conditions_from_meta", 1, NULL, input)) && selectedInput("conditions_from_meta", 1, NULL, input) != "No Selection"){
    facts <- levels(factor(metadata[,selectedInput("conditions_from_meta", 1, NULL, input)]))
    facts <- facts[facts != "" & facts != "NA"]
    if (length(facts) < 2) {
      showNotification("There must be more than 2 groups in the selected condition.", 
                       type = "error")
      updateSelectInput(session, paste0("conditions_from_meta", 1), selected="No Selection" )
    }
  }
  
  # update selected of the covariate
  # if condition is selected, dont let the same column selected as covariate, then remove
  metadata_columns <- colnames(metadata)[2:ncol(metadata)]
  metadata_columns <- metadata_columns[!metadata_columns %in% selectedInput("conditions_from_meta", 1, NULL, input)]
  if(!is.null(selectedInput("conditions_from_meta", 1, NULL, input)) && !is.null(selectedInput("covariate", 1, NULL, input)) &&
     selectedInput("conditions_from_meta", 1, NULL, input) != "No Selection"){
    if(selectedInput("conditions_from_meta", 1, NULL, input) %in% selectedInput("covariate", 1, NULL, input)){
      selected_meta <- selectedInput("covariate", 1, NULL, input)
      selected_meta <- selected_meta[!selected_meta %in% selectedInput("conditions_from_meta", 1, NULL, input)]
      updateSelectInput(session, paste0("covariate", 1), choices = c(metadata_columns), 
                        selected = c(selected_meta))
      showNotification("Condition column is included in covariate list, removing!", 
                       type = "error")
      return(to_return)
    } 
  }
  
  # check appropriateness of covariates
  if (!is.null(selectedInput("covariate", 1, NULL, input))){
    
    # establish metadata with selected samples, conditions and covariates
    selected_samples <- c(selected1,selected2)
    match_selected_samples <- match(selected_samples, colnames(Dataset))
    selected_covariate_names <- selectedInput("covariate",1, NULL, input)
    selected_covariates_metadata <- metadata[match_selected_samples,selected_covariate_names, drop = FALSE]
    selected_treatment <- rep(c("Cond1","Cond2"), c(length(selected1), length(selected2)))
    
    # loop over all covariates, flag covariates to be removed
    flag_selected <- rep(F,length(selected_covariate_names))
    for(kk in 1:ncol(selected_covariates_metadata)){
      covariate <- selected_covariates_metadata[,kk]
      
      # check if there are at least two groups in metadata
      if (sum(is.na(covariate)) > 0) {
        showNotification("Covariate shouldnt have an NA or empty values in any selected sample.",
                         type = "error")
        flag_selected[kk] <- T
      }
      
      # check if there are at least two groups in metadata
      if (length(unique(covariate)) < 2) {
        showNotification("There must be at least 2 groups in the selected covariate.",
                         type = "error")
        flag_selected[kk] <- T
      }
      
      # check if there are confounding covariates with treatment
      covariate_vs_treatment <- table(covariate, selected_treatment)
      if (any(covariate_vs_treatment == 0)) {
        showNotification("Each condition should have at least one of all covariate groups.",
                         type = "error")
        flag_selected[kk] <- T
      }
    }
    
    # remove covariates if an inappropriate flag is found
    if(any(flag_selected)){
      new_covariate_names <- selected_covariate_names[!flag_selected]
      selected_covariate_names <- selectedInput("covariate",1, NULL, input)
      selected_covariate_names <- selected_covariate_names[selected_covariate_names %in% new_covariate_names]
      updateSelectInput(session, paste0("covariate", 1), selected= c(selected_covariate_names))   
    }
    
  }
  
  return(to_return)
}

#' getMetaSelector
#'
#' Return the sample selection box using meta data table.
#' 
#' @param metadata meta data table
#' @param session session
#' @param input input params 
#' @param n the box number
#'
#' @examples
#'      x <- getMetaSelector()
#'     
#' @export
#'    
getMetaSelector <- function (metadata = NULL, session = NULL, input = NULL, n = 0) 
{
  if (!is.null(metadata)) {
    df <- metadata
    col_count <- length(colnames(df))
    list(HTML("<hr style=\"color: white; border:solid 1px white;\">"), 
         br(), column(12, 
                      selectInput(paste0(session$ns("conditions_from_meta"), n), 
                                  label = "Select Meta", choices = as.list(c("No Selection", colnames(df)[2:col_count])), 
                                  multiple = FALSE, 
                                  selected = selectedInput(session$ns("conditions_from_meta"), n, "Selection 2", input))))
  }
}

#' getConditionSelectorFromMeta
#' 
#' Selects user input conditions to run in DE analysis from metadata.
#' 
#' @param metadata meta data table
#' @param session session 
#' @param input input
#' @param index index
#' @param num num
#' @param choices choices 
#' @param selected selected
#'
#' @examples
#'      x <- getConditionSelectorFromMeta()
#'     
#' @export
#'     
getConditionSelectorFromMeta <- function (metadata = NULL, session = NULL, input = NULL, index = 1, num = 0, 
                                          choices = NULL, selected = NULL) {
  if (is.null(metadata)) return(NULL)
  a <- list(column(6, selectInput(paste0(session$ns("condition"), num),
                                  label = paste0("Condition ", num), choices = choices,
                                  multiple = TRUE, selected = selected)))
  
  if (!is.null(metadata)){
    selected_meta <- selectedInput(session$ns("conditions_from_meta"), 
                                   index, NULL, input)
    
    if (is.null(selected_meta)) selected_meta <- "No Selection"
    
    if (selected_meta != "No Selection"){
      old_selection <- ""
      
      if (!is.null(input[[paste0(session$ns("condition"), num)]])) {
        selected <- input[[paste0(session$ns("condition"), num)]]
      }
      grps <- unique(metadata[selected_meta])
      grps <- grps[grps!="NA"]
      grps<- grps[!is.na(grps) ]
      if (length(grps) == 2) {
        meta_choices_all <- NULL
        if (!is.null(selected_meta))
          meta_choices_all <- get_conditions_given_selection(metadata, selected_meta)
        if(old_selection != selected_meta){
          if(typeof(meta_choices_all) == "character"){
            meta_choices <- list("There must be exactly 2 groups.")
          } else{
            meta1 <- meta_choices_all[[2 - (num %% 2)]]
            meta_choices <- unlist(meta1, recursive=FALSE)
          }
          selected <- meta_choices
        }
      }else{
        if(!is.null(input[[paste0(session$ns("group"), num)]])){
          selected <- metadata[metadata[,selected_meta] == input[[paste0(session$ns("group"), num)]], 1]
        }
      }
      
      a <- list(column(6, selectInput(paste0(session$ns("condition"), num), 
                                      label = paste0("Condition ", num), 
                                      choices = choices, multiple = TRUE, selected = selected)))
    }
  }
  return(a)
}

#  getMethodDetails
#'
#' get the detail boxes after DE method selected 
#'
#' @param num panel that is going to be shown
#' @param session session
#' @param input user input
#' 
#' @examples
#'     x <- getMethodDetails()
#'     
#' @export
#' 
getMethodDetails <- function(num = 0, input = NULL) {
  if (num > 0)
    list(
      # conditionalPanel(
      #   (condition <- paste0("input.demethod",num," == 'DESeq2'")),
      tags$div(id = "demethod_deseq2",
               column(2,selectInput("fitType", 
                                    getIconLabel("Fit Type", message = "fitting type for dispersion estimate"), 
                                    c("parametric", "local", "mean"), "parametric")),
               column(2,selectInput("betaPrior", 
                                    getIconLabel("Beta Prior", message = "Use a zero-mean normal prior. Default for Wald statistics"), 
                                    c(FALSE, TRUE), FALSE)),
               column(2,selectInput("testType", 
                                    getIconLabel("Test Type", message = "Test statistics used in DESeq2"), 
                                    c("LRT", "Wald"), "Wald")),
               column(2,selectInput("shrinkage", 
                                    getIconLabel("Shrinkage", message = "The method of shrinkage estimate for log2FC"), 
                                    c("None", "apeglm", "ashr", "normal"), "None"))
               ),
      tags$div(id = "demethod_edger",
               column(2,selectInput("edgeR_normfact", 
                                 getIconLabel("Normalization", message = "Normalization method used in EdgeR"), 
                                 c("TMM","RLE","upperquartile","none"), "TMM")),
               column(2,textInput("dispersion", "Dispersion", "0")),
               column(2,selectInput("edgeR_testType", 
                                 getIconLabel("Test Type", message = "Type of model fitting for gene abundances"), 
                                 c("exactTest", "glmLRT"), "exactTest"))
      ),
      tags$div(id = "demethod_limma",
               column(2,selectInput("limma_normfact", 
                                 getIconLabel("Normalization", message = "Normalization method used in Limma"), 
                                 c("TMM","RLE","upperquartile","none"), "TMM")),
               column(2,selectInput("limma_fitType", 
                                 getIconLabel("Fit Type", message = "Type of model fitting: 'ls' for least squares, 'robust' for robust regression"), 
                                 c("ls", "robust"), "ls")),
               column(2,selectInput("normBetween", 
                                 getIconLabel("Norm. Bet. Arrays", message = "Normalization Between Array"), 
                                 c("none", "scale", "quantile", "cyclicloess",
                                        "Aquantile", "Gquantile", "Rquantile","Tquantile"), "none")),
               column(2,selectInput("datatype", 
                                 getIconLabel("Data Type", message = "Either Count or microarray data"), 
                                 c("count", "microarray"), "count"))
      )
  )
}

#  getMethodDetails
#'
#' get the detail boxes after DE method selected 
#'
#' @param num panel that is going to be shown
#' @param session session
#' @param input user input
#' 
#' @examples
#'     x <- getMethodDetails()
#'     
#' @export
#' 
getMethodDetails_old <- function(num = 0, session = NULL, input = NULL) {
  if (num > 0)
    list(
      conditionalPanel(
        (condition <- paste0("input.demethod",num," == 'DESeq2'")),
        # ns = session$ns, 
        getSelectInputBox(session$ns("fitType"), 
                          getIconLabel("Fit Type", message = "fitting type for dispersion estimate"), 
                          num, c("parametric", "local", "mean"), 
                          selectedInput(session$ns("testType"), num, "parametric", input), 2),
        getSelectInputBox(session$ns("betaPrior"), 
                          getIconLabel("Beta Prior", message = "Use a zero-mean normal prior. Default for Wald statistics"), 
                          num, c(FALSE, TRUE), 
                          selectedInput("betaPrior", num, FALSE, input),2),
        getSelectInputBox(session$ns("testType"), 
                          getIconLabel("Test Type", message = "Test statistics used in DESeq2"), 
                          num, c("LRT", "Wald"),  
                          selectedInput(session$ns("testType"), num, "Wald", input)),
        getSelectInputBox(session$ns("shrinkage"), 
                          getIconLabel("Shrinkage", message = "The method of shrinkage estimate for log2FC"), 
                          num, c("None", "apeglm", "ashr", "normal"),
                          selectedInput(session$ns("shrinkage"), num, "None", input))),
      conditionalPanel(
        (condition <- paste0("input.demethod",num," == 'EdgeR'")),
        # ns = session$ns, 
        getSelectInputBox(session$ns("edgeR_normfact"), 
                          getIconLabel("Normalization", message = "Normalization method used in EdgeR"), 
                          num, c("TMM","RLE","upperquartile","none"), 
                          selectedInput(session$ns("edgeR_normfact"), num, "TMM", input), 2),
        column(2,textInput(paste0(session$ns("dispersion"), num), "Dispersion", 
                           value = selectedInput(session$ns("dispersion"), num, "0", input))),
        getSelectInputBox(session$ns("edgeR_testType"), 
                          getIconLabel("Test Type", message = "Type of model fitting for gene abundances"), 
                          num, c("exactTest", "glmLRT"), 
                          selectedInput(session$ns("edgeR_testType"), num, "exactTest", input))),
      conditionalPanel(
        (condition <- paste0("input.demethod",num," ==  'Limma'")),
        # ns = session$ns, 
        getSelectInputBox(session$ns("limma_normfact"), 
                          getIconLabel("Normalization", message = "Normalization method used in Limma"), 
                          num, c("TMM","RLE","upperquartile","none"), 
                          selectedInput(session$ns("limma_normfact"), num, "TMM", input), 2),
        getSelectInputBox(session$ns("limma_fitType"), 
                          getIconLabel("Fit Type", message = "Type of model fitting: 'ls' for least squares, 'robust' for robust regression"), 
                          num, c("ls", "robust"), 
                          selectedInput(session$ns("limma_fitType"), num, "ls", input)),
        getSelectInputBox(session$ns("normBetween"), 
                          getIconLabel("Norm. Bet. Arrays", message = "Normalization Between Array"), 
                          num, c("none", "scale", "quantile", "cyclicloess",
                                 "Aquantile", "Gquantile", "Rquantile","Tquantile"),
                          selectedInput(session$ns("normBetween"), num, "none", input)),
        getSelectInputBox(session$ns("datatype"), 
                          getIconLabel("Data Type", message = "Either Count or microarray data"), 
                          num, c("count", "microarray"),
                          selectedInput(session$ns("datatype"), num, "count", input))
        
      ),
      
      br())
}

#' getIterMethodDetails
#'
#' get the detail boxes after interative DE method selected 
#'
#' @param num panel that is going to be shown
#' @param session session
#' @param input user input
#' 
#' @examples
#'     x <- getIterMethodDetails()
#'     
#' @export
#' 
getIterMethodDetails <- function(num = 0, input = NULL) {
  if (num > 0)
    list(
      column(2,selectInput("scoremethod", 
                        getIconLabel("Score Method", message = "Method for estimating Membership scores"), 
                        c("Silhouette", "NNLS-based"), "Silhouette")),
      column(2,textInput("minscore", 
                         getIconLabel("Min. Score", message = "Threshold for acceptable self-membership scores. Score < Min.Score indicates dismemberment of the sample, if 'auto' then threshold is set automatically (with respect to condition imbalance)"), 
                         value = 0.5)),
      column(2,selectInput("iterde",
                        getIconLabel("DE Selection Method", message = "Set of conditions for selecting DE genes on each iteration"), 
                        c("Top n Stat.", "Log2FC+Padj"), "Log2FC+Padj")),
      tags$div(id = "iterde_topstat",
               column(2,textInput("topstat",
                                  getIconLabel("Top n Stat.", message = "DE genes with top n DE Analysis statistics. Might be different for each DE method"), 
                                  value = 100))
      ),
      tags$div(id = "iterde_logfoldchange",
               column(2,textInput("logfoldchange",
                                  getIconLabel("Log2FC", message = "Minimum log fold change"), 
                                  value = 0.5)),
               column(2,textInput("padj", 
                                  getIconLabel("P-adj value", message = "Maximum adjusted p-value"),
                                  value = 0.05))
      )
    )
}

###
# Condition Selector for Compositional Profiling ####
###

#' selectScRNAConditions
#'
#' Selects user input conditions used in bulk RNA deconvolution.
#'
#' @param scdata, used single cell dataset 
#' @param session, session
#' @param choicecounter, choicecounter to add multiple comparisons
#' @param input, input params
#'
#' @examples
#'     x<- selectScRNAConditions()
#'     
#' @export
#' 
selectScRNAConditions<-function(scdata = NULL,
                                session = NULL,
                                input = NULL) {
  if(is.null(scdata)) return(NULL)
  
  # meta data
  metadata <- pData(scdata)
  character_columns <- colnames(metadata)[!sapply(as.data.frame(metadata),is.numeric)]
  metadata <- metadata[,character_columns]
  sample_idents_ind <- grepl("sample|donor|patient", tolower(character_columns))
  sample_ident <- ifelse(any(sample_idents_ind), character_columns[sample_idents_ind][1], character_columns[1])
  
  # meta data selection pane
  to_return <- list(
    # column(12, 
    #        p(strong("Note:")," Here, bulk RNA Deconvolution is used to profiling bulk lesional and non-lesional vitiligo samples with reference skin cell types of Vitiligo (Keratinocyte, Melanocyte, t-cells etc.). ",
    #          strong("Membership scores"), " and estimated cell type proportions are concurrantly visualized.")
    # ),
    column(12,
           
           # Select Cell Types that will be used to deconvolute samples
           getScRNAMetaSelector(metadata, session, input),
           getIdentSelectorFromMeta(metadata, session, input, choices = "No Selection"),
           
           # subset single cells to deconvolute distinct samples
           # to_return_phenotype,
           # column(12,actionButtonDE(session$ns("add_btn"), "Add New Comparison",styleclass = "primary"),
           #       actionButtonDE(session$ns("rm_btn"), "Remove", styleclass = "primary")),
    ),
    column(12,
           column(2,
                  selectInput(session$ns('deconvolute_methods'), 
                              label = getIconLabel("Method", message = "select a deconvolution method"),
                              choices = c("MuSIC","BisqueRNA","SCDC"),
                              selected = "MuSIC")),
           column(2,
                  selectInput(session$ns('deconvolute_samples'), 
                              label = getIconLabel("Samples", message = "select metadata column with samples of origins of cells"),
                              choices = character_columns,
                              selected = sample_ident)),
           column(2,
                  selectInput(session$ns("deconvolute_norm"), 
                              label = getIconLabel("Normalization", message = "select a normalization method for the bulk data"),
                              choices = c("none", "CPM+TMM", "TMM", "MRN", "RLE", "upperquartile"),
                              selected = "none")),
           column(2,
                  selectInput(session$ns("allgenes"), 
                              label = getIconLabel("Use All Genes ?", message = "Do you want to use all genes available for the analysis"),
                              choices = c("Yes","No"),
                              selected = "Yes"))
    ), 
    column(12,
           # conditionalPanel(
           #   (condition <- "input.allgenes == 'No'"),
           #     ns = session$ns,
           column(2,
                  textInput(session$ns("top_genes"), 
                            label = getIconLabel("Top N Markers", message = "# of markers for each cel type, only for cell marker genes option"),
                            value = "100")),
           column(2,
                  textInput(session$ns('logFC'), 
                            label = getIconLabel("LogFC", message = "Minimum log2FC for a marker"),
                            value = "1")
           ),
           column(2,
                  textInput(session$ns('padj'), 
                            label = getIconLabel("Padj-value", message = "Padj-value"),
                            value = "0.05")
           ),
           column(2,
                  textInput(session$ns('pct1'), 
                            label = getIconLabel("Pct.1", message = "Min. Perc. of non-zero expressing cells of target cell type"),
                            value = "0.5")
           ),
           column(2,
                  textInput(session$ns('pct2'), 
                            label = getIconLabel("Pct.2", message = "Max. Perc. of non-zero expressing cells in other cell types"),
                            value = "1")
           ),
           #)
    )
  )
  
  # # Update Phenotypes
  # for(i in 1:nc){
  #   if (!is.null(selectedInput(session$ns("phenotype_from_meta"), i, NULL, input))){
  #       updateSelectInput(session, paste0(session$ns("phenotypecase_from_meta"), 1), 
  #                         choices=unique(metadata[,selectedInput(session$ns("phenotype_from_meta"), i, NULL, input)]))
  #     }
  # }
  
  return(to_return)
}

#' getScRNAMetaSelector
#'
#' Return the sample selection box using meta data table of scRNA data.
#' 
#' @param metadata meta data table
#' @param session session
#' @param input input params 
#' @param n the box number
#'
#' @examples
#'      x <- getScRNAMetaSelector()
#'     
#' @export
#'    
getScRNAMetaSelector <- function (metadata = NULL, session = NULL, input = NULL, n = 0) 
{
  if (!is.null(metadata)) {
    df <- metadata
    col_count <- length(colnames(df))
    default <- ifelse("CellType" %in% colnames(df), "CellType", "Selection 2")
    list(HTML("<hr style=\"color: white; border:solid 1px white;\">"), 
         br(), column(12, 
                      selectInput(paste0(session$ns("conditions_from_meta"), n), 
                                  label = getIconLabel("Select Annotation", message = "select metadata columns with cell types or cell annotations"),
                                  choices = as.list(c("No Selection", colnames(df)[2:col_count])), 
                                  multiple = FALSE, 
                                  selected = selectedInput(session$ns("conditions_from_meta"), n, default, input))))
  }
}

#' getIdentSelectorFromMeta
#' 
#' Select identification from single cell metadata.
#'
#' @param metadata meta data table
#' @param session session
#' @param input input 
#' @param choices choices
#'
#' @examples
#'      x <- getIdentSelectorFromMeta()
#'     
#' @export
#'      
getIdentSelectorFromMeta <- function (metadata = NULL, session = NULL, input = NULL, 
                                      choices = NULL){
  if (!is.null(metadata)) {
    selected_meta <- selectedInput(session$ns("conditions_from_meta"), 0, NULL, input)
    if (is.null(selected_meta)) 
      selected_meta <- "No Selection"
    if(selected_meta != "No Selection"){
      choices <- unique(metadata[,selected_meta])
    }
    a <- list(column(12, selectInput(session$ns("condition"), 
                                     label = getIconLabel("Identifications", message = "select cell types or cell annotations"),
                                     choices = choices, multiple = TRUE, selected = choices[1])))
  }
}

###
# Condition Selector for Comparative Profiling ####
###

#' selectProfilingConditions
#'
#' Select metadata field to use DE results from reference data
#'
#' @param Dataset, used dataset 
#' @param metadata, metadatatable to select from metadata
#' @param marker_table, the marker table of the profile data
#' @param session, session
#' @param input, input params
#'
#' @examples
#'     x<- selectProfilingConditions()
#'     
#' @export
#' 
selectProfilingConditions<-function(Dataset = NULL,
                                    metadata = NULL,
                                    marker_table = NULL,
                                    session = NULL,
                                    input = NULL) {
  if (is.null(Dataset)) return(NULL)
  
  to_return <- list(
    # column(12, 
    #        p(strong("Note:")," Here, we use similarity measures based on silhouette measure and non-negative least squares to calculate ", strong("the membership score of Vitiligo samples using another reference
    #          bulk Vitiligo dataset,"), " revealing similarities of lesional and non-lesional samples across datasets.")
    # ),
    column(6,
           getProfileSelector(metadata, marker_table, session, input)
    ),
    column(12, style='padding-bottom:25px;',
           column(2,  
                  strong(style ='border-bottom: 1px solid #000000; padding-bottom:5px', 
                         "Scoring Parameters:"))
    ),
    column(12, 
           getProfilingMethodDetails(1, session, input),
    )
  )
  
  return(to_return)
}

#' getProfileSelector
#'
#' Return the series selection box using profiling meta data table.
#' 
#' @param metadata meta data table
#' @param marker_table marker_table
#' @param session session
#' @param input input params 
#' @param n the box number
#'
#' @examples
#'      x <- getProfileSelector()
#'     
#' @export
#'     
getProfileSelector <- function (metadata = NULL, marker_table = NULL, session = NULL, input = NULL, n = 0) 
{
  if(!is.null(metadata)){
    default <- unique(marker_table$Level)[1]
    list(HTML("<hr style=\"color: white; border:solid 1px white;\">"), 
         br(), column(12, 
                      selectInput(paste0(session$ns("conditions_from_meta"), n), 
                                  label = getIconLabel("Select Meta", message = "select metadata column with reference conditions"),
                                  choices = as.list(c("No Selection", unique(marker_table$Level))), 
                                  multiple = FALSE, 
                                  selected = selectedInput(session$ns("conditions_from_meta"), n, default, input)))
    )
  }
}

#' getProfilingMethodDetails
#'
#' get the detail boxes after interative DE method selected 
#'
#' @param num panel that is going to be shown
#' @param session session
#' @param input user input
#' 
#' @examples
#'     x <- getProfilingMethodDetails()
#'     
#' @export
#' 
getProfilingMethodDetails <- function(num = 0, session = NULL, input = NULL) {
  if (num > 0)
    list(
      getSelectInputBox(session$ns("scoremethod"), 
                        getIconLabel("Score Method", message = "Method for estimating Membership scores"), 
                        num, c("Silhouette", "NNLS-based"),
                        selectedInput(session$ns("scoremethod"), num, "Silhouette", input)),
      column(2,
             selectInput(session$ns("profiling_norm"), 
                         label = getIconLabel("Normalization", message = "select a normalization method for the bulk data"),
                         choices = c("none", "MRN", "TMM", "RLE", "upperquartile"),
                         selected = "none")),
      column(2,
             textInput(session$ns("top_genes"), 
                       label = getIconLabel("Top N Markers", message = "# of markers for each cel type, only for cell marker genes option"),
                       value = "100")),
      column(2,
             textInput(session$ns('logFC'), 
                       label = getIconLabel("LogFC", message = "Minimum log2FC for a gene"),
                       value = "1")
      ),
      column(2,
             textInput(session$ns('padj'), 
                       label = getIconLabel("Padj-value", message = "Padj-value"),
                       value = "0.05")
      )
    )
}