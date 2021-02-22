#' debrowsercondselect
#'
#' Condition selection
#' This is not a module. Module construction didn't used here, just use it 
#' as functions not in a module.
#' 
#' @param input, input variables
#' @param output, output objects
#' @param session, session 
#' @param data, count data
#' @param metadata, metadata
#' @return main plot
#'
#' @return panel
#' @export
#'
#' @examples
#'     x <- debrowsercondselect()
#'
debrowsercondselect <- function(input = NULL, output = NULL, session = NULL, data = NULL, metadata = NULL) {
  if (is.null(data)) return(NULL)
  choicecounter <- reactiveVal(0)
  output$conditionSelector <- renderUI({
    selectConditions(data, metadata, choicecounter(), session, input)
  })
  
  observeEvent(input$add_btn, {
    choicecounter(choicecounter() + 1)
  })
  observeEvent(input$rm_btn, {
    if (choicecounter() > 0) 
      choicecounter(choicecounter() - 1)
  })
  list(cc = choicecounter)
}

#' condSelectUI
#' Creates a panel to select samples for each condition
#'
#' @return panel
#' @examples
#'     x <- condSelectUI()
#'
#' @export
#'
condSelectUI<- function () {
  list(
      shinydashboard::box(title = "Comparison Selection",
                          solidHeader = TRUE, status = "info",  width = NULL, height = NULL, collapsible = TRUE,
                          fluidRow(
                            uiOutput("conditionSelector"),
                            column(12,actionButtonDE("add_btn", "Add New Comparison",styleclass = "primary"),
                                   actionButtonDE("rm_btn", "Remove", styleclass = "primary"),
                                   getHelpButton("method", "http://debrowser.readthedocs.io/en/master/deseq/deseq.html"),
                                   # conditionalPanel(condition = ("output.condReady>0"),
                                   #                  actionButtonDE("startDE", "Start DE", styleclass = "primary"))
                                   actionButtonDE("startDE", "Start DE", styleclass = "primary")
                                   )
                          ))
  )
}


#' selectConditions
#'
#' Selects user input conditions, multiple if present, to be
#' used in DESeq.
#'
#' @param Dataset, used dataset 
#' @param metadata, metadatatable to select from metadata
#' @param choicecounter, choicecounter to add multiple comparisons
#' @param session, session
#' @param input, input params
#' @note \code{selectConditions}
#' @return the panel for go plots;
#'
#' @examples
#'     x<- selectConditions()
#'
#' @export
#'
selectConditions<-function(Dataset = NULL,
                           metadata = NULL,
                           choicecounter = NULL,
                           session = NULL,
                           input = NULL) {
  if (is.null(Dataset)) return(NULL)
  
  selectedSamples <- function(num){
    if (is.null(input[[paste0("condition", num)]]))
      getSampleNames(colnames(Dataset), num %% 2 )
    else
      input[[paste0("condition", num)]]
  }
  nc <- choicecounter
  
  if (nc >= 0) {
    allsamples <- getSampleNames( colnames(Dataset), "all" )
    
    lapply(seq_len(nc), function(i) {
      
      selected1 <- selectedSamples(2 * i - 1)
      selected2 <- selectedSamples( 2 * i )
      to_return <- list(column(12, getMetaSelector(metadata = metadata, input=input, n = i),
                               getConditionSelectorFromMeta(metadata, input, i,
                                                            (2 * i - 1), allsamples, selected1),
                               getConditionSelectorFromMeta(metadata, input, i,
                                                            (2 * i), allsamples, selected2)
      ),
      column(12, 
             column(1, helpText(" ")),
             getSelectInputBox("demethod", "DE Method", i, 
                               c("DESeq2", "EdgeR", "Limma"),
                               selectedInput("demethod", i, "DESeq2", input)),
             getMethodDetails(i, input)),
      column(12, 
             column(1, helpText(" ")),
             getSelectInputBox("scoremethod", "Score Method", i, 
                               c("Silhouette", "NNLS-based"),
                               selectedInput("scoremethod", i, "Silhouette", input)))
      )
      if (!is.null(selectedInput("conditions_from_meta", 
                                 i, NULL, input)) && selectedInput("conditions_from_meta", 
                                                                   i, NULL, input) != "No Selection"){
        facts <- levels(factor(metadata[,selectedInput("conditions_from_meta", 
                                                       i, NULL, input)]))
        facts <- facts[facts != "" & facts != "NA"]
        if (length(facts) != 2) {
          showNotification("There must be exactly 2 groups in the selected condition. 
                         Please use NA or space to remove extra sample groups from metadata selection.", 
                           type = "error")
          updateSelectInput(session, paste0("conditions_from_meta", i), selected="No Selection" )
        }
      }
      return(to_return)
    })
  }
}

#' prepDataContainer
#'
#' Prepares the data container that stores values used within DESeq.
#'
#' @param data, loaded dataset
#' @param counter, the number of comparisons
#' @param input, input parameters
#' @return data
#' @export
#'
#' @examples
#'     x <- prepDataContainer()
#'
prepDataContainer <- function(data = NULL, counter=NULL, 
                              input = NULL) {
  if (is.null(data)) return(NULL)
  
  inputconds <- reactiveValues(demethod_params = list(), conds = list(), dclist = list())
  iterinputconds <- reactiveValues(demethod_params = list(), conds = list(), dclist = list())
  
  inputconds$conds <- list()
  for (cnt in seq(1:(2*counter))){
    inputconds$conds[cnt] <- list(isolate(input[[paste0("condition",cnt)]]))
  }
  #Get parameters for each method
  inputconds$demethod_params <- NULL
  for (cnt in seq(1:counter)){
    if (isolate(input[[paste0("demethod",cnt)]]) == "DESeq2"){
      inputconds$demethod_params[cnt] <- paste(
        isolate(input[[paste0("demethod",cnt)]]),
        isolate(input[[paste0("fitType",cnt)]]),
        isolate(input[[paste0("betaPrior",cnt)]]),
        isolate(input[[paste0("testType",cnt)]]),
        isolate(input[[paste0("shrinkage",cnt)]]), 
        isolate(input[[paste0("scoremethod",cnt)]]), sep=",")
    }
    else if (isolate(input[[paste0("demethod",cnt)]]) == "EdgeR"){
      inputconds$demethod_params[cnt]<- paste(
        isolate(input[[paste0("demethod",cnt)]]),
        isolate(input[[paste0("edgeR_normfact",cnt)]]),
        isolate(input[[paste0("dispersion",cnt)]]),
        isolate(input[[paste0("edgeR_testType",cnt)]]), 
        isolate(input[[paste0("scoremethod",cnt)]]), sep=",")
    }
    else if (isolate(input[[paste0("demethod",cnt)]]) == "Limma"){
      inputconds$demethod_params[cnt] <- paste(
        isolate(input[[paste0("demethod",cnt)]]),
        isolate(input[[paste0("limma_normfact",cnt)]]),
        isolate(input[[paste0("limma_fitType",cnt)]]),
        isolate(input[[paste0("normBetween",cnt)]]), 
        isolate(input[[paste0("scoremethod",cnt)]]), sep=",")
    }
  }
  
  for (i in seq(1:counter))
  {
    conds <- c(rep(paste0("Cond", 2*i-1), 
                   length(inputconds$conds[[2*i-1]])), 
               rep(paste0("Cond", 2*i), length(inputconds$conds[[2*i]])))
    cols <- c(paste(inputconds$conds[[2*i-1]]), 
              paste(inputconds$conds[[2*i]]))
    params <- unlist(strsplit(inputconds$demethod_params[i], ","))
    withProgress(message = 'Running DE Algorithms', detail = inputconds$demethod_params[i], value = 0, {
      initd <- callModule(debrowserdeanalysis, paste0("DEResults",i), data = data, 
                          columns = cols, conds = conds, params = params)
      if (!is.null(initd$dat()) && nrow(initd$dat()) > 1){
        inputconds$dclist[[i]] <- list(conds = conds, cols = cols, init_data=initd$dat(), 
                                       demethod_params = inputconds$demethod_params[i])
      }else{
        return(NULL)
      }
      iterinitd <- callModule(debrowseriterdeanalysis, paste0("DEResults",i), data = data, 
                          columns = cols, conds = conds, params = params)
      if (!is.null(iterinitd$dat()) && nrow(iterinitd$dat()) > 1){
        iterinputconds$dclist[[i]] <- list(conds = conds, cols = cols, init_data=iterinitd$dat(), 
                                       demethod_params = inputconds$demethod_params[i])
      }else{
        return(NULL)
      }
      incProgress(1/counter)
    })
  }
  
  if(length(inputconds$dclist) <1) return(NULL)
  
  return(list(de = inputconds$dclist, iterde = iterinputconds$dclist))
}