#' dprofilercondselect
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
#'     x <- dprofilercondselect()
#'
dprofilercondselect <- function(input = NULL, output = NULL, session = NULL, data = NULL, metadata = NULL) {
  if (is.null(data)) return(NULL)
  
  output$conditionSelector <- renderUI({
    selectConditions(data, metadata, session, input)
  })
  
  return(NULL)
}

#' condSelectUI
#' 
#' Creates a panel to select samples for each condition
#'
#' @return panel
#' @examples
#'     x <- condSelectUI()
#'
#' @export
#'
condSelectUI<- function(){
  list(
      shinydashboard::box(title = "Comparison Selection",
                          solidHeader = TRUE, status = "info",  width = NULL, height = NULL, collapsible = TRUE,
                          fluidRow(
                            uiOutput("conditionSelector"),
                            column(12,
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
                           session = NULL,
                           input = NULL) {
  if (is.null(Dataset)) return(NULL)
  
  selectedSamples <- function(num){
    if (is.null(input[[paste0("condition", num)]]))
      getSampleNames(colnames(Dataset), num %% 2 )
    else
      input[[paste0("condition", num)]]
  }

  allsamples <- getSampleNames( colnames(Dataset), "all" )
  selected1 <- selectedSamples(1)
  selected2 <- selectedSamples(2)
  to_return <- list(
    column(12, getMetaSelector(metadata = metadata, input=input, n = 1),
           getConditionSelectorFromMeta(metadata, input, 1,
                                        1, allsamples, selected1),
           getConditionSelectorFromMeta(metadata, input, 1,
                                        2, allsamples, selected2)
    ),
    column(12, style='padding-bottom:15px;',
           column(12,strong("DE Analysis Parameters:"))
           ),
    column(12,
           getSelectInputBox("demethod", "DE Method", 1, 
                             c("DESeq2", "EdgeR", "Limma"),
                             selectedInput("demethod", 1, "DESeq2", input)),
           getMethodDetails(1, input)),
    column(12, style='padding-bottom:15px;',
           column(12,strong("Scoring Parameters:"))
    ),
    column(12, 
           getIterMethodDetails(1, input),
    )
  )
  
  if (!is.null(selectedInput("conditions_from_meta", 
                             1, NULL, input)) && selectedInput("conditions_from_meta", 
                                                               1, NULL, input) != "No Selection"){
    facts <- levels(factor(metadata[,selectedInput("conditions_from_meta", 
                                                   1, NULL, input)]))
    facts <- facts[facts != "" & facts != "NA"]
    if (length(facts) != 2) {
      showNotification("There must be exactly 2 groups in the selected condition. 
                         Please use NA or space to remove extra sample groups from metadata selection.", 
                       type = "error")
      updateSelectInput(session, paste0("conditions_from_meta", 1), selected="No Selection" )
    }
  }
  
  return(to_return)
}

#' getIterMethodDetails
#'
#' get the detail boxes after DE method selected 
#'
#' @param num, panel that is going to be shown
#' @param input, user input
#' @examples
#'     x <- getIterMethodDetails()
#'
#' @export
#'
#'
getIterMethodDetails <- function(num = 0, input = NULL) {
  if (num > 0)
    list(
      getSelectInputBox("scoremethod", "Score Method", num, 
                        c("Silhouette", "NNLS-based"),
                        selectedInput("scoremethod", num, "Silhouette", input)),
      column(2,textInput(paste0("minscore", num), "Min. Score", 
                         value = isolate(selectedInput("minscore", 
                                                       num, "0.5", input)))),
      getSelectInputBox("iterde_norm", "Normalization", num, 
                        c("TMM","RLE","upperquartile","none"), 
                        selectedInput("iterde_norm", num, "TMM", input), 2),
      getSelectInputBox("iterde", "DE Selection Method", num, 
                        c("Stat.", "Log2FC+Padj"),
                        selectedInput("iterde", num, "Log2FC+Padj", input)),
      conditionalPanel(
        (condition <- paste0("input.iterde",num," == 'Stat.'")),
        column(2,textInput(paste0("topstat", num), "Top Stat", 
                           value = isolate(selectedInput("topstat", 
                                                         num, "100", input) )))
      ),
      conditionalPanel(
        (condition <- paste0("input.iterde",num," == 'Log2FC+Padj'")),
        column(1,textInput(paste0("logfoldchange", num), "Log2FC",
                           value = isolate(selectedInput("logfoldchange",
                                                         num, "0.5", input) ))),
        column(1,textInput(paste0("padj", num), "Padj",
                           value = isolate(selectedInput("padj",
                                                         num, "0.05", input) )))
      )
    )
}

#' prepDataContainer
#'
#' Prepares the data container that stores values used within DESeq.
#'
#' @param data, loaded dataset
#' @param input, input parameters
#' @return data
#' @export
#'
#' @examples
#'     x <- prepDataContainer()
#'
prepDataContainer <- function(data = NULL, input = NULL) {
  if (is.null(data)) return(NULL)
  
  inputconds <- reactiveValues(demethod_params = list(), conds = list(), dclist = list())
  cnt <- 1
  inputconds$conds <- list()
  inputconds$conds[[1]] <- isolate(input[[paste0("condition",1)]])
  inputconds$conds[[2]] <- isolate(input[[paste0("condition",2)]])
  
  #Get parameters for each method
  inputconds$demethod_params <- NULL
  if (isolate(input[[paste0("demethod",cnt)]]) == "DESeq2"){
    inputconds$demethod_params[cnt] <- paste(
      isolate(input[[paste0("demethod",cnt)]]),
      isolate(input[[paste0("fitType",cnt)]]),
      isolate(input[[paste0("betaPrior",cnt)]]),
      isolate(input[[paste0("testType",cnt)]]),
      isolate(input[[paste0("shrinkage",cnt)]]), sep=",")
  }
  else if (isolate(input[[paste0("demethod",cnt)]]) == "EdgeR"){
    inputconds$demethod_params[cnt]<- paste(
      isolate(input[[paste0("demethod",cnt)]]),
      isolate(input[[paste0("edgeR_normfact",cnt)]]),
      isolate(input[[paste0("dispersion",cnt)]]),
      isolate(input[[paste0("edgeR_testType",cnt)]]), 
      "", sep=",")
  }
  else if (isolate(input[[paste0("demethod",cnt)]]) == "Limma"){
    inputconds$demethod_params[cnt] <- paste(
      isolate(input[[paste0("demethod",cnt)]]),
      isolate(input[[paste0("limma_normfact",cnt)]]),
      isolate(input[[paste0("limma_fitType",cnt)]]),
      isolate(input[[paste0("normBetween",cnt)]]), 
      "", sep=",")
  }
  
  # condition inputs for iterative de analysis
  inputconds$demethod_params[cnt] <- paste(inputconds$demethod_params[cnt],
                                           isolate(input[[paste0("scoremethod",cnt)]]),
                                           isolate(input[[paste0("minscore",cnt)]]),
                                           isolate(input[[paste0("iterde_norm",cnt)]]), 
                                           isolate(input[[paste0("iterde",cnt)]]),
                                           isolate(input[[paste0("logfoldchange",cnt)]]),
                                           isolate(input[[paste0("padj",cnt)]]),
                                           isolate(input[[paste0("topstat",cnt)]]), sep = ","
                                           )
  
  conds <- c(rep(paste0("Cond", 1), 
                 length(inputconds$conds[[1]])), 
             rep(paste0("Cond", 2), length(inputconds$conds[[2]])))
  cols <- c(paste(inputconds$conds[[1]]), 
            paste(inputconds$conds[[2]]))
  params <- unlist(strsplit(inputconds$demethod_params[1], ","))
  withProgress(message = 'Running DE Algorithms', detail = inputconds$demethod_params[1], value = 0, {
    initd <- callModule(dprofilerdeanalysis, "deresults", data = data, 
                        columns = cols, conds = conds, params = params)
    if (!is.null(initd$dat()) && nrow(initd$dat()) > 1){
      inputconds$dclist[[1]] <- list(conds = conds, cols = cols, 
                                     init_dedata=initd$dat(),
                                     DEgenes = initd$DEgenes,
                                     init_iterdedata = initd$iterdat(),
                                     IterDEgenes = initd$IterDEgenes,
                                     score = initd$score,
                                     demethod_params = inputconds$demethod_params[1])
    } else {
      return(NULL)
    }
    incProgress(1)
  })
  
  if(length(inputconds$dclist) <1) return(NULL)
  
  return(inputconds$dclist[[1]])
}