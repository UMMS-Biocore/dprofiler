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
           getSelectInputBox("demethod", 
                             getIconLabel("DE Method", message = "Method for DE Analysis"),
                             1, c("DESeq2", "EdgeR", "Limma"),
                             selectedInput("demethod", 1, "DESeq2", input)),
           getMethodDetails(1, input))
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

#getMethodDetails
#'
#' get the detail boxes after DE method selected 
#'
#' @param num, panel that is going to be shown
#' @param input, user input
#' @examples
#'     x <- getMethodDetails()
#'
#' @export
#'
#'
getMethodDetails <- function(num = 0, input = NULL) {
  if (num > 0)
    list(
      conditionalPanel(
        (condition <- paste0("input.demethod",num," == 'DESeq2'")),
        getSelectInputBox("fitType", 
                          getIconLabel("Fit Type", message = "fitting type for dispersion estimate"), 
                          num, c("parametric", "local", "mean"), 
                          selectedInput("testType", num, "parametric",
                                        input), 2),
        getSelectInputBox("betaPrior", 
                          getIconLabel("Beta Prior", message = "Use a zero-mean normal prior. Default for Wald statistics"), 
                          num, c(FALSE, TRUE), 
                          selectedInput("betaPrior", num,
                                        FALSE, input),2),
        getSelectInputBox("testType", 
                          getIconLabel("Test Type", message = "Test statistics used in DESeq2"), 
                          num, c("LRT", "Wald"),  
                          selectedInput("testType", num, "LRT", input)),
        getSelectInputBox("shrinkage", 
                          getIconLabel("Shrinkage", message = "The method of shrinkage estimate for log2FC"), 
                          num, c("None", "apeglm", "ashr", "normal"),
                          selectedInput("shrinkage", num, "None", input))),
      conditionalPanel(
        (condition <- paste0("input.demethod",num," == 'EdgeR'")),
        getSelectInputBox("edgeR_normfact", 
                          getIconLabel("Normalization", message = "Normalization method used in EdgeR"), 
                          num, c("TMM","RLE","upperquartile","none"), 
                          selectedInput("edgeR_normfact", num, "TMM", input), 2),
        column(2,textInput(paste0("dispersion", num), "Dispersion", 
                           value = isolate(selectedInput("dispersion", 
                                                         num, "0", input) ))),
        getSelectInputBox("edgeR_testType", 
                          getIconLabel("Test Type", message = "Type of model fitting for gene abundances"), 
                          num, c("exactTest", "glmLRT"), 
                          selectedInput("edgeR_testType", num,
                                        "exactTest", input))),
      conditionalPanel(
        (condition <- paste0("input.demethod",num," ==  'Limma'")),
        getSelectInputBox("limma_normfact", getIconLabel("Normalization", message = "Normalization method used in Limma"), 
                          num, c("TMM","RLE","upperquartile","none"), 
                          selectedInput("limma_normfact", num, "TMM", input), 2),
        getSelectInputBox("limma_fitType", 
                          getIconLabel("Fit Type", message = "Type of model fitting: 'ls' for least squares, 'robust' for robust regression"), 
                          num, c("ls", "robust"), 
                          selectedInput("limma_fitType", num, "ls", input)),
        getSelectInputBox("normBetween", 
                          getIconLabel("Norm. Bet. Arrays", message = "Normalization Between Array"), 
                          num, c("none", "scale", "quantile", "cyclicloess",
                            "Aquantile", "Gquantile", "Rquantile","Tquantile"),
                          selectedInput("normBetween", num, "none", input))),
      br())
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
      getSelectInputBox("scoremethod", 
                        getIconLabel("Score Method", message = "Method for estimating Membership scores"), 
                        num, c("Silhouette", "NNLS-based"),
                        selectedInput("scoremethod", num, "Silhouette", input)),
      column(2,textInput(paste0("minscore", num), 
                         getIconLabel("Min. Score", message = "Threshold for acceptable self-membership scores. Score < Min.Score indicates dismemberment of the sample"), 
                         value = isolate(selectedInput("minscore", 
                                                       num, "0.5", input)))),
      getSelectInputBox("iterde_norm", 
                        getIconLabel("Normalization", message = "Normalization method prior to scoring"), 
                        num, c("TMM","RLE","upperquartile","none"), 
                        selectedInput("iterde_norm", num, "TMM", input), 2),
      getSelectInputBox("iterde",
                        getIconLabel("DE Selection Method", message = "Set of conditions for selecting DE genes on each iteration"), 
                        num, c("Top n Stat.", "Log2FC+Padj"),
                        selectedInput("iterde", num, "Log2FC+Padj", input)),
      conditionalPanel(
        (condition <- paste0("input.iterde",num," == 'Top n Stat.'")),
        column(2,textInput(paste0("topstat", num), 
                           getIconLabel("Top n Stat.", message = "DE genes with top n DE Analysis statistics. Might be different for each DE method"), 
                           value = isolate(selectedInput("topstat", 
                                                         num, "100", input) )))
      ),
      conditionalPanel(
        (condition <- paste0("input.iterde",num," == 'Log2FC+Padj'")),
        column(1,textInput(paste0("logfoldchange", num),
                           getIconLabel("Log2FC", message = "Minimum log fold change"), 
                           value = isolate(selectedInput("logfoldchange",
                                                         num, "0.5", input) ))),
        column(1,textInput(paste0("padj", num), 
                           getIconLabel("P-adj value", message = "Maximum adjusted p-value"),
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
  withProgress(message = 'Running Heterogeneity Detection', detail = inputconds$demethod_params[1], value = 0, {
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


#' getIconLabel
#' 
#' creates a label with information icon
#'
#' @param label label
#' @param message message
#'
getIconLabel <- function(label = NULL, message = NULL){
  
  icon_label <- tags$span(label, 
    tags$i(
      class = "glyphicon glyphicon-info-sign", 
      style = "color:#0072B2;",
      title = message
    )
  )
  return(icon_label)
}