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
  
  if(class(data) == "ExpressionSet"){
    metadata <- pData(data)
  }
  
  output$conditionSelector <- renderUI({
    if(class(data) == "ExpressionSet"){
      selected <- selectScRNAConditions(data, session, input)
    } else {
      selected <- selectConditions(data, metadata, session, input)
    }
    selected
  })
  
  observeEvent(input[["conditions_from_meta0"]],{
    if(class(data) == "ExpressionSet"){
      if(!is.null(input[["conditions_from_meta0"]])){
        if(input[["conditions_from_meta0"]]!="No Selection"){
          updateSelectInput(session, "condition", choices = unique(metadata[,input[["conditions_from_meta0"]]]),
                            selected = unique(metadata[,input[["conditions_from_meta0"]]]))
        }
      }
    }
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
                                   actionButtonDE("startDE", "Start", styleclass = "primary")
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
           getConditionSelectorFromMeta(metadata, session = session, input, 1,
                                        1, allsamples, selected1),
           getConditionSelectorFromMeta(metadata, session = session, input, 1,
                                        2, allsamples, selected2)
    ),
    column(12, style='padding-bottom:25px;',
           column(2,  
                  strong(style ='border-bottom: 1px solid #000000; padding-bottom:5px', 
                         "Scoring Parameters:"))
    ),
    column(12, 
           getIterMethodDetails(1, session, input),
    ),
    column(12, style='padding-bottom:25px',
           column(2,  
                  strong(style ='border-bottom: 1px solid #000000; padding-bottom:5px', 
                         "DE Analysis Parameters:"))
           ),
    column(12,
           getSelectInputBox(session$ns("demethod"), 
                             getIconLabel("DE Method", message = "Method for DE Analysis"),
                             1, c("DESeq2", "EdgeR", "Limma"),
                             selectedInput(session$ns("demethod"), 1, "DESeq2", input)),
           getMethodDetails(1, session, input))
  )
  
  if (!is.null(selectedInput(session$ns("conditions_from_meta"), 1, NULL, input)) && 
      selectedInput(session$ns("conditions_from_meta"), 1, NULL, input) != "No Selection"){
    facts <- levels(factor(metadata[,selectedInput(session$ns("conditions_from_meta"),
                                                   1, NULL, input)]))
    facts <- facts[facts != "" & facts != "NA"]
    if (length(facts) != 2) {
      showNotification("There must be exactly 2 groups in the selected condition. 
                         Please use NA or space to remove extra sample groups from metadata selection.", 
                       type = "error")
      updateSelectInput(session, paste0(session$ns("conditions_from_meta"), 1), selected="No Selection" )
    }
  }
  
  return(to_return)
}

#' selectScRNAConditions
#'
#' Selects user input conditions, multiple if present, to be
#' used in DESeq.
#'
#' @param scdata, used single cell dataset 
#' @param metadata, metadatatable to select from metadata
#' @param session, session
#' @param input, input params
#'
#' @examples
#'     x<- selectScRNAConditions()
#'
selectScRNAConditions<-function(scdata = NULL,
                                session = NULL,
                                input = NULL) {
  
  # meta data
  metadata <- pData(scdata)
  character_columns <- colnames(metadata)[!sapply(as.data.frame(metadata),is.numeric)]
  metadata <- metadata[,character_columns]
  
  # meta data selection pane
  to_return <- list(
    column(12,
           getMetaSelector(metadata, session, input),
           getIdentSelectorFromMeta(metadata, session, input, choices = "No Selection"),
           column(3,
                  selectInput(session$ns('deconvolute_samples'), "Samples", choices = colnames(metadata),
                              selected = colnames(metadata)[1])
                  ),
           column(3,
                  selectInput(session$ns("deconvolute_genes"), "DE genes", selected = c("Homogeneous Conditions"),
                              choices = c("Heterogeneous Conditions","Homogeneous Conditions","Marker Genes"))),
           # conditionalPanel(
           #   (condition <- "input.deconvolute_genes == 'Marker Genes'"),
             column(3,
                    textInput(session$ns("top_genes"), label = "Top N Markers", value = "1000"))
           #)
    )
    
  )
  return(to_return)
}

getMetaSelector <- function (metadata = NULL, session = NULL, input = NULL, n = 0) 
{
  if (!is.null(metadata)) {
    df <- metadata
    col_count <- length(colnames(df))
    list(HTML("<hr style=\"color: white; border:solid 1px white;\">"), 
         br(), column(10, 
                      selectInput(paste0(session$ns("conditions_from_meta"), n), 
                                  label = "Select Meta", choices = as.list(c("No Selection", colnames(df)[2:col_count])), 
                                  multiple = FALSE, 
                                  selected = selectedInput(session$ns("conditions_from_meta"), n, "Selection 2", input))))
  }
}

getConditionSelectorFromMeta <- function (metadata = NULL, session = NULL, input = NULL, index = 1, num = 0, 
                                                   choices = NULL, selected = NULL) 
{
  a <- list(column(6, selectInput(paste0(session$ns("condition"), num), 
                                  label = paste0("Condition ", num), choices = choices, 
                                  multiple = TRUE, selected = selected)))
  if (!is.null(metadata)) {
    selected_meta <- selectedInput(session$ns("conditions_from_meta"), 
                                   index, NULL, input)
    if (is.null(selected_meta)) 
      selected_meta <- "No Selection"
    if (selected_meta != "No Selection") {
      old_selection <- ""
      if (!is.null(input[[paste0(session$ns("condition"), num)]])) {
        selected <- input[[paste0(session$ns("condition"), num)]]
      }
      meta_choices_all <- NULL
      if (!is.null(selected_meta)) 
        meta_choices_all <- get_conditions_given_selection(metadata, selected_meta)
      if (old_selection != selected_meta) {
        if (typeof(meta_choices_all) == "character") {
          meta_choices <- list("There must be exactly 2 groups.")
        }
        else {
          meta1 <- meta_choices_all[[2 - (num%%2)]]
          meta_choices <- unlist(meta1, recursive = FALSE)
        }
        selected <- meta_choices
      }
      a <- list(column(6, selectInput(paste0(session$ns("condition"), num), 
                                      label = paste0("Condition ", num), 
                                      choices = choices, multiple = TRUE, selected = selected)))
    }
  }
  return(a)
}

getIdentSelectorFromMeta <- function (metadata = NULL, session = NULL, input = NULL, 
                                      choices = NULL, selected = NULL){
  if (!is.null(metadata)) {
    selected_meta <- selectedInput(session$ns("conditions_from_meta"), 0, NULL, input)
    if (is.null(selected_meta)) 
      selected_meta <- "No Selection"
    if(selected_meta != "No Selection"){
      choices <- unique(metadata[,selected_meta])
    }
    a <- list(column(12, selectInput(session$ns("condition"), 
                                    label = "Identifications", 
                                    choices = choices, multiple = TRUE, selected = choices[1])))
    }
}

#  getMethodDetails
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
getMethodDetails <- function(num = 0, session = NULL, input = NULL) {
  if (num > 0)
    list(
      conditionalPanel(
        (condition <- paste0("input.demethod",num," == 'DESeq2'")),
        ns = session$ns, 
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
        ns = session$ns, 
        getSelectInputBox(session$ns("edgeR_normfact"), 
                          getIconLabel("Normalization", message = "Normalization method used in EdgeR"), 
                          num, c("TMM","RLE","upperquartile","none"), 
                          selectedInput(session$ns("edgeR_normfact"), num, "TMM", input), 2),
        column(2,textInput(paste0(session$ns("dispersion"), num), "Dispersion", 
                           value = isolate(selectedInput(session$ns("dispersion"), 
                                                         num, "0", input) ))),
        getSelectInputBox(session$ns("edgeR_testType"), 
                          getIconLabel("Test Type", message = "Type of model fitting for gene abundances"), 
                          num, c("exactTest", "glmLRT"), 
                          selectedInput(session$ns("edgeR_testType"), num, "exactTest", input))),
      conditionalPanel(
        (condition <- paste0("input.demethod",num," ==  'Limma'")),
        ns = session$ns, 
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
getIterMethodDetails <- function(num = 0, session = NULL, input = NULL) {
  if (num > 0)
    list(
      getSelectInputBox(session$ns("scoremethod"), 
                        getIconLabel("Score Method", message = "Method for estimating Membership scores"), 
                        num, c("Silhouette", "NNLS-based"),
                        selectedInput(session$ns("scoremethod"), num, "Silhouette", input)),
      column(2,textInput(paste0(session$ns("minscore"), num), 
                         getIconLabel("Min. Score", message = "Threshold for acceptable self-membership scores. Score < Min.Score indicates dismemberment of the sample"), 
                         value = isolate(selectedInput(session$ns("minscore"), 
                                                       num, "0.5", input)))),
      # getSelectInputBox(session$ns("iterde_norm"), 
      #                   getIconLabel("Normalization", message = "Normalization method prior to scoring"), 
      #                   num, c("TMM","RLE","upperquartile","none"), 
      #                   selectedInput(session$ns("iterde_norm"), num, "TMM", input), 2),
      getSelectInputBox(session$ns("iterde"),
                        getIconLabel("DE Selection Method", message = "Set of conditions for selecting DE genes on each iteration"), 
                        num, c("Top n Stat.", "Log2FC+Padj"),
                        selectedInput(session$ns("iterde"), num, "Log2FC+Padj", input)),
      conditionalPanel(
        (condition <- paste0("input.iterde",num," == 'Top n Stat.'")),
        ns = session$ns,
        column(2,textInput(paste0(session$ns("topstat"), num), 
                           getIconLabel("Top n Stat.", message = "DE genes with top n DE Analysis statistics. Might be different for each DE method"), 
                           value = isolate(selectedInput(session$ns("topstat"), 
                                                         num, "100", input) )))
      ),
      conditionalPanel(
        (condition <- paste0("input.iterde",num," == 'Log2FC+Padj'")), 
        ns = session$ns,
        column(1,textInput(paste0(session$ns("logfoldchange"), num),
                           getIconLabel("Log2FC", message = "Minimum log fold change"), 
                           value = isolate(selectedInput(session$ns("logfoldchange"),
                                                         num, "0.5", input) ))),
        column(1,textInput(paste0(session$ns("padj"), num), 
                           getIconLabel("P-adj value", message = "Maximum adjusted p-value"),
                           value = isolate(selectedInput(session$ns("padj"),
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
#' @param session session
#' @return data
#' @export
#'
#' @examples
#'     x <- prepDataContainer()
#'
prepDataContainer <- function(data = NULL, input = NULL, session = NULL) {
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
      isolate(input[[paste0("datatype",cnt)]]), sep=",")
  }
  
  # condition inputs for iterative de analysis
  inputconds$demethod_params[cnt] <- paste(inputconds$demethod_params[cnt],
                                           isolate(input[[paste0("scoremethod",cnt)]]),
                                           isolate(input[[paste0("minscore",cnt)]]),
                                           # isolate(input[[paste0("iterde_norm",cnt)]]), 
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
  withProgress(message = 'Running Computational Profiling', 
               detail = paste0("DEmethod: ", params[1], " ScoringMethod: ", params[6]),
               value = 0, {
    initd <- callModule(dprofilerdeanalysis, "deresults", data = data, 
                        columns = cols, conds = conds, params = params, parent_session = session)
    if (!is.null(initd$dat()) && nrow(initd$dat()) > 1){
      inputconds$dclist[[1]] <- list(conds = conds, cols = cols, 
                                     init_dedata=initd$dat(),
                                     DEgenes = initd$DEgenes,
                                     init_iterdedata = initd$iterdat(),
                                     IterDEgenes = initd$IterDEgenes,
                                     score = initd$score,
                                     # deconvolute_genes= initd$deconvolute_genes,
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