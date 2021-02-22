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
debrowsercondselectIter <- function(input = NULL, output = NULL, session = NULL, data = NULL, metadata = NULL) {
  if (is.null(data)) return(NULL)
  choicecounter <- reactiveVal(0)
  output$conditionSelectorIter <- renderUI({
    selectConditionsIter(data, metadata, choicecounter(), session, input)
  })
  
  observeEvent(input$add_btn_iter, {
    choicecounter(choicecounter() + 1)
  })
  observeEvent(input$rm_btn_iter, {
    if (choicecounter() > 0) 
      choicecounter(choicecounter() - 1)
  })
  list(cc = choicecounter)
}


#' condSelectIterUI
#' Creates a panel to select samples for each condition in the 
#' Iterative DE Analysis Section
#'
#' @return panel
#' @examples
#'     x <- condSelectIterUI()
#'
#' @export
#'
condSelectIterUI<- function () {
  list(
    shinydashboard::box(title = "Comparison Selection",
                        solidHeader = TRUE, status = "info",  width = NULL, height = NULL, collapsible = TRUE,
                        fluidRow(
                          uiOutput("conditionSelectorIter"),
                          column(12,actionButtonDE("add_btn_iter", "Add New Comparison",styleclass = "primary"),
                                 actionButtonDE("rm_btn_iter", "Remove", styleclass = "primary"),
                                 getHelpButton("method", "http://debrowser.readthedocs.io/en/master/deseq/deseq.html"),
                                 # conditionalPanel(condition = ("output.condReady>0"),
                                 #                  actionButtonDE("startDE", "Start DE", styleclass = "primary"))
                                 # actionButtonDE("startDE", "Start DE", styleclass = "primary")
                          )
                        ))
  )
}

#' selectConditionsIter
#'
#' Selects user input conditions, multiple if present, to be
#' used in DESeq for Iterative DE Analysis
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
#'     x<- selectConditionsIter()
#'
#' @export
#'
selectConditionsIter <-function(Dataset = NULL,
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
             getSelectInputBox("iterativemethod", "Iterative Method", i, 
                               c("NNLS", "Silhouette"),
                               "Silhouette")
      )
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