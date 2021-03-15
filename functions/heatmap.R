#' dprofilerheatmap
#'
#' Heatmap module to create interactive heatmaps and get selected list from
#' a heatmap
#' @param input, input variables
#' @param output, output objects
#' @param session, session 
#' @param expdata, a matrix that includes expression values
#' @return heatmapply plot
#'
#' @examples
#'     x <- dprofilerheatmap()
#'
#' @export
#'
#'
dprofilerheatmap <- function( input, output, session, expdata = NULL){
  if(is.null(expdata)) return(NULL)
  output$heatmap <- renderPlotly({
    shinyjs::onevent("mousemove", "heatmap", js$getHoverName(session$ns("hoveredgenename")))
    shinyjs::onevent("click", "heatmap", js$getHoverName(session$ns("hoveredgenenameclick")))
    
    withProgress(message = 'Drawing Heatmap', detail = "interactive", value = 0, {
      runHeatmap(input, session, orderData())
    })
  })
  output$heatmap2 <- renderPlot({
    withProgress(message = 'Drawing Heatmap', detail = "non-interactive", value = 0, {
      runHeatmap2(input, session, orderData())
    })
  })
  heatdata <- reactive({
    cld <- prepHeatData(expdata, input)
    if (input$kmeansControl)
    {
      res <- niceKmeans(cld, input)
      cld <- res$clustered
    }
    cld
  })
  
  button <- reactiveVal(FALSE)
  orderData <- reactive({
    newclus <- heatdata()
    if (input$changeOrder && isolate(button()) && !is.null(input$clusterorder)){
      newclus <- changeClusterOrder(isolate(input$clusterorder), newclus)
    }
    button(FALSE)
    newclus
  })
  observeEvent(input$changeOrder,{
    button(TRUE)
  })
  output$heatmapUI <- renderUI({
    if (is.null(input$interactive)) return(NULL)
    column(4,
      shinydashboard::box(
        collapsible = TRUE, title = session$ns("Heatmap"), status = "primary", 
        solidHeader = TRUE, width = NULL,
        draggable = TRUE, getPlotArea(input, session)
      )
    )
  })
  
  hselGenes <- reactive({
    if (is.null(input$selgenenames)) return("")
    unlist(strsplit(input$selgenenames, split=","))
  })
  shg <- reactive({
    if (is.null(input$hoveredgenename)) return("")
    js$getSelectedGenes(session$ns("heatmap"), session$ns("selgenenames"))
    input$hoveredgenename
  })
  observe({
    if(!input$changeOrder)
      updateTextInput(session, "clusterorder", value = paste(seq(1:input$knum), collapse=","))

    if (is.null(shg()))
      js$getSelectedGenes()
  })
  shgClicked <- reactive({
    if (is.null(input$hoveredgenenameclick) || input$hoveredgenenameclick == "") 
      return(input$hoveredgenename)
    input$hoveredgenenameclick
  })
  
  list(shg = (shg), shgClicked=(shgClicked), selGenes=(hselGenes), getSelected = (orderData))
}


getPlotArea <- function (input = NULL, session = NULL) 
{
  if (is.null(input)) 
    return(NULL)
  ret <- c()
  if (input$interactive) {
    ret <- plotlyOutput(session$ns("heatmap"))
  }
  else {
    ret <- plotOutput(session$ns("heatmap2"))
  }
  ret
}