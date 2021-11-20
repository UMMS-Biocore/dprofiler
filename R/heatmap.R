#' dprofilerheatmap
#'
#' Heatmap module to create interactive heatmaps and get selected list from
#' a heatmap
#' 
#' @param input, input variables
#' @param output, output objects
#' @param session, session 
#' @param expdata, a matrix that includes expression values
#'
#' @examples
#'     x <- dprofilerheatmap()
#'     
#' @export
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


#' getPlotArea
#'
#' a version of debrowser::getPlotArea function with automatic width and length
#' 
#' @param input input variables
#' @param session session
#'
#' @examples
#'      x <- getPlotArea()
#'     
#' @export
#' 
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

#' heatmapJScode
#'
#' heatmap JS code for selection functionality. Adapted from debrowser::heatmapJScode()
#'
#' @return JS Code
#' 
#' @examples
#'     x <- heatmapJScode()
#'     
#' @export
#' 
heatmapJScode <- function() {        
  'shinyjs.getHoverName = function(params){
    
    var defaultParams = {
    controlname : "hoveredgenename"
    };
    params = shinyjs.getParams(params, defaultParams);
    var out = ""
    
    if (typeof  document.getElementsByClassName("nums")[0] != "undefined"){
    if (typeof  document.getElementsByClassName("nums")[0].querySelectorAll("tspan.line")[0] != "undefined"){
    out = document.getElementsByClassName("nums")[0].querySelectorAll("tspan.line")[0].innerHTML.match("row: (.*)")[1]
    $("#deresults-heatmap").attr("gname", out)
    }
    }
    Shiny.onInputChange(params.controlname, $("#deresults-heatmap").attr("gname"));
    }
    shinyjs.resetInputParam = function(params){
        var defaultParams = {
                controlname : "hoveredgenename"
        };
        params = shinyjs.getParams(params, defaultParams);
        console.log(params.controlname)
        Shiny.onInputChange(params.controlname, "");
    }
    shinyjs.getSelectedGenes = function(params){
    var defaultParams = {
    plotId : "heatmap",
    controlname : "selgenenames"
    };
    params = shinyjs.getParams(params, defaultParams);
    var count = document.getElementById(params.plotId).querySelectorAll("g.y2tick").length
    var start = 0
    var out = ""
    
    for (i = start; i < count; i++)
    {
        if (typeof document.getElementById(params.plotId).querySelectorAll("g.y2tick")[i] != "undefined"){
        out += document.getElementById(params.plotId).querySelectorAll("g.y2tick")[i].innerHTML.match(">(.*)</text>")[1]  + ","
        }
    }
    Shiny.onInputChange(params.controlname, out);
    }'
}

#' getJSLine
#'
#' heatmap JS code for selection functionality. Adapted from debrowser::getJSLine()
#'
#' @return JS Code
#' 
#' @examples
#'     x <- getJSLine()
#'     
#' @export
#' 
getJSLine <- function() 
{
  list(shinyjs::useShinyjs(), shinyjs::extendShinyjs(text = heatmapJScode(), 
                                                     functions = c("getHoverName", "getSelectedGenes",
                                                                   "resetInputParam")))
}