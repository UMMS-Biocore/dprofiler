dprofilerboxmainplot <- function (input = NULL, output = NULL, session = NULL, data = NULL, 
                                  cols = NULL, conds = NULL, key = NULL) 
{
  if (is.null(data)) 
    return(NULL)
  output$BoxMain <- renderPlotly({
    getBoxMainPlot(data, cols, conds, key, title = "", 
                   input)
  })
  output$BoxMainUI <- renderUI({
    column(4,
      shinydashboard::box(collapsible = TRUE, title = "Gene Box Plots", 
                          status = "primary", solidHeader = TRUE, width = NULL, 
                          draggable = TRUE, plotlyOutput(session$ns("BoxMain")))
    )
  })
}

dprofilerbarmainplot <- function (input, output, session, data = NULL, cols = NULL, conds = NULL, 
          key = NULL) 
{
  if (is.null(data)) 
    return(NULL)
  output$BarMainUI <- renderUI({
    column(4,
      shinydashboard::box(collapsible = TRUE, title = "Gene Bar Plot", 
                          status = "primary", solidHeader = TRUE, width = NULL, 
                          draggable = TRUE, plotlyOutput(session$ns("BarMain")))
    )
  })
  output$BarMain <- renderPlotly({
    getBarMainPlot(data, cols, conds, key, title = "", 
                   input = input)
  })
}