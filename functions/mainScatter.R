#' dprofilermainplot
#'
#' Module for a scatter, volcano and ma plots that are going to be used 
#' as a mainplot in dprofiler
#' 
#' @param input, input variables
#' @param output, output objects
#' @param session, session 
#' @param data, a matrix that includes expression values
#' @return main plot
#'
#' @return panel
#' @export
#'
#' @examples
#'     x <- dprofilermainplot()
#'
dprofilermainplot <- function(input = NULL, output = NULL, session = NULL, data = NULL) {
  if (is.null(data$init_dedata)) return(NULL)
  
  # Heterogeneous conditions mainplot
  plotdata_de <-  reactive({
    plotData(data$init_dedata, input)
  })
  
  # Homogeneous conditions mainplot  
  plotdata_iterde <-  reactive({
    plotData(data$init_iterdedata, input)
  })
  
  observe({
    
    # apply filters DE gene main plots 
    iterde_data <- plotdata_iterde()$data
    de_data <- plotdata_de()$data
    iterde_data <- applyFiltersNew(iterde_data[,!colnames(iterde_data) %in% "Legend"], input)
    de_data <- applyFiltersNew(de_data[,!colnames(de_data) %in% "Legend"], input)
    
    # Homogeneous and Heterogeneous conditions main plots
    getMainPlot(input, output, session, "mainde", 6, de_data)
    getMainPlot(input, output, session, "mainiterde", 6, iterde_data)
    
    # apply filter for Initial and Final DEgenes
    iterde_data_foriter <- applyFiltersIter(iterde_data[,!colnames(iterde_data) %in% "Legend"], input)
    de_data_foriter <- applyFiltersIter(de_data[,!colnames(de_data) %in% "Legend"], input)
    
    # Initial, Overlapping, and Final Genes main plots 
    getMainPlot(input, output, session, "maininitial", 4, de_data_foriter, 
                which_genes = "Initial", DEgenes = data$DEgenes, IterDEgenes = data$IterDEgenes)
    getMainPlot(input, output, session, "mainoverlap", 4, iterde_data_foriter, 
                which_genes = "Overlap", DEgenes = data$DEgenes, IterDEgenes = data$IterDEgenes)
    getMainPlot(input, output, session, "mainfinal", 4, iterde_data_foriter, 
                which_genes = "Final", DEgenes = data$DEgenes, IterDEgenes = data$IterDEgenes)
    
  })
  
  # Control Panel 
  output$mainPlotControlsUI <- renderUI({
    if (input$mainplot == "scatter"){
      x <- paste0('log10 Norm. Mean(Read Counts) in cond1')
      y <- paste0('log10 Norm. Mean(Read Counts) in cond2')
    }else if  (input$mainplot == "volcano"){
      x <- "log2FC"
      y <- "-log10padj"
    }else {
      x <- "A"
      y <- "M"
    }
    list(
      textInput(session$ns('xlab'),'x label', x),
      textInput(session$ns('ylab'),'y label', y),
      checkboxInput(session$ns('labelsearched'), 'Label searched points', value = FALSE),
      conditionalPanel(paste0("input['",session$ns("labelsearched"), "']"),
                       colourpicker::colourInput(session$ns("labelcolor"), "Label colour", "black"),
                       selectInput(session$ns("labelsize"), "Label Size", choices=c(6:30), selected=14))
    )
  })
  
  selectedPoint <- reactive({
    eventdata <- event_data("plotly_click", source = session$ns("source"))
    if (is.null(eventdata)){
      eventdata <- event_data("plotly_hover", source = session$ns("source"))
    }
    key <- ""
    if (!is.null(eventdata$key))
      key <- as.vector(unlist(eventdata$key))
    return(key)
  })
  
  getSelected  <- reactive({
    keys <- NULL
    selGeneList <- event_data("plotly_selected", source = session$ns("source"))
    if (is.null(selGeneList$key)) return (NULL)
    keys <- as.vector(unlist(selGeneList$key))
    return(keys)
  })
  
  list(shg=(selectedPoint), shgClicked=(selectedPoint), selGenes=(getSelected))
}


#' getMainPlot
#' 
#' a wrapper function for all main plots in Dprofiler
#'
#' @param input input
#' @param output output
#' @param session session
#' @param mainname main plot name
#' @param width shiny box width
#' @param plotdata plot data
#' @param which_genes if it is true, DEgenes and IterDEgenes will be used for subseting
#' @param DEgenes list of initial DE genes
#' @param IterDEgenes list of Final DE genes
#'
#' @return
#' @export
#'
#' @examples
getMainPlot <- function(input = NULL, output = NULL, session = NULL, mainname = NULL, 
                        width = NULL, plotdata = NULL,
                        which_genes = NULL, DEgenes = NULL, IterDEgenes = NULL){
 
  mainnameplot <- paste0(mainname, "plot")
  output[[mainnameplot]] <- renderUI({
    list(
      column(width,
             shinydashboard::box(
               collapsible = TRUE, 
               title = paste(which_genes, "Genes Main Plots", seo = " "), 
               status = "primary", 
               solidHeader = TRUE, width = NULL, draggable = TRUE, 
               column(12,
                plotlyOutput(session$ns(mainname))
               )
             ))
    )
  })
  output[[mainname]] <- renderPlotly({
    if(!is.null(which_genes)){
      if(which_genes == "Initial"){
        genes <- setdiff(DEgenes, IterDEgenes)
        if(length(genes) == 0) genes <- DEgenes
      } else if(which_genes == "Overlap"){
        genes <- intersect(DEgenes, IterDEgenes)
        if(length(genes) == 0) genes <- IterDEgenes
      } else {
        genes <- setdiff(IterDEgenes, DEgenes)
        if(length(genes) == 0) genes <- IterDEgenes
      }
      data <- plotdata[genes,]
    } else {
      data <- plotdata
    }
    mainScatter(input, data, session$ns("source"))
  })
}

#' mainScatter
#'
#' Creates the main scatter, volcano or MA plot to be displayed within the main
#' panel. A version of debrowser's mainScatterNew function with automated width and height
#' 
#' @param input, input params
#' @param data, dataframe that has log2FoldChange and log10padj values
#' @param source, for event triggering to select genes
#' @return scatter, volcano or MA plot
#'
#' @examples
#'     
#'     x <- mainScatter()
#'
#' @export
#'
mainScatter <- function(input = NULL, data = NULL, source = NULL) {
  if ( is.null(data) ) return(NULL)
  
  data <- na.omit(data)
  p <-   plot_ly(source = source, data=data, x=~x, y=~y, key=~key, alpha = 0.8,
                 color=~droplevels(factor(Legend, levels = c("NS","Up","Down"))), 
                 colors=getLegendColors(getLevelOrder(unique(data$Legend))), 
                 type="scatter", mode = "markers",
                 text=~paste("<b>", ID, "</b><br>",
                             "<br>", "padj=", format.pval(padj, digits = 2), " ",
                             "-log10padj=", round(log10padj, digits = 2),
                             "<br>", "log2FC=", round(log2FoldChange, digits = 2), " ",
                             "foldChange=", round(foldChange, digits = 2),
                             "<br>", sep = " ")) %>%
    plotly::layout(xaxis = list(title = input$xlab),
                   yaxis = list(title = input$ylab), 
                   autosize = TRUE)
  
  if (!is.null(input$labelsearched) && input$labelsearched == TRUE){
    searched_genes <- data[(data$Legend == "GS"),]
    a <- list()
    for (i in seq_len(nrow(searched_genes))) {
      m <- searched_genes[i, ]
      a[[i]] <- list(
        x = m$x,
        y = m$y,
        text = rownames(m),
        color = 'blue',
        xref = "x",
        yref = "y",
        showarrow = TRUE,
        arrowhead = 0.5,
        ax = 20,
        ay = -40,
        font = list(color = input$labelcolor,
                    face = 2,
                    size = input$labelsize)
      )
    }
    
    p <- p %>%  plotly::layout(annotations = a)
  }
  if (!is.null(input$svg) && input$svg == TRUE)
    p <- p %>% config(toImageButtonOptions = list(format = "svg"))
  p$elementId <- NULL
  
  return(p)
}

#' addDataCols
#' 
#' add data columns to de results
#'
#' @param data data 
#' @param de_res DE results
#' @param cols columns
#' @param conds conditions
#'
#' @return
#' @export
#'
#' @examples
addDataCols <- function (data = NULL, de_res = NULL, cols = NULL, conds = NULL) 
{
  if (is.null(data) || (nrow(de_res) == 0 && ncol(de_res) == 
                        0)) 
    return(NULL)
  norm_data <- data[, cols]
  coldata <- prepGroup(conds, cols)
  mean_cond_first <- getMean(norm_data, as.vector(coldata[coldata$group == 
                                                            levels(coldata$group)[1], "libname"]))
  mean_cond_second <- getMean(norm_data, as.vector(coldata[coldata$group == 
                                                             levels(coldata$group)[2], "libname"]))
  m <- cbind(rownames(de_res), norm_data[de_res$gene, cols], 
             log10(unlist(mean_cond_second) + 1), log10(unlist(mean_cond_first) + 1), 
             de_res[rownames(de_res), c("padj", "log2FoldChange", "pvalue", "stat")], 
             2^de_res[rownames(de_res), "log2FoldChange"], -1 * log10(de_res[rownames(de_res), "padj"]))
  colnames(m) <- c("ID", cols, "y", "x", 
                   "padj", "log2FoldChange", "pvalue", 
                   "stat", "foldChange", "log10padj")
  m <- as.data.frame(m)
  m$padj[is.na(m[paste0("padj")])] <- 1
  m$pvalue[is.na(m[paste0("pvalue")])] <- 1
  
  if(!is.null(de_res$Comparison)){
    m <- data.frame(m, Comparison = de_res[rownames(de_res), "Comparison"])
  }
  
  m
}

mainPlotControlsUI <- function (id)
{
  ns <- NS(id)
  list(shinydashboard::menuItem(" Plot Type", startExpanded = TRUE,
                                radioButtons(ns("mainplot"), "Main Plots:",
                                             c(Scatter = "scatter", VolcanoPlot = "volcano",
                                               MAPlot = "maplot"))),
       shinydashboard::menuItem("Main Options", startExpanded = TRUE,
                                sliderInput(ns("backperc"), "Background Data(%):", min = 10, max = 100, value = 100, sep = "", animate = FALSE),
                                conditionalPanel(condition <- paste0("input['", ns("mainplot"), "'] == 'volcano'"),
                                                 sliderInput(ns("log10padjCutoff"),"Log10 padj value cutoff:",
                                                             min = 2, max = 100, value = 60, sep = "", animate = FALSE)),
                                uiOutput(ns("mainPlotControlsUI"))))
}


applyFiltersIter <- function (data = NULL, input = NULL) 
{
  if (is.null(data)) 
    return(NULL)
  m <- data
  m$Legend[m$log2FoldChange > 0 ] <- "Up"
  m$Legend[m$log2FoldChange < 0 ] <- "Down"
  return(m)
}

getLegendColors <-  function(Legend = c("Up", "Down", "NS")) 
{
  colors <- c()
  for (i in seq(1:length(Legend))) {
    if (Legend[i] == "Up") {
      colors <- c(colors, "#ff0000")
    }
    else if (Legend[i] == "Down") {
      colors <- c(colors, "#0000ff")
    }
    else if (Legend[i] == "NS") {
      colors <- c(colors, "#808080")
    }
    else if (Legend[i] == "GS") {
      colors <- c(colors, "#008000")
    }
  }
  colors
}