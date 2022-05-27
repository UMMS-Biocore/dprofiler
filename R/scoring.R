#' getScoreDetails
#'
#' get score details of Iterative DE analysis
#' 
#' @param output output
#' @param session session
#' @param plotname the name of the plot
#' @param DEscores Cross Condition Scores of impure conditions
#' @param IterDEscores Cross Condition Scores of pure conditions
#' @param OutlierDEscores Intra Condition Scores of impure conditions
#' @param OutlierIterDEscores Intra Condition Scores of pure conditions
#'
#' @examples
#'      x <- getScoreDetails()
#'     
#' @export
#' 
getScoreDetails <- function(output = NULL, session = NULL, 
                                 plotname = NULL, DEscores = NULL, IterDEscores = NULL, OutlierDEscores = NULL, OutlierIterDEscores = NULL) {
  if (is.null(DEscores)) return(NULL)
  
  output[[plotname]] <- renderPlotly({
    dat <- rbind(
      data.frame(IterDEscores, beforeafter = "Before Profiling", type = "Cross Cond."),
      data.frame(DEscores, beforeafter = "After Profiling", type = "Cross Cond."),
      data.frame(OutlierIterDEscores, beforeafter = "Before Profiling", type = "Intra Cond."),
      data.frame(OutlierDEscores, beforeafter = "After Profiling", type = "Intra Cond.")
    )
    dat$Scores <- as.numeric(dat$Scores)
    p <- ggplot(data=dat, aes(x = reorder(Samples,Scores), y = Scores)) +
      geom_bar(aes(fill = Conds), stat="identity") + 
      facet_grid(type ~ beforeafter) + 
      scale_y_continuous(limits=c(0, 1)) + 
      xlab("Samples") +
      theme(axis.text.x = element_text(angle = 45),
            axis.title.x=element_blank(),
            axis.title.y=element_blank()) 

    p <- ggplotly(p)
    
    # Manage plotly labels
    for(i in 1:length(p$x$data)){
      temp <- p$x$data[[i]]$text
      temp <- gsub("Conds: ", "", temp)
      temp <- gsub("reorder\\(Samples, Scores\\):", "Sample:", temp)
      temp <- gsub("Scores:", "Score:", temp)
      p$x$data[[i]]$text <- temp
    }
    
    p
  })
}

#' getIterDESummary
#' 
#' get summary of Iterative DE Analysis results and venn diagram
#'
#' @param output output 
#' @param session session
#' @param vennname venn diagram name
#' @param summaryname summary table name
#' @param deres DE results
#' @param params DE parameters
#'
#' @examples
#'      x <- getIterDESummary()
#'     
#' @export
#'   
getIterDESummary <- function(output = NULL, session = NULL, vennname = NULL, summaryname = NULL, 
                             deres = NULL, params = NULL){
  if (is.null(output)) return(NULL)
  
  output[[summaryname]] <- renderUI({
   style = "padding-right: 10px"
   texts <- tags$div(
      h4("Summary"),
      tags$table(
        tags$tr(
          tags$td(style = style, p(strong("# of Iterations:"), deres$NumberofIters)),
          tags$td(style = style, p(strong("# of Initial DE genes:"), length(deres$DEgenes)))
        ),
        tags$tr(
          tags$td(style = style, p(strong("# of Hetero. Samples:"), length(deres$cleaned_columns))),
          tags$td(style = style, p(strong("# of Final DE genes:"), length(deres$IterDEgenes)))
        )
      )
   )
  })
  
  # Scoring parameters 
  scoresummaryname <- paste0(summaryname,"Iter")
  output[[scoresummaryname]] <- renderUI({
    
    style = "padding-right: 10px"
    if(params[9] == "Stat."){
      additional_texts <- paste(p(strong("# of Top Statistics:"), params[11]), sep = " ")
    } else {
      additional_texts <- paste(p(strong("Log2FC:"), params[9], strong("P-adj:"),  params[10]), sep = " ")
    }
    texts <- tags$div(
      h4("Scoring Parameters"),
      tags$table(
        tags$tr(
          tags$td(style = style, p(strong("Score Method:"), params[6])),
          tags$td(style = style, p(strong("Selection Method:"),  params[8])),
        ),
        tags$tr(
          tags$td(style = style, p(strong("Min. Score:"),  params[7])),
          tags$td(style = style, HTML(additional_texts))
        )
      )
    )
    texts
  })
  
  # DE parameters
  desummaryname <- paste0(summaryname,"DE")
  
  output[[desummaryname]] <- renderUI({
    
    style = "padding-right: 10px"
    if(params[1] == "DESeq2"){
      param_text <- c("DE method:", "Fit Type:", "Beta Prior:", "Test Type:", "Shrinkage:")
    } else if(params[1] == "EdgeR") {
      param_text <- c("DE method:", "Normalization:", "Dispersion:", "Test Type:", "")
    } else{
      param_text <- c("DE method:", "Normalization:", "Fit Type:", "Norm.Bet.Arrays:", "")
    }

    texts <- tags$div(
      h4("DE parameters"),
      tags$table(
        tags$tr(
          tags$td(style = style, p(strong(param_text[1]), params[1])),
          tags$td(style = style, p(strong(param_text[3]),  params[3])),
        ),
        tags$tr(
          tags$td(style = style, p(strong(param_text[2]),  params[2])),
          tags$td(style = style, p(strong(param_text[4]),  params[4])),
        ),
        if(params[1] == "DESeq2"){
          tags$tr(
            tags$td(style = style, p(strong("Normalization:"),  params[5]))
          )
        }
      )
    )
    texts
  })
  
  x <- list(
    A = deres$DEgenes, 
    B = deres$IterDEgenes
  )
  output[[vennname]] <- renderPlot({
    display_venn(x, category.names = c("Init. DE genes" , "Final DE genes"),
                 fill = c("#999999", "#E69F00"))
  })
}

#' display_venn
#' 
#' display venn diagram
#'
#' @param x a list of elements for each groups
#' @param category.names names of venn diagram categories
#' @param ... addtional parameters passed to draw.pairwise.venn function. 
#'
#' @examples
#'      x <- display_venn()
#'     
#' @export
#'    
display_venn <- function(x = NULL, category.names = NULL, ...){
  if (is.null(x)) return(NULL)
  
  grid.newpage()
  venn_object <- draw.pairwise.venn(area1 = length(x[[1]]), 
                                    area2 = length(x[[2]]), 
                                    cross.area = length(intersect(x[[1]], x[[2]])), 
                                    category = category.names, 
                                    ind = FALSE, ...)
  grid.draw(venn_object)
}