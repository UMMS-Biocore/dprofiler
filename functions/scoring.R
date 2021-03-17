#' getScoreDetails
#'
#' get score details of columns
#' 
#' @param output output
#' @param session session
#' @param plotname the name of the plot
#' @param DEscores Scores of impure conditions
#' @param IterDEscores Scores of pure conditions
#'
#' @return
#' @export
#'
#' @examples
getScoreDetails <- function(output = NULL, session = NULL, 
                                 plotname = NULL, DEscores = NULL, IterDEscores = NULL) {
  output[[plotname]] <- renderPlotly({
    dat <- rbind(
      data.frame(IterDEscores, type = "Scores of Pure Conditions"),
      data.frame(DEscores, type = "Scores of Impure Conditions")
    )
    dat$Scores <- as.numeric(dat$Scores)
    p <- ggplot(data=dat, aes(x = reorder(Samples,Scores), y = Scores)) +
      geom_bar(aes(fill = Conds), stat="identity") + 
      facet_grid(. ~ type) + 
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
#' @return
#' @export
#'
#' @examples
getIterDESummary <- function(output = NULL, session = NULL, vennname = NULL, summaryname = NULL, 
                             deres = NULL, params = NULL){
  
  # Iteration Summary
  # output[[summaryname]] <- renderTable({
  #   result <- rbind(deres$NumberofIters, length(deres$cleaned_columns), 
  #                   length(deres$DEgenes), length(deres$IterDEgenes))
  #   rownames(result) <- c("# of Iterations", "# of Hetero. Samples", 
  #                         "# of Initial DE genes", "# of Final DE genes ")
  #   colnames(result) <- "Value"
  #   result
  # }, digits = 0, rownames = TRUE, align = "lc")
  output[[summaryname]] <- renderUI({
    texts <- paste(
      h4("Iter. DE Summary"),
      p(strong("# of Iterations:"), deres$NumberofIters),
      p(strong("# of Hetero. Samples:"), length(deres$cleaned_columns)), 
      p(strong("# of Initial DE genes:"), length(deres$DEgenes)),
      p(strong("# of Final DE genes:"), length(deres$IterDEgenes)), 
      sep = " "
      )
    HTML(texts)
  })
  
  # Scoring parameters 
  scoresummaryname <- paste0(summaryname,"Iter")
  # output[[scoresummaryname]] <- renderTable({
  #   result <- params[6:9]
  #   if(params[9] == "Stat."){
  #     result <- matrix(c(result, params[10:11]),ncol = 1)
  #     rownames(result) <- c("Score Method", "Min. Score", "Normalization","Selection Method", 
  #                           "Log2FC", "P-adj")
  #   } else {
  #     result <- matrix(c(result, params[12]),ncol = 1)
  #     rownames(result) <- c("Score Method", "Min. Score", "Normalization","Selection Method", 
  #                           "# of Top Statistics")
  #   }
  #   colnames(result) <- "Value"
  #   result
  # }, digits = 0, rownames = TRUE, align = "lc")
  output[[scoresummaryname]] <- renderUI({
    texts <- paste(
      h4("Scoring Parameters"),
      p(strong("Score Method:"), params[6]),
      p(strong("Min. Score:"),  params[7]),
      p(strong("Normalization:"),  params[8]),
      p(strong("Selection Method:"),  params[9]),
      sep = " "
    )
    if(params[9] == "Stat."){
      texts <- paste(texts,
                     p(strong("# of Top Statistics:"), params[12]),
                     sep = " "
      )
    } else {
      texts <- paste(texts,
                     p(strong("Log2FC:"), params[10], strong("P-adj:"),  params[11]),
                     sep = " "
      )
    }
    HTML(texts)
  })
  
  # DE parameters
  desummaryname <- paste0(summaryname,"DE")
  # output[[desummaryname]] <- renderTable({
  #   if(params[1] == "DESeq2"){
  #     result <- matrix(params[1:5],ncol = 1)
  #     rownames(result) <- c("DE Method", "Fit Type", "Beta Prior","Test Type", "Shrinkage")
  #   } else if(params[1] == "EdgeR") {
  #     result <- matrix(params[1:4],ncol = 1)
  #     rownames(result) <- c("DE Method", "Normalization", "Dispersion", "Test Type")
  #   } else{
  #     result <- matrix(params[1:4],ncol = 1)
  #     rownames(result) <- c("DE Method", "Normalization", "Fit Type", "Norm.Bet.Arrays")
  #   }
  #   colnames(result) <- "Value"
  #   result
  # }, digits = 0, rownames = TRUE, align = "lc")
  
  output[[desummaryname]] <- renderUI({
    texts <- paste(h4("DE Summary"),
                   p(strong("DE Method:"), params[1]), 
                   sep = " ")
    
    if(params[1] == "DESeq2"){
      texts <- paste(texts,
        p(strong("Fit Type:"), params[2]),
        p(strong("Beta Prior:"), params[3]), 
        p(strong("Test Type:"), params[4]),
        p(strong("Shrinkage:"), params[5]),
        sep = " "
      )
    } else if(params[1] == "EdgeR") {
      texts <- paste(texts,
                     p(strong("Normalization:"), params[2]),
                     p(strong("Dispersion:"), params[3]), 
                     p(strong("Test Type:"), params[4]),
                     sep = " "
      )
    } else{
      texts <- paste(texts,
                     p(strong("Normalization:"), params[2]),
                     p(strong("Fit Type:"), params[3]), 
                     p(strong("Norm.Bet.Arrays:"), params[4]),
                     sep = " "
      )
    }
    HTML(texts)
  })
  
  x <- list(
    A = deres$DEgenes, 
    B = deres$IterDEgenes
  )
  output[[vennname]] <- renderPlot({
    display_venn(x,
                 category.names = c("Init. DE genes" , "Final DE genes"),
                 fill = c("#999999", "#E69F00"))
  })
}

#' display_venn
#' 
#' display venn diagram
#'
#' @param x a list of elements for each groups
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
display_venn <- function(x, ...){
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}