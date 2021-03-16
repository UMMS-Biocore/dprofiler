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

getIterDESummary <- function(output = NULL, session = NULL, vennname = NULL, summaryname = NULL, 
                             deres = NULL, params = NULL){
  
  # Iteration Summary
  output[[summaryname]] <- renderTable({
    result <- rbind(deres$NumberofIters, length(deres$cleaned_columns), 
                    length(deres$DEgenes), length(deres$IterDEgenes))
    rownames(result) <- c("# of Iterations", "# of Hetero. Samples", 
                          "# of Initial DE genes", "# of Final DE genes ")
    colnames(result) <- "Value"
    result
  }, digits = 0, rownames = TRUE, align = "lc")
  
  # Scoring parameters 
  scoresummaryname <- paste0(summaryname,"Iter")
  output[[scoresummaryname]] <- renderTable({
    result <- params[6:9]
    if(params[9] == "Stat."){
      result <- matrix(c(result, params[10:11]),ncol = 1)
      rownames(result) <- c("Score Method", "Min. Score", "Normalization","Selection Method", 
                            "Log2FC", "P-adj")
    } else {
      result <- matrix(c(result, params[12]),ncol = 1)
      rownames(result) <- c("Score Method", "Min. Score", "Normalization","Selection Method", 
                            "# of Top Statistics")
    }
    colnames(result) <- "Value"
    result
  }, digits = 0, rownames = TRUE, align = "lc")
  
  # DE parameters
  desummaryname <- paste0(summaryname,"DE")
  output[[desummaryname]] <- renderTable({

    if(params[1] == "DESeq2"){
      result <- matrix(params[1:5],ncol = 1)
      rownames(result) <- c("DE Method", "Fit Type", "Beta Prior","Test Type", "Shrinkage")
    } else if(params[1] == "EdgeR") {
      result <- matrix(params[1:4],ncol = 1) 
      rownames(result) <- c("DE Method", "Normalization", "Dispersion", "Test Type")
    } else{
      result <- matrix(params[1:4],ncol = 1)
      rownames(result) <- c("DE Method", "Normalization", "Fit Type", "Norm.Bet.Arrays")
    }
    colnames(result) <- "Value"
    result
  }, digits = 0, rownames = TRUE, align = "lc")
  
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

display_venn <- function(x, ...){
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}