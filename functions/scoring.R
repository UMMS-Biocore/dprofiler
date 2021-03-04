getScoreDetails <- function(output = NULL, session = NULL, 
                                 plotname = NULL, DEscores = NULL, IterDEscores = NULL) {
  output[[plotname]] <- renderPlotly({
    dat <- rbind(
      data.frame(DEscores, type = "Homogeneous"),
      data.frame(IterDEscores, type = "Heterogeneous")
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
  })
}

getIterDESummary <- function(output = NULL, session = NULL, vennname = NULL, summaryname = NULL, 
                             deres = NULL){
  output[[summaryname]] <- renderTable({
    result <- rbind(deres$NumberofIters, length(deres$cleaned_columns), 
                    length(deres$DEgenes), length(deres$IterDEgenes))
    rownames(result) <- c("# of Iterations", "# of Hetero. Samples", 
                          "# of Initial DE genes", "# of Final DE genes ")
    colnames(result) <- "Value"
    result
  }, digits = 0, rownames = TRUE, align = "lc")
  
  x <- list(
    A = deres$DEgenes, 
    B = deres$IterDEgenes
  )
  output[[vennname]] <- renderPlot({
    display_venn(x,
                 category.names = c("Before" , "After"),
                 fill = c("#999999", "#E69F00"))
  })
}

display_venn <- function(x, ...){
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}