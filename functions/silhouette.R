getSilhouetteDetails <- function(output = NULL, session = NULL, 
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

getVennDiagram <- function(output = NULL, session = NULL, plotname = NULL, 
                           IterDEgenes = NULL, DEgenes= NULL){
  x <- list(
    A = DEgenes, 
    B = IterDEgenes
  )
  output[[plotname]] <- renderPlot({
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