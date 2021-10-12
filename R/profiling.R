#' dprofilerprofiling
#'
#' Module to perform and visualize deconvolution results.
#'
#' @param input, input variables
#' @param output, output objects
#' @param session, session
#' @param dc, results of iterative DE Analysis
#' @param profiledata expression matrix of profile data
#' @param profilemetadata metadata of profile data
#' @param parent_session parent session
#'
#' @examples
#'     x <- dprofilerprofiling()
#'
dprofilerprofiling <- function(input = NULL, output = NULL, session = NULL, dc = NULL, profiledata = NULL,
                               profilemetadata = NULL, parent_session = NULL) {
  if(is.null(dc)) return(NULL)

  # get details
  prof_details <- getProfilingDEparameter(profiledata, profilemetadata, session, input)
  
  # DE Results
  deres <- reactive({
    waiter_show(html = spin_ring(), color = transparent(.5))
    withProgress(message = 'Running Profiling', value = 0, {
      deresults <- runMultipleDE(profiledata, prof_details$cols, prof_details$conds, prof_details$params)
    })
    updateTabsetPanel(session = parent_session, "ProfilingResults", "profilingscores")
    waiter_hide()
    deresults
  })
  
  # Apply Filters to DE Results
  prepDat <- reactive({
    applyFiltersNew(addDataCols(profiledata, deres(), prof_details$cols, prof_details$conds), input)
  })
  
  # # Membership Scores
  # scores <- reactive({
  #     getProfileScores(dc, profiledata[deres()$], profilemetadata)
  # })
  
  # Observe for Tables and Plots
  observe({
    
    # prepare DE tables
    dat <-  prepDat()[prepDat()$padj < prof_details$params[11] & abs(prepDat()$log2FoldChange) > prof_details$params[10],]
    dat2 <- removeCols(c("ID", "x", "y","Legend", "Size"), dat)
    
    # DE Results
    getTableDetails(output, session, "DEResults", dat2, modal=FALSE)

    # get scores
    scores <- getProfileScores(dc, dat, profiledata, profilemetadata)
    getProfileScoreDetails(output, session, "MembershipScores", scores)
    
  })
  
  return(NULL)
}

#' getProfilingUI
#' 
#' Creates a panel to visualize Profiling results
#'
#' @param id, namespace id
#' 
#' @examples
#'     x <- getProfilingUI("profiling")
#'
getProfilingUI <- function (id) {
  ns <- NS(id)
  list(
    tabBox(id = "ProfilingResults",
           width = NULL,
           tabPanel(title = "Conditions",
                    shinydashboard::box(title = "Comparison Selection",
                                        solidHeader = TRUE, status = "info",  width = NULL, height = NULL, collapsible = TRUE,
                                        fluidRow(
                                          uiOutput(ns("conditionSelector")),
                                          column(12,
                                                 actionButtonDE("startprofiling", "Start", styleclass = "primary")
                                          )
                                        ))
           ),
           tabPanel(title = "DE Analysis",
                    fluidRow(
                      shinydashboard::box(title = "Differentially Expressed Genes",
                                          solidHeader = T, status = "info",  width = 9, collapsible = TRUE,
                                          uiOutput(ns("DEResults"))
                      ),
                    )
           ),
           tabPanel(title = "Profiling Results",
                    fluidRow(
                      shinydashboard::box(title = "Membership Scores",
                                          solidHeader = T, status = "info",  width = 9, collapsible = TRUE,
                                          plotlyOutput(ns("MembershipScores"))
                      )
                    ),
                    value = "profilingscores"
           )
    )
  )
}


#' getProfileScoreDetails
#' 
#' Details and scores of the profiling analysis
#'
#' @param output output 
#' @param session session
#' @param plotname name of the plot
#' @param scores profiling scores
#'
#' @examples
#'     x <- getProfileScoreDetails()
#'     
getProfileScoreDetails <- function(output = NULL, session = NULL, plotname = NULL, scores = NULL) {
  if(is.null(output)) return(NULL)
  
  output[[plotname]] <- renderPlotly({
    p <- ggplot(data=scores, aes(x = reorder(Samples,Scores), y = Scores)) +
      geom_bar(aes(fill = profileConds), stat="identity") + 
      facet_grid(. ~ Conds) + 
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


#' getProfileScores
#' 
#' Calculate the profiling scores
#'
#' @param dc DE results 
#' @param dat expression matrix of the main data 
#' @param profiledata expression matrix of the profiling data
#' @param profilemetadata metadata of the profiling data
#'
#' @examples
#'      x <- getProfileScores()
#' 
getProfileScores <- function(dc = NULL, dat = NULL, profiledata = NULL, profilemetadata = NULL){
  if(is.null(dc)) return(NULL)
  
  # unique genes of profile data
  unique_genes <- unique(dat$ID)
  profiledata <- profiledata[unique_genes,]
    
  # get data 
  remaining_columns <- dc$cols[!dc$cols %in% dc$cleaned_columns]
  conds <- dc$conds[!dc$cols %in% dc$cleaned_columns]
  data <- dc$init_dedata[dc$IterDEgenes, remaining_columns]
  data <- t(apply(data,1,scale))
  colnames(data) <- remaining_columns
  
  # get overlapping genes
  genes <- intersect(rownames(data), rownames(profiledata))
  data <- data[rownames(data) %in% genes,]
  profiledata <- profiledata[rownames(profiledata) %in% genes,]
  colnames_profiledata <- colnames(profiledata)
  profiledata <- t(apply(profiledata,1,scale))
  colnames(profiledata) <- colnames_profiledata
  
  # get scores
  datax_spear <- (cor(profiledata, data, method = "spearman") + 1)/2
  rownames(datax_spear) <- colnames(profiledata)
  colnames(datax_spear) <-  colnames(data)
  datax_spear_sim <- 1- datax_spear
  scores <- (external_silhouette(conds, datax_spear) + 1)/2
  rownames(scores) <- unique(conds)
  scores <- melt(scores)
  
  # get score data frame
  scores <- data.frame(Conds = scores[,1], Samples = scores[,2], Scores = scores[,3],
                       profileConds = rep(profilemetadata$treatment, each = 2))
  
  return(scores)
}

#' getExpressionProfiles
#'
#' Given the gene expression matrix with multiple conditions, generate expression profiles.
#' 
#' @param deres DE results
#' @param data expression data
#' @param columns columns
#' @param conds conditions
#'
#' @examples
#'      x <- getExpressionProfiles()
#'      
getExpressionProfiles <- function(deres = NULL, data = NULL, columns = NULL, conds = NULL){
  if(is.null(deres)) return(NULL)
  
  # get remaining columns and data
  remaining_columns <- (columns != deres$cleaned_columns)
  if(length(remaining_columns) > 0){
    columns <- columns[remaining_columns]
    conds <- conds[remaining_columns]
  }
  data <- data[deres$IterDEgenes, columns]
  
  # calculate profiles
  profiles <- aggregate(t(data),list(conds),mean)
  colnames_profiles <- profiles[,1]
  profiles <- data.frame(t(as.matrix(profiles[,-1])))
  colnames(profiles) <- colnames_profiles
  
  return(profiles)
}


#' external_silhouette
#' 
#' Calculate silhouette measure for non-symmetric distance matrices. 
#' Suitable for comparison of expression profiles across different datasets and experiments 
#'
#' @param cluster cluster labels of profile data
#' @param dist2 non-symmetric similarity matrix
#'
#' @examples
#'      x <- external_silhouette()
external_silhouette <- function(cluster = NULL, dist2 = NULL){
  if(is.null(cluster)) return(NULL)
  
  cls <- levels(cluster)
  allmeans <- apply(dist2,1,function(x){
    aggdata <- aggregate(x,list(cluster),mean)
    temp <- aggdata[,1]
    aggdata <- aggdata[,-1,drop=FALSE]
    rownames(aggdata) <- temp
    return(aggdata$x)
  })
  allsil <- NULL
  for(i in 1:nrow(allmeans)){
    a <- allmeans[i,,drop=FALSE]
    b <- apply(allmeans[-i,,drop=FALSE],2,min)
    sil <- (b-a)/pmax(a,b)
    allsil <- rbind(allsil,sil)
  }
  rownames(allsil) <- cls
  return(allsil)
}

#' getProfilingDEparameter
#' 
#' Prepare DE analysis parameters for the profiling data
#'
#' @param data data
#' @param metadata metadata
#' @param session session
#' @param input input 
#'
#' @examples
#'      x <- getProfilingDEparameter()
#'      
getProfilingDEparameter <- function(data = NULL, metadata = NULL, session = NULL, input = NULL) {
  
    if (is.null(data)) return(NULL)
    
    inputconds <- reactiveValues(demethod_params = list(), conds = list(), dclist = list())
    cnt <- 1
    # inputconds$conds <- list()
    # inputconds$conds[[1]] <- isolate(input[[paste0("condition",1)]])
    # inputconds$conds[[2]] <- isolate(input[[paste0("condition",2)]])

    #Get parameters for each method
    inputconds$demethod_params <- NULL
    if (isolate(input[[paste0("demethod",cnt)]]) == "DESeq2"){
      inputconds$demethod_params[cnt] <- paste(
        isolate(input[[paste0("demethod",cnt)]]),
        isolate(input[[paste0("fitType",cnt)]]),
        isolate(input[[paste0("betaPrior",cnt)]]),
        isolate(input[[paste0("testType",cnt)]]),
        isolate(input[[paste0("shrinkage",cnt)]]), sep=",")
    }
    else if (isolate(input[[paste0("demethod",cnt)]]) == "EdgeR"){
      inputconds$demethod_params[cnt]<- paste(
        isolate(input[[paste0("demethod",cnt)]]),
        isolate(input[[paste0("edgeR_normfact",cnt)]]),
        isolate(input[[paste0("dispersion",cnt)]]),
        isolate(input[[paste0("edgeR_testType",cnt)]]), 
        "", sep=",")
    }
    else if (isolate(input[[paste0("demethod",cnt)]]) == "Limma"){
      inputconds$demethod_params[cnt] <- paste(
        isolate(input[[paste0("demethod",cnt)]]),
        isolate(input[[paste0("limma_normfact",cnt)]]),
        isolate(input[[paste0("limma_fitType",cnt)]]),
        isolate(input[[paste0("normBetween",cnt)]]), 
        "", sep=",")
    }
    
    # condition inputs for iterative de analysis
    inputconds$demethod_params[cnt] <- paste(inputconds$demethod_params[cnt],
                                             isolate(input[[paste0("scoremethod",cnt)]]),
                                             isolate(input[[paste0("minscore",cnt)]]),
                                             isolate(input[[paste0("iterde_norm",cnt)]]), 
                                             isolate(input[[paste0("iterde",cnt)]]),
                                             isolate(input[[paste0("logfoldchange",cnt)]]),
                                             isolate(input[[paste0("padj",cnt)]]),
                                             isolate(input[[paste0("topstat",cnt)]]), sep = ","
    )
    
    # conds <- c(rep(paste0("Cond", 1), 
    #                length(inputconds$conds[[1]])), 
    #            rep(paste0("Cond", 2), length(inputconds$conds[[2]])))
    # cols <- c(paste(inputconds$conds[[1]]), 
    #           paste(inputconds$conds[[2]]))
    conds <- metadata[[input$conditions_from_meta1]]
    cols <- colnames(data)
    
    params <- unlist(strsplit(inputconds$demethod_params[1], ","))
    
    return(list(conds = conds, cols = cols, params = params))
}