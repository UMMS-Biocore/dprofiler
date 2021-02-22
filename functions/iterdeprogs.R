#' debrowseriterdeanalysis
#'
#' Module to perform and visualize DE results.
#' 
#' @param input, input variables
#' @param output, output objects
#' @param session, session 
#' @param data, a matrix that includes expression values
#' @param columns, columns
#' @param conds, conditions
#' @param params, de parameters
#' @return DE panel 
#' @export
#'
#' @examples
#'     x <- debrowseriterdeanalysis()
#'
debrowseriterdeanalysis <- function(input = NULL, output = NULL, session = NULL, 
                                data = NULL, columns = NULL, conds = NULL, params = NULL) {
  if(is.null(data)) return(NULL)
  iterderes <- reactive({
    runIterDE(data, columns, conds, params)
  })
  iterprepDat <- reactive({
    applyFiltersNew(addDataCols(data, iterderes()$results, columns, conds), input)
  })
  observe({
    iterdat <-  iterprepDat()[iterprepDat()$Legend == input$legendradio,]
    iterdat2 <- removeCols(c("ID", "x", "y","Legend", "Size"), iterdat)
    getTableDetails(output, session, "IterDEResults", iterdat2, modal = FALSE)
    
    getTableDetails(output, session, "MembershipScores", iterderes()$score, modal = FALSE)
  })
  list(dat = iterprepDat)
}

#' getIterDEResultsUI
#' Creates a panel to visualize DE results
#'
#' @param id, namespace id
#' @return panel
#' @examples
#'     x <- getIterDEResultsUI("batcheffect")
#'
#' @export
#'
getIterDEResultsUI<- function (id) {
  ns <- NS(id)
  list(
    fluidRow(
      shinydashboard::box(title = "Results",
                          solidHeader = T, status = "info",  width = 12, collapsible = TRUE,
                          fluidRow(
                            column(12,
                                   uiOutput(ns("IterDEResults"))
                            ),
                            # actionButtonDE("goMain", "Go to Main Plots", styleclass = "primary")
                          )
      ),
      shinydashboard::box(title = "Membership Scores",
                          solidHeader = T, status = "info",  width = 12, collapsible = TRUE,
                          fluidRow(
                            column(12,
                                   uiOutput(ns("MembershipScores"))
                            )
                          )
      )
    )
  )
}


# runIterDE <- function(data = NULL, columns = NULL, conds = NULL, params = NULL) {
#   if (is.null(data)) return(NULL)
#   de_res <- NULL
# 
#   print(conds)
#   print(columns)
#   
#   data.frame(de_res)
# }

#' runIterDE
#' 
#' Run Iterative DE algorithms on the selected parameters.  Output is
#' to be used for the interactive display.
#' 
#' @param data, A matrix that includes all the expression raw counts,
#'     rownames has to be the gene, isoform or region names/IDs
#' @param columns, is a vector that includes the columns that are going
#'     to be analyzed. These columns has to match with the given data.
#' @param conds, experimental conditions. The order has to match
#'     with the column order
#' @param params, all params for the DE methods
#' @return de results
#' 
#' @export
#' 
#' @examples
#'     x <- runDE()
#'     
runIterDE <- function(data = NULL, columns = NULL, conds = NULL, params = NULL){
  
  if (is.null(data)) return(NULL)
  data <- data[,columns]
  
  # set variables and cutoff values
  cleaned_columns <- NULL
  DEgenes <- NULL
  data_tmm <- getNormalizedMatrix(data, method="TMM")
  log2FC <- 0.5
  padj <- 0.05
  iter <- 100
  
  # iteration until convergence
  for(i in 1:iter){
    
    # select subset of genes and columns
    cur_columns <- setdiff(columns, cleaned_columns)
    cur_data <- data[, columns %in% cur_columns]
    cur_data_tmm <- data_tmm[, columns %in% cur_columns]
    cur_conds <- conds[columns %in% cur_columns]
    
    # DE analysis
    results <- runDE(data = cur_data, columns = cur_columns, conds = cur_conds, params = params)
    results$padj[is.na(results$padj)] <- 1
    DEgenes_new <- rownames(results[abs(results$log2FoldChange) > log2FC & results$padj < padj, ])
    DEgenes <- union(DEgenes, DEgenes_new) 
    print(length(DEgenes))
    
    # DEgenes of the current data
    data_de <- cur_data_tmm[rownames(cur_data) %in% DEgenes, ]
    
    # Calculate silhouette with spearman
    datax_spear <- (cor(data_de, method = "spearman") + 1)/2
    datax_spear_sim <- 1- datax_spear
    sil_spear <- silhouette(as.integer(factor(cur_conds)),as.dist(datax_spear_sim))
    
    # detect impure columns
    exclude_list <- which(sil_spear[,3] < 0)
    
    # if there are no impure columns, stop
    if(length(exclude_list) == 0 | i == iter){
      break
    }
    
    # update deleted columns
    cleaned_columns <- c(cleaned_columns, cur_columns[exclude_list])
    print(length(cleaned_columns))
  }
  
  # # select subset of genes and columns
  # cur_columns <- setdiff(columns, cleaned_columns)
  # cur_data <- data[, columns %in% cur_columns]
  # cur_covariate <- covariate[columns %in% cur_columns]
  # cur_conds <- conds[columns %in% cur_columns]
  # 
  # # Final DE analysis 
  # if(!is.null(covariate)){
  #   results <- getDEgenes(data = cur_data, conds = cur_conds, 
  #                         covariate = cur_covariate, DEtest = DEtest)
  # } else {
  #   results <- getDEgenes(data = cur_data, conds = cur_conds, DEtest = DEtest)
  # }
  # DEgenes_new <- results$DEgenes
  # DEgenes <- cbind(union(DEgenes, DEgenes_new))
  
  # Final Silhouette scoring of all columns
  data_de <- data_tmm[rownames(data_tmm) %in% DEgenes,]
  datax_spear <- (cor(data_de, method = "spearman") + 1)/2
  datax_spear_sim <- 1- datax_spear
  datax_spear_sim_clean <- datax_spear_sim[,columns %in% cur_columns]
  score <- custom_silhouette(cur_conds, datax_spear_sim_clean)
  score <- data.frame(score, conds = conds, sample = colnames(data_tmm))

  return(list(results = results, score = score))
}


custom_silhouette <- function(x,dist){
  
  # Classes and Membership
  unique_class <- unique(x)
  membership <- NULL
  
  # Check Silhouette of all columns
  for(i in 1:nrow(dist)){
    
    avg_dist <- aggregate(as.data.frame(dist[i,-i]),list(x[-i]),mean)
    
    sil <- rep(0,nrow(avg_dist))
    for(j in 1:nrow(avg_dist)){
      a <- avg_dist[j,2]
      b <- min(avg_dist[-j,2])
      sil[j] <- (b-a)/max(a,b)
    }
    membership <- rbind(membership, sil[match(avg_dist$Group.1,unique_class)])
    
  }
  membership <- (membership+1)/2
  colnames(membership) <- unique_class
  return(membership)
}