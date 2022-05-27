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
#' @param session, the session
#' @param verbose, TRUE if you want results to be printed
#' @param visualize_prefix, the prefix of the PCA plot
#' 
#' @examples
#'     x <- runIterDE()
#'     
#' @export
#' 
runIterDE <- function(data = NULL, columns = NULL, conds = NULL, params = NULL, session = NULL, verbose = FALSE, visualize_prefix = "plot"){
  
  if (is.null(data)) return(NULL)
  data <- data[,columns]
  
  # parameters 
  ScoreMethod <- params[6]
  threshold <- params[7]
  iterde <- params[8]
  log2FC <- as.numeric(params[9])
  padj <- as.numeric(params[10])
  topstat <- as.numeric(params[11])

  # set variables and cutoff values
  cleaned_columns <- NULL
  DEgenes_new <- NULL
  DEgenes <- NULL
  g <- NULL
  iter <- 100
  
  # check for initial set of cleaned columns using cooks distance
  cur_columns <- setdiff(columns, cleaned_columns)
  cur_data <- data[, columns %in% cur_columns]
  cur_conds <- conds[columns %in% cur_columns]
  outlier_analysis <- getOutlierScores(data = cur_data, columns = cur_columns, conds = cur_conds)
  cleaned_columns <- outlier_analysis$cleaned_columns
  remaining_columns <- columns[!columns %in% cleaned_columns]
  outlier_scores <- outlier_analysis$Score
  
  # check iterations
  if(!is.null(session)){
    setProgress(value = (0 %% 11)/10,
                message = paste("Computational Profiling:  Iteration", 0, sep = " "),
                detail = paste("# of Removed Samples: ", length(cleaned_columns), sep = " "))
  }
  
  # iteration until convergence
  for(i in 1:iter){
    
    # select subset of genes and columns
    cur_columns <- setdiff(columns, cleaned_columns)
    cur_data <- data[, columns %in% cur_columns]
    cur_conds <- conds[columns %in% cur_columns]
    
    # adjust scoring threshold 
    if(threshold == "auto"){
      cut_threshold <- table(cur_conds)/length(cur_conds)
      param_threshold <- cut_threshold[cur_conds]
    } else{
      param_threshold <- as.numeric(threshold)
    }

    # DE analysis
    results <- runDE(data = cur_data, columns = cur_columns, conds = cur_conds, params = params)
    results$padj[is.na(results$padj)] <- 1
    if(iterde == "Stat."){
      Topresults <- results[order(results$stat,decreasing = TRUE)[1:topstat],]
    } else{
      Topresults <- results[abs(results$log2FoldChange) > log2FC & results$padj < padj, ]
    }
    DEgenes <- rownames(Topresults)
    DEgenes_new <- union(DEgenes_new, DEgenes) 

    # DEgenes of the current data, and normalize
    # data_de <- getNormalizedMatrix(cur_data, method = "TMM")
    data_de <- cur_data[rownames(cur_data) %in% DEgenes_new, ]
    
    # score with silhouette method
    if(ScoreMethod == "Silhouette"){
      
      # Calculate silhouette with spearman
      datax_spear <- (cor(data_de, method = "spearman") + 1)/2
      datax_spear_sim <- 1- datax_spear
      sil_spear <- silhouette(as.integer(factor(cur_conds)),as.dist(datax_spear_sim))
      
      # detect impure columns
      score <- (sil_spear[,3] + 1)/2
      exclude_list <- which(score < param_threshold) 
      
    }
    
    # score with NNLS method
    if(ScoreMethod == "NNLS-based"){
      
      # Expression Profiles
      profiles <- aggregate(t(data_de),list(cur_conds),mean)
      profiles <- t(profiles[,-1])
      
      # Deconvulate points based on expression profiles of conditions
      score <- NULL
      for(j in 1:ncol(data_de)){
        samples_DCV <- nnls(profiles,data_de[,j])
        norm_coef <- samples_DCV$x/sum(samples_DCV$x)
        score <- rbind(score,norm_coef)
      }
      colnames(score) <- unique(cur_conds)
      rownames(score) <- cur_columns
      
      # detect impure samples
      exclude_list <- NULL
      score_new <- NULL
      for(j in 1:nrow(score)){
        est_cond <- score[j,cur_conds[j]==unique(cur_conds)]
        if(est_cond < param_threshold) exclude_list <- c(exclude_list,j)
        score_new <- c(score_new, est_cond)
      }
      score <- score_new
    }
    
    # record the first score
    if(i == 1){
      score <- format(round(score, 3), nsmall = 3)
      score <- data.frame(Samples = colnames(cur_data), Conds = cur_conds, Scores = score)
      DEResults <- results
      DEscore <- score
      NonIterDEgenes <- DEgenes
    }
    
    # if there are no impure columns, stop
    if(length(exclude_list) == 0 | i == iter){
      NumberofIters <- i
      break
    }
    
    # update deleted columns
    cleaned_columns <- c(cleaned_columns, cur_columns[exclude_list])
    remaining_columns <- columns[!columns %in% cleaned_columns]
    
    # if there are less than two samples on each side, stop
    remaining_conds <- conds[!columns %in% cleaned_columns]
    remaining_conds <- table(remaining_conds)
    if(any(remaining_conds < 3)){
      NumberofIters <- i
      break
    }
    
    # print iteration results
    if(verbose){
      cat("# of Removed Samples: ", length(cleaned_columns), "\n")
      cat("# New DE genes:", length(DEgenes), "\n")
    }
    
    # visualize PCA
    if(!is.null(visualize_prefix)){
      datax_pr <- apply(data_de,1,scale)
      prtemp <- prcomp(datax_pr)
      datax_pr <- prtemp$x
      cur_cleaned <- ifelse(cur_columns %in% cur_columns[exclude_list], "remove", "keep")
      ggplot(data.frame(datax_pr), 
             ggplot2::aes(x = PC1, y = PC2, colour = cur_cleaned, shape = cur_conds), diag = "blank") + geom_point() + labs(title = paste0("# of Removed Samples: ", length(cleaned_columns)))
      ggsave(paste0(visualize_prefix, i, ".jpeg"), device = "jpeg", width = 7, height = 5)
    }
    
    # check iterations
    if(!is.null(session)){
      setProgress(value = (i %% 11)/10,
                  message = paste("Computational Profiling:  Iteration", i, sep = " "),
                  detail = paste("# New DE genes:", length(DEgenes),
                                 "# of Removed Samples: ", length(cleaned_columns), sep = " "))
    }
  }
  IterDEgenes <- DEgenes

  return(list(IterDEResults = results, DEResults = DEResults, IterDEgenes = IterDEgenes, DEgenes = NonIterDEgenes, 
              cleaned_columns = cleaned_columns, remaining_columns = remaining_columns, NumberofIters = NumberofIters, g = g))
}

#' getOutlierScores
#' 
#' Calculate membership scores given the iterative DE results
#'
#' @param deres DE results
#' @param data data
#' @param columns samples 
#' @param conds conditions
#' @param threshold threshold for score
#'
#' @examples
#'      x <- getOutlierScores()
#'     
#' @export
#'  
getOutlierScores <- function(data = NULL, columns = NULL, conds = NULL, threshold = 0.5){
  
  # set DESeq model
  data <- data[,columns]
  coldata <- prepGroup(conds, columns)
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = data,
                                colData = coldata,
                                design= ~ group)
  dds <- DESeq2::DESeq(dds)
  
  # get outlier scores
  cooks_dds <- apply(SummarizedExperiment::assays(dds)[["cooks"]], 2, function(x) quantile(x, probs = 0.75, na.rm = TRUE))
  Score <- 1-pf(cooks_dds, 2, ncol(data) - 2)
  cleaned_columns <- columns[Score < threshold]
  Score <- data.frame(Samples = colnames(data), Conds = conds, Scores = Score)
  
  return(list(cleaned_columns = cleaned_columns, Score = Score))
}

#' getFinalScores
#' 
#' Calculate membership scores given the iterative DE results
#'
#' @param deres DE results
#' @param data data
#' @param columns samples 
#' @param conds conditions
#' @param params DE params
#' @param ManualDEgenes manually inserted list of DE genes
#' @param TopStat the number of top statistics
#'
#' @examples
#'      x <- getFinalScores()
#'     
#' @export
#'  
getFinalScores <- function(deres = NULL, data = NULL, columns = NULL, conds = NULL, params = NULL, 
                           ManualDEgenes = NULL, TopStat = NULL){
 
  if (is.null(data)) return(NULL)
  data <- data[,columns]
  
  # parameters
  ScoreMethod <- params[6]
  TopStat <- as.numeric(TopStat)
  
  # DE genes
  if(is.na(TopStat)){
    if(is.null(ManualDEgenes)){
      DEgenes <- deres$DEgenes
      IterDEgenes <- deres$IterDEgenes
    } else {
      DEgenes <- IterDEgenes <- read.table(file = ManualDEgenes$datapath)[,1]
    }
  } else {
    if(TopStat > 5){
      TopStat <- as.numeric(TopStat)
      DEgenes <- rownames(deres$DEResults[order(deres$DEResults$stat, decreasing = TRUE)[1:TopStat],])
      IterDEgenes <- rownames(deres$IterDEResults[order(deres$IterDEResults$stat, decreasing = TRUE)[1:TopStat],]) 
    } else {
      DEgenes <- deres$DEgenes
      IterDEgenes <- deres$IterDEgenes
    }
  }
  
  # set variables and cutoff values
  cleaned_columns <- deres$cleaned_columns

  # select subset of genes and columns
  cur_columns <- setdiff(columns, cleaned_columns)
  cur_data <- data[, columns %in% cur_columns]
  cur_data <- data[, columns %in% cur_columns]
  cur_conds <- conds[columns %in% cur_columns]
  
  # Final Scores of Heterogeneous Groups
  # data_de <- getNormalizedMatrix(data, method = "TMM")
  data_de <- data[rownames(data) %in% DEgenes, ]
  
  if(ScoreMethod == "Silhouette"){
    
    # Spearman correlations
    datax_spear <- (cor(data_de, method = "spearman") + 1)/2
    datax_spear_sim <- 1- datax_spear
    DEscore <- custom_silhouette(conds, datax_spear_sim)
    
  }
  if(ScoreMethod == "NNLS-based"){
    
    # Expression Profiles
    profiles <- data_de
    profiles <- aggregate(t(profiles),list(conds),mean)
    profiles <- t(profiles[,-1])
    
    # Deconvulate points based on expression profiles of conditions
    DEscore <- NULL
    for(j in 1:ncol(data_de)){
      samples_DCV <- nnls(profiles,data_de[,j])
      norm_coef <- samples_DCV$x/sum(samples_DCV$x)
      DEscore <- rbind(DEscore,norm_coef)
    }
    colnames(DEscore) <- unique(conds)
    rownames(DEscore) <- columns
    
  }

  DEscore_new <- NULL
  for(i in 1:nrow(DEscore)){
    DEscore_new <- c(DEscore_new,
                     DEscore[i,conds[i]==unique(conds)])
  }
  DEscore <- data.frame(Samples = colnames(data), Conds = conds, Scores = DEscore_new)
  
  # Final Scores of Heterogeneous groups
  # data_de <- getNormalizedMatrix(data, method = "TMM")
  data_de <- data[rownames(data) %in% IterDEgenes, ]
  if(ScoreMethod == "Silhouette"){
    
    # Spearman correlations
    datax_spear <- (cor(data_de, method = "spearman") + 1)/2
    datax_spear_sim <- 1- datax_spear
    datax_spear_sim_clean <- datax_spear_sim[,columns %in% cur_columns]
    IterDEscore <- custom_silhouette(cur_conds, datax_spear_sim_clean)
    
  }
  if(ScoreMethod == "NNLS-based"){
    
    # Expression Profiles
    profiles <- data_de[,colnames(data_de) %in% cur_columns]
    profiles <- aggregate(t(profiles),list(cur_conds),mean)
    profiles <- t(profiles[,-1])
    
    # Deconvulate points based on expression profiles of conditions
    IterDEscore <- NULL
    for(j in 1:ncol(data_de)){
      samples_DCV <- nnls(profiles,data_de[,j])
      norm_coef <- samples_DCV$x/sum(samples_DCV$x)
      IterDEscore <- rbind(IterDEscore,norm_coef)
    }
    colnames(IterDEscore) <- unique(conds)
    rownames(IterDEscore) <- columns
    
  }
  # IterDEscore <- format(round(IterDEscore, 3), nsmall = 3)
  IterDEscore_new <- NULL
  for(i in 1:nrow(IterDEscore)){
    IterDEscore_new <- c(IterDEscore_new,
                         IterDEscore[i,conds[i]==unique(conds)])
  }
  IterDEscore <- data.frame(Samples = colnames(data), Conds = conds, Scores = IterDEscore_new)

  return(list(IterDEscore = IterDEscore, DEscore = DEscore))
}

#' custom_silhouette
#' 
#' A custom silhouette function for calculating the silhouette measure of each point with respect
#' to each condition
#'
#' @param x conditions
#' @param dist distance matrix
#'
#' @examples
#'      x <- custom_silhouette()
#'     
#' @export
#'       
custom_silhouette <- function(x = NULL,dist = NULL){
  if (is.null(x)) return(NULL)
  
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