#' dprofilerprofiling
#'
#' Module to perform and visualize deconvolution results.
#'
#' @param input, input variables
#' @param output, output objects
#' @param session, session
#' @param dc, results of iterative DE Analysis
#' @param metadatatable, the meta data table
#' @param profiledata expression matrix of profile data
#' @param profilemetadata metadata of profile data
#' @param profmarkertable profiling DE genes table
#' @param parent_session parent session
#'
#' @examples
#'     x <- dprofilerprofiling()
#'     
#' @export
#' 
dprofilerprofiling <- function(input = NULL, output = NULL, session = NULL, parent_session = NULL, data = NULL) {

  startProf <- data$startProf
  data <- data$selectedData
  
  ###
  # Reactive Values ####
  ###
  
  ProfResults <- reactiveValues(score = NULL)
  ScoreTable <- reactiveVal()
  Summarise <- reactive(input$summarisescores)
  
  ###
  # Start Profiling ####
  ###
  
  # # metadata
  # sample_column_ind <- which(apply(metadatatable, 2, function(x) sum(x %in% columns) == length(columns)))
  # sample_id <- colnames(metadatatable)[sample_column_ind]
  # # metadatatable <- metadatatable[metadatatable[,sample_id] %in% columns,]
  # metadatatable <- metadatatable[match(columns, metadatatable[,sample_id]),]
  
  observeEvent(startProf(), {
    waiter_show(html = spin_ring(), color = transparent(.5))
    withProgress(message = 'Running Comparative Profiling', value = 0, {
      prof_genes <- getAllProfGenes(data()$markers, norm_counts(), input)
      ProfResults$score <- getProfileScores(norm_counts(), data()$meta, prof_genes, data()$prof_count, data()$prof_meta, input)
    })
    waiter_hide()
    showTab("menutabs","discover", session = parent_session)
    showTab("ProfilingResults","profilingscores", session = parent_session)
    updateTabsetPanel(session = parent_session, "ProfilingResults", "profilingscores")
  })
  
  # counts, normalization if asked
  norm_counts <- reactive({
    results <- getNormalizedMatrix(data()$count, method = isolate(input$profiling_norm))
    results
  })
  
  ScoreTable <- reactive({
    if(is.null(ProfResults$score)) return(NULL)
    
    scoretable <- ProfResults$score
    scoretable <- scoretable$cross_scores
    
    # if(!is.null(data()$score)){
    #   scoretable <- cbind(data()$score, scoretable)
    # } 
    data.frame(Sample = rownames(scoretable), scoretable[,!colnames(scoretable) %in% "Sample"], row.names = NULL)
  })
  
  ###
  # Heatmap Events ####
  ###
  
  # prepare heat data
  data_de_tmm <- reactive({
    prof_genes <- getAllProfGenes(data()$markers, data()$count, input)
    if(is.null(prof_genes)) return(NULL)
    heatdata <- prepHeatData(data()$count, input)
    as.matrix(heatdata[prof_genes,])
  })

  output$heatmap <- renderPlotly({
    if(!is.null(data_de_tmm())){
      if(nrow(data_de_tmm()) > 1 & ncol(data_de_tmm()) > 1)
        withProgress(message = 'Drawing Heatmap', detail = "interactive", value = 0, {
          runHeatmap(input, session, data_de_tmm())
        })
    }
  })

  output$heatmap2 <- renderPlot({
    if(!is.null(data_de_tmm())){
      if(nrow(data_de_tmm()) > 1 & ncol(data_de_tmm()) > 1)
        withProgress(message = 'Drawing Heatmap', detail = "non-interactive", value = 0, {
          runHeatmap2(input, session, data_de_tmm())
        })
    }
  })

  output$heatmapUI <- renderUI({
    if (is.null(input$interactive)) return(NULL)
    column(6,
           shinydashboard::box(
             collapsible = TRUE, title = "Heatmaps", status = "info",
             solidHeader = TRUE, width = NULL,
             draggable = TRUE,
             column(12,getPlotArea(input, session)),
           )
    )
  })

  # Condition Selector for Reference Profile Data
  output$identSelector <- renderUI({

    if(!is.null(data()$prof_meta)){
      character_columns <- colnames(data()$prof_meta)[!sapply(as.data.frame(data()$prof_meta),is.numeric)]
      character_columns <- character_columns[sapply(as.data.frame(data()$prof_meta[,character_columns]),function(x) length(unique(x)) < round(length(x)/4))]
      selected_columns <- character_columns[sapply(as.data.frame(data()$prof_meta[,character_columns]),function(x) length(unique(x)) > 1 && length(unique(x)) < 5)]
      if(is.null(selected_columns)) selected_columns <- character_columns
      list(
        column(5,
               selectInput(session$ns("selectcondpca"), label = "Select Ident",
                           choices = character_columns, selected = selected_columns[1]))
      )
    }
  })
  
  ###
  # Window Settings ####
  ###

  # hide initial tabs
  hideTab("ProfilingResults","profilingscores", session = parent_session)
  
  # Observe for Tables and Plots
  observe({
    
    # get scores
    getProfileScoreDetails(output, session, "MembershipScores", ScoreTable(), modal = FALSE, highlight = TRUE)

    # visualize data
    # getPCAEmbedding(output, "uploadSummaryPCA", profilemetadata, ident = profileConds)
    getCCAEmbedding(output, "uploadSummaryPCA", ProfResults$score$ccametadata, input)
  })
  
  return(list(Summarise = Summarise, ScoreTable = ScoreTable))
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
#' @export
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
                                                 actionButtonDE(ns("startprofiling"), "Start", styleclass = "primary")
                                          )
                                        ))
           ),
           tabPanel(title = "Profiling Results",
                    # column(12, 
                    #        p(strong("Note:")," Here, we use similarity measures based on silhouette measure and non-negative least squares to calculate ", strong("the membership score of Vitiligo samples using another reference
                    #          bulk Vitiligo dataset,"), " revealing similarities of lesional and non-lesional samples across datasets.")
                    # ),
                    fluidRow(
                      shinydashboard::box(title = "Membership Scores",
                                          solidHeader = T, status = "info",  width = 12, collapsible = TRUE,
                                          DT::dataTableOutput(ns("MembershipScores")),
                                          column(2,actionButtonDE(ns("summarisescores"), "Summarise Scores", styleclass = "primary", 
                                                                  style = "width: 100%; margin-top: 25px; margin-left: 0px"))
                      )
                    ),
                    fluidRow(
                      uiOutput(ns("heatmapUI")),
                      shinydashboard::box(title = "Joint Embedding CCA",
                                          solidHeader = T, status = "info",  width = 6, collapsible = TRUE,
                                          fluidRow(
                                            column(12,div(style = 'overflow: scroll',
                                                          plotOutput(ns("uploadSummaryPCA")))
                                            )                                          
                                          )
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
#' @param output, output
#' @param session, session
#' @param tablename, table name
#' @param data, matrix data
#' @param modal, if it is true, the matrix is going to be in a modal
#' @param highlight if it is true, numerical columns are highlighted
#'
#' @examples
#'     x <- getProfileScoreDetails()
#'     
#' @export
#'     
getProfileScoreDetails <- function(output = NULL, session = NULL, tablename = NULL, data = NULL, 
                                   modal = NULL, highlight = FALSE){
  if (is.null(data)) return(NULL)
  output[[tablename]] <- DT::renderDataTable({
    if (!is.null(data)){
      dttable <- DT::datatable(data, extensions = 'Buttons',
                               rownames = FALSE,
                               options = list( server = TRUE,
                                               dom = "Blfrtip",
                                               buttons = 
                                                 list("copy", list(
                                                   extend = "collection"
                                                   , buttons = c("csv", "excel", "pdf")
                                                   , text = "Download"
                                                 ) ), # end of buttons customization
                                               
                                               # customize the length menu
                                               lengthMenu = list( c(10, 20,  50, -1) # declare values
                                                                  , c(10, 20, 50, "All") # declare titles
                                               ), # end of length Menu customization
                                               pageLength = 10))
      numeric_names <- colnames(data[,sapply(as.data.frame(data), is.numeric), drop = FALSE])
      numeric_names <- numeric_names[!numeric_names %in% "Reads"]
      dttable <- dttable %>% DT::formatRound(numeric_names, digits=3)
      if(highlight){
        colours <- rainbow(length(numeric_names))
        for(i in 1:length(numeric_names)){
          dttable <-  dttable %>% DT::formatStyle(numeric_names[i],
                                                  background = DT::styleColorBar(c(0,1), colours[i]))
        }
      } 
      dttable
    }
  })
}

#' getAllProfGenes
#'
#' @param prof_marker_table comparative profiling marker table
#' @param data submitted bulk data
#' @param input input 
#'
#' @examples
#'     x <- getAllProfGenes()
#'     
#' @export
#'  
getAllProfGenes <- function(prof_marker_table = NULL, data = NULL, input = NULL){
  if (is.null(prof_marker_table)) return(NULL)
  
  # pull those genes that are in bulk data
  prof_table_genes <- unique(prof_marker_table$gene)
  common_genes <- intersect(prof_table_genes, rownames(data))
  prof_marker_table <- prof_marker_table[prof_marker_table$gene %in% common_genes,]
  
  # select cell types, and other conditions
  prof_marker_table <- prof_marker_table[prof_marker_table$Level %in% isolate(input$conditions_from_meta0),]
  prof_marker_table <- prof_marker_table[abs(prof_marker_table$log2FoldChange) > as.numeric(isolate(input$logFC)) & 
                                           prof_marker_table$padj < as.numeric(isolate(input$padj)),]
  
  # take out additional genes
  prof_marker_table <- prof_marker_table[!grepl("^AC[0-9]|^AL[0-9]|^AP[0-9]|^MT-", prof_marker_table$gene),]
  
  # grep cell type specific markers
  top_n_markers <- isolate(input$top_genes)
  prof_marker_table_up <- prof_marker_table[prof_marker_table$log2FoldChange > as.numeric(isolate(input$logFC)), ]
  prof_up_genes <- prof_marker_table_up$gene[order(prof_marker_table_up$log2FoldChange, decreasing = TRUE)[1:top_n_markers]]
  prof_marker_table_down <- prof_marker_table[prof_marker_table$log2FoldChange < -1*as.numeric(isolate(input$logFC)), ]
  prof_down_genes <- prof_marker_table_down$gene[order(prof_marker_table_down$log2FoldChange, decreasing = FALSE)[1:top_n_markers]]
  prof_genes <- c(prof_up_genes, prof_down_genes)
  
  return(prof_genes)
}

#' getProfileScores
#'
#' Calculate the profiling scores
#'
#' @param data Bulk data
#' @param degenes DE genes
#' @param profiledata expression matrix of the profiling data
#' @param profilemetadata metadata of the profiling data
#' @param input input
#'
#' @examples
#'      x <- getProfileScores()
#'
#' @export
#'
getProfileScores <- function(data = NULL, metadata = NULL, degenes = NULL, profiledata = NULL, profilemetadata = NULL, input = NULL){
  if(is.null(data)) return(NULL)
  
  # get params
  # params <- params$params
  # ScoreMethod = if (!is.null(params[1]))
  #   params[1]
  # profileCondsName = if (!is.null(params[2]))
  #   params[2]
  ScoreMethod <- isolate(input[[paste0("scoremethod",1)]])
  profileCondsName <- strsplit(isolate(input[[paste0("conditions_from_meta",0)]]), split = "_")[[1]][1]
  
  # normalize
  data[data == 0] <- 1
  profiledata[profiledata == 0] <- 1
  sizefactor <- matrix(rep(1,nrow(data)), nrow = nrow(data)) %*% matrix(colSums(data), nrow = 1)
  data <- (data/sizefactor)*1000000
  sizefactor <- matrix(rep(1,nrow(profiledata)), nrow = nrow(profiledata)) %*% matrix(colSums(profiledata), nrow = 1)
  profiledata <- (profiledata/sizefactor)*1000000
  
  # remove low counts and log normalize
  max_count <- apply(data,1,max)
  data <- data[which(max_count > 10),, drop = FALSE]
  data <- log(data)
  profiledata <- log(profiledata)
  
  # prep data and get overlapping genes
  profiledata <- profiledata[degenes,]
  genes <- intersect(rownames(data), rownames(profiledata))
  data <- data[rownames(data) %in% genes,, drop = FALSE]
  profiledata <- profiledata[rownames(profiledata) %in% genes,]
  profileConds <- profilemetadata[[profileCondsName]]
  
  # joint metadata
  columns <- colnames(data)
  sample_column_ind <- which(apply(metadata, 2, function(x) sum(x %in% columns) == length(columns)))
  sample_id <- colnames(metadata)[sample_column_ind]
  metadata1 <- metadata[,c(sample_id,profileCondsName), drop = FALSE]
  colnames(metadata1) <- NULL
  columns <- colnames(profiledata)
  sample_column_ind <- which(apply(profilemetadata, 2, function(x) sum(x %in% columns) == length(columns)))
  sample_id <- colnames(profilemetadata)[sample_column_ind]
  metadata2 <- profilemetadata[,c(sample_id,profileCondsName), drop = FALSE]
  colnames(metadata2) <- NULL
  jointmetadata <- as.data.frame(rbind(cbind(as.matrix(metadata1), "Target"), 
                                       cbind(as.matrix(metadata2), "Reference")))
  colnames(jointmetadata) <- c(sample_id, profileCondsName, "Dataset")
  
  # CCA joint embedding reduction
  num.cc <- min(ncol(data), ncol(profiledata))-1
  data_scaled <- t(apply(data, 1, scale))
  profiledata_scaled <-  t(apply(profiledata, 1, scale))
  mat3 <- crossprod(x = data_scaled, y = profiledata_scaled)
  cca.svd <- irlba(A = mat3, nv = num.cc)
  cca.data <- rbind(cca.svd$u, cca.svd$v)
  colnames(x = cca.data) <- paste0("CC", 1:num.cc)
  ncca <- cumsum(cca.svd$d^2)/sum(cca.svd$d^2)
  ncca <- length(ncca) - sum(ncca > 0.8) + 1
  if(ncca < 2) ncca <- 2
  cca.data <- t(apply(cca.data, 1, function(x) x/sqrt(sum(x^2))))
  ccametadata <- data.frame(cca.data[,1:2], jointmetadata)
  targetcca.data <- cca.data[1:ncol(data_scaled),1:ncca]
  profcca.data <- cca.data[(ncol(data_scaled)+1):(ncol(data_scaled)+ncol(profiledata_scaled)),1:ncca]
  
  # get scores
  if(ScoreMethod=="Silhouette"){
    
    # datax_spear <- (cor(profiledata, data, method = "spearman") + 1)/2
    datax_spear <- dist2(profcca.data, targetcca.data)
    rownames(datax_spear) <- colnames(profiledata)
    colnames(datax_spear) <-  colnames(data)
    datax_spear <- t(datax_spear)
    cross_scores <- (external_silhouette(profileConds, datax_spear) + 1)/2
    cross_scores <- t(cross_scores)
    colnames(cross_scores) <- paste0("Cross_",colnames(cross_scores))
    
    # outlier scores
    background_spear <- as.matrix(dist(profcca.data))
    outlier_scores <- (outlier_silhouette(profileConds, datax_spear, background_spear) + 1)/2
    outlier_scores <- t(outlier_scores)
    colnames(outlier_scores) <- paste0("Inter_",colnames(outlier_scores))
    
  }
  # if(ScoreMethod=="NNLS-based"){
  #
  #   # Expression Profiles
  #   profiles <- aggregate(t(profiledata),list(profileConds), mean)
  #   profile_names <- profiles$Group.1
  #   profiles <- t(profiles[,-1])
  #
  #   # Deconvulate points based on expression profiles of conditions
  #   scores <- NULL
  #   for(j in 1:ncol(data)){
  #     samples_DCV <- nnls(profiles,data[,j])
  #     norm_coef <- samples_DCV$x/sum(samples_DCV$x)
  #     scores <- rbind(scores,norm_coef)
  #   }
  #   # change this later
  #   colnames(scores) <- profile_names
  #   rownames(scores) <- colnames(data)
  #   # scores <- melt(scores)
  # }
  
  return(list(cross_scores = cross_scores, outlier_scores = outlier_scores, ccametadata = ccametadata))
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
#' @export
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
#'     
#' @export
#' 
external_silhouette <- function(cluster = NULL, dist2 = NULL){
  if(is.null(cluster)) return(NULL)
  
  cls <- unique(cluster)
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
  aggdata <- aggregate(dist2[1,],list(cluster),mean)
  rownames(allsil) <- aggdata$Group.1
  return(allsil)
}

#' outlier_silhouette
#' 
#' Calculate outlier scores for non-symmetric distance matrices. 
#' Suitable for comparison of expression profiles across different datasets and experiments 
#' detecing outlying samples compared to reference bulk data sets
#'
#' @param cluster cluster labels of profile data
#' @param dist2 non-symmetric similarity matrix
#' @param dist symmmetric similarity matrix
#'
#' @examples
#'      x <- outlier_silhouette()
#'     
#' @export
#' 
outlier_silhouette <- function(cluster = NULL, dist2 = NULL, dist = NULL){
  if(is.null(cluster)) return(NULL)
  
  # background distribution for each cluster
  back_dist <- lapply(unique(cluster), function(x){
    cur_dist <- dist[,cluster == x]
    cur_dist[which.min(rowMeans(cur_dist)),]
  })
  names(back_dist) <- unique(cluster)
  
  # t.test for p-value
  alltest <- apply(dist2,1,function(x){
    result <- sapply(back_dist, function(y){
      1-t.test(x, y, var.equal = FALSE)$p.value 
    })
    return(result)
  })
  
  return(alltest)
}

#' getProfilingParameter
#' 
#' Prepare DE analysis parameters for the profiling data
#'
#' @param session session
#' @param input input 
#'
#' @examples
#'      x <- getProfilingParameter()
#'     
#' @export
#' 
getProfilingParameter <- function(session = NULL, input = NULL) {
  if (is.null(session)) return(NULL)
  
  params <- c("NULL","NULL")
  params[1] <- isolate(input[[paste0("scoremethod",1)]])
  params[2] <- strsplit(isolate(input[[paste0("conditions_from_meta",0)]]), split = "_")[[1]][1]
  
  return(list(params = params))
}

#' getPCAEmbedding
#'
#' get single cell RNA Embedding
#'
#' @param output output
#' @param summary summary output name
#' @param data single cell ExpressionSet Object
#' @param ident column in scRNA metadata to visualize
#'
#' @examples
#'       x <- getPCAEmbedding()
#'
#' @export
#'
getPCAEmbedding <- function (output = NULL, summary = NULL, data = NULL, ident = NULL)
{
  if (is.null(output))
    return(NULL)
  meta.data <- data
  meta.data$Group <- meta.data[,ident]
  PCs <- colnames(meta.data)[grepl("^PC[1-9]",colnames(meta.data))]
  var.perc <- sapply(PCs,function(x) as.numeric(strsplit(x, split = "_")[[1]][2]))
  colnames(meta.data)[grepl("^PC[1-9]",colnames(meta.data))] <- c("PC1","PC2")
  if(!is.null(data)){
    output[[summary]] <- renderPlot({
      ggplot(meta.data, ggplot2::aes(x = PC1, y = PC2, colour = Group), diag = "blank") + geom_point() + xlab(paste0("PC1 (%", var.perc[1],")")) + ylab(paste0("PC2 (%", var.perc[2],")"))
    })
  }
}

#' getCCAEmbedding
#'
#' get single cell RNA Embedding
#'
#' @param output output
#' @param summary summary output name
#' @param data single cell ExpressionSet Object
#' @param input input
#'
#' @examples
#'       x <- getCCAEmbedding()
#'
#' @export
#'
getCCAEmbedding <- function (output = NULL, summary = NULL, data = NULL, input = NULL)
{
  if (is.null(data))
    return(NULL)
  ident <- strsplit(isolate(input[[paste0("conditions_from_meta",0)]]), split = "_")[[1]][1]
  meta.data <- data
  meta.data$Group <- meta.data[,ident]
  meta.data$Highlight <- ifelse(meta.data$Dataset == "Target", meta.data[,3], "")
  if(!is.null(data)){
    output[[summary]] <- renderPlot({
      ggplot(meta.data, ggplot2::aes(x = CC1, y = CC2, colour = Group, label = Highlight), diag = "blank") + geom_point() + xlab(paste0("CC1")) + ylab(paste0("CC2")) + geom_point(alpha = 0.6) + geom_text()
    })
  }
}

#' simulateProfilingData
#'
#' simulates a reference and target data for testing Comparative Profiling
#'
#' @param compcodeR.param parameters for the compcodeR package, simulating reference and target data
#' @param overlappingDE the percentage of the number of overlapping DE genes for reference and target data
#' @param seed seed value for random number generation
#'
#' @import compcodeR
#' 
#' @examples
#'       x <- simulateProfilingData()
#'
#' @export
#'
simulateProfilingData <- function (compcodeR.param = NULL, overlappingDE = NULL, seed = 1)
{
  if (is.null(compcodeR.param))
    return(NULL)
  
  # set seed 
  set.seed(seed)
  
  # generate reference data with some number of DE genes
  cat("Generating Reference Data \n")
  reference_DE <- generateSyntheticData("reference_data", 
                                        n.vars = compcodeR.param$n.vars, 
                                        n.diffexp = compcodeR.param$n.diffexp, 
                                        samples.per.cond = compcodeR.param$samples.per.cond_reference,
                                        fraction.upregulated = compcodeR.param$fraction.upregulated)
  reference_metadata <- reference_DE@sample.annotations
  reference_data <- reference_DE@count.matrix
  reference_genemetadata <- reference_DE@variable.annotations
  reference_upDEgenes <- rownames(reference_genemetadata)[reference_genemetadata$upregulation == 1]
  reference_downDEgenes <- rownames(reference_genemetadata)[reference_genemetadata$downregulation == 1]
  
  # generate target data with some number of DE genes
  cat("Generating Target Data \n")
  target_DE <- generateSyntheticData("reference_data", 
                                     n.vars = compcodeR.param$n.vars, 
                                     n.diffexp = round(compcodeR.param$n.diffexp*overlappingDE), 
                                     samples.per.cond = compcodeR.param$samples.per.cond_target,
                                     fraction.upregulated = compcodeR.param$fraction.upregulated)
  target_metadata <- target_DE@sample.annotations
  target_data <- target_DE@count.matrix
  target_genemetadata <- target_DE@variable.annotations
  target_upDEgenes <- rownames(target_genemetadata)[target_genemetadata$upregulation == 1]
  target_downDEgenes <- rownames(target_genemetadata)[target_genemetadata$downregulation == 1]
  
  # # switch target up and down DE genes to match with DE genes of reference data
  cat("Matching Overlapping Genes \n")
  target_genes <- rownames(target_data)
  reference_genes <- rownames(reference_data)
  new_target_genes <- target_genes
  new_target_genes[target_genes %in% target_upDEgenes] <- reference_upDEgenes[1:length(target_upDEgenes)]
  new_target_genes[target_genes %in% target_downDEgenes] <- reference_downDEgenes[1:length(target_downDEgenes)]
  new_target_genes <- c(reference_upDEgenes[1:length(target_upDEgenes)], 
                        reference_downDEgenes[1:length(target_downDEgenes)],
                        target_genes[!target_genes %in% c(reference_upDEgenes[1:length(target_upDEgenes)],
                                                          reference_downDEgenes[1:length(target_downDEgenes)])])
  rownames(target_data) <- new_target_genes
  
  return(list(reference = list(data = reference_data, metadata = reference_metadata, 
                               featuredata = reference_genemetadata, DEgenes = c(reference_upDEgenes, reference_downDEgenes)),
              target = list(data = target_data, metadata = target_metadata, 
                            featuredata = target_genemetadata, DEgenes = c(target_upDEgenes, target_downDEgenes))))
}