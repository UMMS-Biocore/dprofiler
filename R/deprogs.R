#' dprofilerdeanalysis
#'
#' Module to perform and visualize Iterative DE results.
#' 
#' @param input, input variables
#' @param output, output objects
#' @param session, session 
#' @param parent_session parent session
#' @param cond_dc data and DE parameters
#'
#' @examples
#'     x <- dprofilerdeanalysis()
#'     
#' @export
#' 
dprofilerdeanalysis <- function(input = NULL, output = NULL, session = NULL, parent_session = NULL, cond_dc = NULL){
  
  ###
  # Set Data ####
  ### 
  
  CountData <- cond_dc$selectedData
  DEparameters <- cond_dc$DEparameters
  startDE <- cond_dc$startDE
  
  ###
  # Reactive Values ####
  ### 
  
  CrossScore <- reactiveVal()
  IntraScore <- reactiveVal()
  ScoreTable <- reactiveVal()
  Summarise <- reactive(input$summarisescores)
  Deconvolute <- reactive(input$gotodeconvolute)
  Profile <- reactive(input$gotoprofile)
   
  ###
  # DE Analysis ####
  ### 
  
  deres <- reactive(runIterDE(CountData()$count, CountData()$meta, DEparameters()$cols, 
                              DEparameters()$conds, DEparameters()$params, session))
  
  # Iterative DE Algorithm
  observeEvent(cond_dc$startDE(), {
    waiter_show(html = spin_ring(), color = transparent(.5))
    withProgress(message = 'Running Computational Profiling',
                 detail = paste0("ScoringMethod: ", DEparameters()$params[7], " DEmethod: ", DEparameters()$params[1]),
                 value = 0, {
                   deres()
                   incProgress(1)
                 })
    waiter_hide()
    
    # hide initial tabs
    showTab("menutabs","discover", session = parent_session)
    showTab("DEAnalysisBox","CompPhenoProfiling", session = parent_session)
    showTab("DEAnalysisBox","ResultsAfterCompPhenoProfiling", session = parent_session)
    showTab("DEAnalysisBox","ResultsBeforeCompPhenoProfiling", session = parent_session)
    updateTabsetPanel(parent_session, "DEAnalysisBox", "CompPhenoProfiling")
  })
  
  # Apply Filters for DE
  prepDat <- reactive({
    if(!is.null(deres())){
      applyFiltersNew(addDataCols(CountData()$count, deres()$DEResults,
                                  DEparameters()$cols, DEparameters()$conds), input)
    }
  })

  # Apply Iter DE Results
  iterprepDat <- reactive({
    if(!is.null(deres())){
      remaining_columns <- deres()$remaining_columns
      conds <- DEparameters()$cols[DEparameters()$cols %in% remaining_columns]
      columns <- DEparameters()$cols[DEparameters()$cols %in% remaining_columns]
      applyFiltersNew(addDataCols(CountData()$count, deres()$IterDEResults,
                                  columns, conds), input) 
    }
  })

  ###
  # Window Settings ####
  ###
  
  # hide initial tabs
  hideTab("DEAnalysisBox","CompPhenoProfiling", session = parent_session)
  hideTab("DEAnalysisBox","ResultsAfterCompPhenoProfiling", session = parent_session)
  hideTab("DEAnalysisBox","ResultsBeforeCompPhenoProfiling", session = parent_session)
  
  ###
  # Main Observable ####
  ###

  # Observe for Tables and Plots
  observe({

    # get cross condition scores and expression profiles
    CrossScore(getFinalScores(deres(), CountData()$count, DEparameters()$cols, DEparameters()$conds, DEparameters()$params,
                              ManualDEgenes = input$manualgenes, TopStat = input$topstat))

    # get intra condition scores and expression profiles
    IntraScore(getOutlierScores(CountData()$count, DEparameters()$cols, DEparameters()$conds, DEparameters()$params))

    # merge cross and intra scores
    ScoreTable(MergeScoreTables(CrossScore()$IterDEscore, IntraScore()$Score))
  
    # prepare DE tables
    dat <-  prepDat()[prepDat()$Legend == input$legendradio,]
    dat2 <- removeCols(c("ID", "x", "y","Legend", "Size"), dat)
    iterdat <-  iterprepDat()[iterprepDat()$Legend == input$legendradio,]
    iterdat2 <- removeCols(c("ID", "x", "y","Legend", "Size"), iterdat)

    # DE Results
    getTableDetails(output, session, "DEResults", dat2, modal=FALSE)
    getTableDetails(output, session, "IterDEResults", iterdat2, modal = FALSE)

    # Membership Scores
    getIterDESummary(output, session, "HomogeneityVenn", "HomogeneitySummary", deres(), DEparameters()$params)
    getScoreDetails(output, session, "HomogeneityScores", CrossScore(), IntraScore())

    # download handler for DE genes
    getDEgenesDownloadButtons(output, session, deres()$DEgenes, deres()$IterDEgenes)
  })
  
  list(count = prepDat, count_iter = iterprepDat, CountData = CountData, DEresults = deres, CrossScore = CrossScore, IntraScore = IntraScore, 
       ScoreTable = ScoreTable, Summarise = Summarise, Deconvolute = Deconvolute, Profile = Profile)
}

#' getDEResultsUI
#' 
#' Creates a panel to visualize DE results
#'
#' @param id, namespace id
#' 
#' @examples
#'     x <- getDEResultsUI("deresults")
#'     
#' @export
#' 
getDEResultsUI<- function (id) {
  ns <- NS(id)
  list(
    tabBox(id = "DEAnalysisBox",
           width = NULL,
           tabPanel(title = "Conditions",
                    condSelectUI()
           ),
           tabPanel(title = "Profiling Summary",
                    fluidRow(
                      shinydashboard::box(title = "Summary",
                                          solidHeader = T, status = "info",  width = 12, collapsible = TRUE,
                                          #column(12, 
                                          #       p(strong("Note:")," Differentially expressed genes and ", strong("Membership Scores")," of samples are calculated to iteratively remove samples with low membership scores. Press ", 
                                          #         strong("Start"), " to score lesional and non-lesional samples which results in", strong("P65_NL being removed."), " Here, P65_NL is a non-lesional Vitiligo sample with a low score, suggesting 
                                          #         that its expression profile may ", strong("not match with the phenotypic profile of a non-lesional Vitiligo sample."))
                                          #),
                                          column(4,uiOutput(ns("HomogeneitySummary"))),
                                          column(4,htmlOutput(ns("HomogeneitySummaryIter"))),
                                          column(4,htmlOutput(ns("HomogeneitySummaryDE")))
                      ),
                      shinydashboard::box(title = "# of DE Genes",
                                          solidHeader = T, status = "info",  width = 6, collapsible = TRUE,
                                          plotOutput(ns("HomogeneityVenn")),
                                          column(12,
                                                 downloadButton(ns("downloadBeforeGenes"), label = "DEgenes before Prof."),
                                                 downloadButton(ns("downloadOverlapGenes"), label = "DE genes (Overlapping)"),
                                                 downloadButton(ns("downloadAfterGenes"), label = "DE genes after Prof."),
                                          )
                      ),
                      shinydashboard::box(title = "Membership Scores",
                                          solidHeader = T, status = "info",  width = 6, collapsible = TRUE,
                                          plotlyOutput(ns("HomogeneityScores")),
                                          column(4,actionButtonDE(ns("summarisescores"), "Summarise Scores", 
                                                                  styleclass = "primary", style = 'margin-top:21px')),
                                          column(4,actionButtonDE(ns("gotodeconvolute"), "Go to Compositional Profiling", 
                                                                  styleclass = "primary", style = 'margin-top:21px')),
                                          column(4,actionButtonDE(ns("gotoprofile"), "Go to Comparative Profiling", 
                                                                  styleclass = "primary", style = 'margin-top:21px')),
                      ),
                      uiOutput(ns("maininitialplot")),                                  
                      uiOutput(ns("mainoverlapplot")),
                      uiOutput(ns("mainfinalplot")),
                      uiOutput(ns("BarMainUI")),
                      uiOutput(ns("BoxMainUI")),
                      uiOutput(ns("heatmapUI"))
                    ),
                    value = "CompPhenoProfiling"
           ),
           tabPanel(title = "Results After Profiling",
                    fluidRow(
                      shinydashboard::box(title = "Differentially Expressed Genes",
                                          solidHeader = T, status = "info",  width = 12, collapsible = TRUE,
                                          uiOutput(ns("IterDEResults"))
                                          
                      ),
                      uiOutput(ns("maindeplot"))
                    ),
                    value = "ResultsAfterCompPhenoProfiling"
           ),
           tabPanel(title = "Results before Profiling",
                    fluidRow(
                      shinydashboard::box(title = "Differentially Expressed Genes",
                                          solidHeader = T, status = "info",  width = 12, collapsible = TRUE,
                                          uiOutput(ns("DEResults"))
                      ),
                      uiOutput(ns("mainiterdeplot"))
                    ),
                    value = "ResultsBeforeCompPhenoProfiling"
           )
    )
  )
}

#' cutOffSelectionUI
#'
#' Gathers the cut off selection for DE analysis
#'
#' @param id, namespace id
#' 
#' @examples
#'     x <- cutOffSelectionUI("cutoff")
#'     
#' @export
#' 
cutOffSelectionUI <- function(id){
  ns <- NS(id)
  list(
    getLegendRadio(id),
    textInput(ns("padj"), "padj value cut off", value = "0.01" ),
    textInput(ns("foldChange"), "foldChange", value = "2" )
  )
}

#' ScoreCutOffSelectionUI
#'
#' Gathers the cut off selection for scoring 
#'
#' @param id namespace id
#' 
#' @examples
#'     x <- ScoreCutOffSelectionUI("cutoff")
#'     
#' @export
#'   
ScoreCutOffSelectionUI <- function(id){
  ns <- NS(id)
  list(
    textInput(ns("topstat"), "Top Stat", value = "" ),
    fileInput(ns("manualgenes"), "Manual DEgenes")
  )
}

#' getDEgenesDownloadButtons
#'
#' Buttons for downloading Initial, overlapping and DE genes
#' 
#' @param output output 
#' @param session session
#' @param DEgenes Initial DE genes
#' @param IterDEgenes Final DE genes
#'
#' @examples
#'     x <- getDEgenesDownloadButtons()
#'     
#' @export
#'    
getDEgenesDownloadButtons <- function(output = NULL, session = NULL,  DEgenes = NULL, IterDEgenes = NULL){
  if(is.null(DEgenes)) return(NULL)
  
  genes <- setdiff(DEgenes,IterDEgenes)
  if(length(genes) == 0) genes <- DEgenes
  output$downloadBeforeGenes <- downloadHandler(
    filename = function() {paste('initial_degenes.txt')},
    content = function(con) {write(genes, con)}
  )
  
  genes <- intersect(DEgenes,IterDEgenes)
  if(length(genes) == 0) genes <- IterDEgenes
  output$downloadOverlapGenes <- downloadHandler(
    filename = function() {paste('overlapping_degenes.txt')},
    content = function(con) {write(genes, con)}
  )
  
  genes <- setdiff(IterDEgenes,DEgenes)
  if(length(genes) == 0) genes <- IterDEgenes
  output$downloadAfterGenes <- downloadHandler(
    filename = function() { paste('final_degenes.txt')},
    content = function(con) {write(genes, con)}
  )
}

#' runDE
#' 
#' Run DE algorithms on the selected parameters. Output is to be used for the interactive display.
#' Adapted from debrowser::runDE
#'
#' @param data A matrix that includes all the expression raw counts, rownames has to be the gene, isoform or region names/IDs.
#' @param metadata metadata
#' @param columns a vector that includes the columns that are going to be analyzed. These columns has to match with the given data.
#' @param conds experimental conditions. The order has to match with the column order
#' @param params all params for the DE methods
#'
#' @examples
#'      x <- runDE()
#'     
#' @export
#'     
runDE <- function (data = NULL, metadata = NULL, columns = NULL, conds = NULL, params = NULL) 
{
  if (is.null(data)) 
    return(NULL)
  de_res <- NULL
  if (startsWith(params[1], "DESeq2")) 
    de_res <- runDESeq2(data, metadata, columns, conds, params)
  else if (startsWith(params[1], "EdgeR")) 
    de_res <- runEdgeR(data, metadata, columns, conds, params)
  else if (startsWith(params[1], "Limma")) 
    de_res <- runLimma(data, metadata, columns, conds, params)
  data.frame(de_res, gene = rownames(de_res))
}

#' runMultipleDE
#' 
#' Run DE algorithms on the selected parameters. Output is to be used for the interactive display.
#' Multiple comparisons version of runDE. 
#'
#' @param data A matrix that includes all the expression raw counts, rownames has to be the gene, isoform or region names/IDs.
#' @param columns a vector that includes the columns that are going to be analyzed. These columns has to match with the given data.
#' @param conds Multiple experimental conditions. 
#' @param params all params for the DE methods
#'
#' @examples
#'      x <- runMultipleDE()
#'     
#' @export
#'     
runMultipleDE <- function (data = NULL, columns = NULL, conds = NULL, params = NULL) 
{
  if (is.null(data)) 
    return(NULL)
  
  combn_treatment <- combn(levels(droplevels(conds)),2)
  multiplede_res <- NULL
  
  for(i in 1:ncol(combn_treatment)){
    comb_treat <- combn_treatment[,i]
    
    cur_conds <- conds[conds %in% comb_treat]
    cur_columns <- columns[conds %in% comb_treat]
    cur_data <- data[,conds %in% comb_treat]
    
    de_res <- NULL
    if (startsWith(params[1], "DESeq2")) 
      de_res <- runDESeq2(cur_data, cur_columns, cur_conds, params)
    else if (startsWith(params[1], "EdgeR")) 
      de_res <- runEdgeR(cur_data, cur_columns, cur_conds, params)
    else if (startsWith(params[1], "Limma")) 
      de_res <- runLimma(cur_data, cur_columns, cur_conds, params)
    data.frame(de_res)
    
    multiplede_res <- rbind(multiplede_res, 
                            data.frame(de_res, 
                                       Comparison = paste0(comb_treat[1], "_vs_", comb_treat[2]),
                                       gene = rownames(de_res)))
  }
  return(multiplede_res)
}

#' runLimma
#' 
#' Run Limma algorithm on the selected conditions. Output is to be used for the interactive display.
#' Adapted from debrowser::runLimma(), a version for non-count RNA data (microarray)
#' 
#' @param data A matrix that includes all the expression raw counts, rownames has to be the gene, isoform or region names/IDs.
#' @param columns is a vector that includes the columns that are going to be analyzed. These columns has to match with the given data.
#' @param conds experimental conditions. The order has to match with the column order
#' @param params normfact: Calculate normalization factors to scale the raw library sizes. Values can be "TMM","RLE","upperquartile","none". fitType, fitting method; "ls" for least squares or "robust" for robust regression normBet: Normalizes expression intensities so that the intensities or log-ratios have similar distributions across a set of arrays. datatype: suggest if the data is a count data or not. 
#'
#' @examples
#'      x <- runLimma()
#'     
#' @export
#'     
runLimma <- function (data = NULL, columns = NULL, conds = NULL, params = NULL) 
{
  if (is.null(data)) 
    return(NULL)
  if (length(params) < 3) 
    params <- strsplit(params, ",")[[1]]
  covariates <- if (!is.null(params[2])) params[2]
  covariates <- strsplit(covariates, split = "\\|")[[1]]
  normfact = if (!is.null(params[3])) 
    params[3]
  fitType = if (!is.null(params[4])) 
    params[4]
  normBet = if (!is.null(params[5])) 
    params[5]
  datatype = if(!is.null((params[6])))
    params[6]
  data <- data[, columns]
  conds <- factor(conds)
  cnum = summary(conds)[levels(conds)[1]]
  tnum = summary(conds)[levels(conds)[2]]
  filtd <- data
  des <- factor(c(rep(levels(conds)[1], cnum), rep(levels(conds)[2], 
                                                   tnum)))
  names(filtd) <- des
  
  if(covariates != "NoCovariate"){
    design <- cbind(Grp1=1,Grp2vs1=des)
    cov_metadata <- metadata[match(columns,metadata$sample),covariates, drop = FALSE]
    for(i in 1:length(covariates)){
      cur_covariate <- cov_metadata[,i]
      cur_covariate <- factor(cur_covariate)
      design <- cbind(design, cur_covariate)
      colnames(design)[length(colnames(design))] <-  paste0("covariate",i)
    }
  } else {
    design <- cbind(Grp1=1,Grp2vs1=des)
  }
  
  if(datatype == "count"){
    dge <- DGEList(counts = filtd, group = des)
    dge <- calcNormFactors(dge, method = normfact, samples = columns)
    v <- voom(dge, design = design, normalize.method = normBet,
              plot = FALSE)
    fit <- lmFit(v, design = design)
  } else {
    fit <- lmFit(filtd, design = design)
  }
  fit <- eBayes(fit)
  options(digits = 4)
  tab <- topTable(fit, coef = 2, number = dim(fit)[1], genelist = fit$genes$NAME)
  res <- data.frame(cbind(tab$logFC, tab$P.Value, tab$adj.P.Val, tab$t))
  colnames(res) <- c("log2FoldChange", "pvalue", "padj", "stat")
  rownames(res) <- rownames(tab)
  return(res)
}


