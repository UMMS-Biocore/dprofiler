#' dprofilerdeanalysis
#'
#' Module to perform and visualize Iterative DE results.
#' 
#' @param input, input variables
#' @param output, output objects
#' @param session, session 
#' @param data, a matrix that includes expression values
#' @param columns, columns
#' @param conds, conditions
#' @param params, de parameters
#' @param parent_session parent session
#' 
#' @return DE panel 
#' @export
#'
#' @examples
#'     x <- dprofilerdeanalysis()
#'
dprofilerdeanalysis <- function(input = NULL, output = NULL, session = NULL, 
                                data = NULL, columns = NULL, conds = NULL, params = NULL, 
                                parent_session = NULL){
    if(is.null(data)) return(NULL)
    
    # Iterative DE Algorithm
    deres <- reactive(runIterDE(data, columns, conds, params, session))
    
    # Expression Profiles
    expression_profiles <- reactive(getExpressionProfiles(deres(), data, columns, conds))
    
    # Apply Filters for DE and Iter DE Results
    prepDat <- reactive({
        applyFiltersNew(addDataCols(data, deres()$DEResults, columns, conds), input)
    })
    
    iterprepDat <- reactive({
        remaining_columns <- (columns != deres()$cleaned_columns)
        if(length(remaining_columns) > 0){
            columns <- columns[remaining_columns]
            conds <- conds[remaining_columns]
        }
        applyFiltersNew(addDataCols(data, deres()$IterDEResults, columns, conds), input)
    })
    
    # # Choose Cell Types and top markers
    # output$deconvolute_genes <- renderUI({
    #     list(selectInput(session$ns("deconvolute_genes"), "", selected = c("Homogeneous Conditions"),
    #                      choices = c("Heterogeneous Conditions","Homogeneous Conditions")))
    # })
    
    # Create a reactive scoring table as well
    score <- reactiveVal()
    deconvolute_genes <- reactiveVal()
    
    # Observe for Tables and Plots
    observe({
        
        # switch to heterogeneity analysis page
        updateTabsetPanel(parent_session, "DEAnalysisBox", "heterogeneity")
        
        temp <- expression_profiles()
        
        # get scores and expression profiles
        score(getFinalScores(deres(), data, columns, conds, params, 
                             ManualDEgenes = input$manualgenes, TopStat = input$topstat))

        # # get deconvolution genes
        # deconvolute_genes({
        #     which_genes <- ifelse(is.null(input$deconvolute_genes), TRUE, input$deconvolute_genes)
        #     if(which_genes == "Homogeneous Conditions"){
        #         deres()$IterDEgenes
        #     } else{
        #         deres()$DEgenes
        #     }
        # })
        
        # prepare DE tables
        dat <-  prepDat()[prepDat()$Legend == input$legendradio,]
        dat2 <- removeCols(c("ID", "x", "y","Legend", "Size"), dat)
        iterdat <-  iterprepDat()[iterprepDat()$Legend == input$legendradio,]
        iterdat2 <- removeCols(c("ID", "x", "y","Legend", "Size"), iterdat)
        
        # DE Results
        getTableDetails(output, session, "DEResults", dat2, modal=FALSE)
        getTableDetails(output, session, "IterDEResults", iterdat2, modal = FALSE)

        # Membership Scores
        getIterDESummary(output, session, "HomogeneityVenn", "HomogeneitySummary", deres(), params)
        getScoreDetails(output, session, "HomogeneityScores", score()$DEscore, score()$IterDEscore)
        
        # download handler for DE genes
        getDEgenesDownloadButtons(output, session, deres()$DEgenes, deres()$IterDEgenes)
    })
    
    list(dat = prepDat, DEgenes = deres()$DEgenes, iterdat = iterprepDat, IterDEgenes = deres()$IterDEgenes,
         score = score, deconvolute_genes = NULL, cleaned_columns = deres()$cleaned_columns)
}

#' getDEResultsUI
#' 
#' Creates a panel to visualize DE results
#'
#' @param id, namespace id
#' @return panel
#' @examples
#'     x <- getDEResultsUI("batcheffect")
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
               tabPanel(title = "Profiling Results",
                        fluidRow(
                            shinydashboard::box(title = "Summary",
                                                solidHeader = T, status = "info",  width = 12, collapsible = TRUE,
                                                column(4,uiOutput(ns("HomogeneitySummary"))),
                                                column(4,htmlOutput(ns("HomogeneitySummaryIter"))),
                                                column(4,htmlOutput(ns("HomogeneitySummaryDE")))
                            ),
                            shinydashboard::box(title = "# of DE Genes",
                                                solidHeader = T, status = "info",  width = 6, collapsible = TRUE,
                                                plotOutput(ns("HomogeneityVenn")),
                                                column(12,
                                                       downloadButton(ns("downloadBeforeGenes"), label = "Initial DEgenes"),
                                                       downloadButton(ns("downloadOverlapGenes"), label = "Overlapping DE genes"),
                                                       downloadButton(ns("downloadAfterGenes"), label = "Final DE genes"),
                                                )
                            ),
                            shinydashboard::box(title = "Membership Scores",
                                                solidHeader = T, status = "info",  width = 6, collapsible = TRUE,
                                                plotlyOutput(ns("HomogeneityScores")),
                                                column(5,actionButtonDE("gotodeconvolute", "Go to Cellular Composition Analysis", 
                                                                        styleclass = "primary", style = 'margin-top:21px')),
                                                column(5,actionButtonDE("gotoprofile", "Go to Comparative Profiling", 
                                                                        styleclass = "primary", style = 'margin-top:21px')),
                                                #column(4,uiOutput(ns("deconvolute_genes")))
                                                
                            ),
                            uiOutput(ns("maininitialplot")),                                  
                            uiOutput(ns("mainoverlapplot")),
                            uiOutput(ns("mainfinalplot")),
                            uiOutput(ns("BarMainUI")),
                            uiOutput(ns("BoxMainUI")),
                            uiOutput(ns("heatmapUI"))
                        ),
                        value = "heterogeneity"
               ),
               tabPanel(title = "Impure (Heterogeneous) Conditions",
                        fluidRow(
                            shinydashboard::box(title = "Differentially Expressed Genes",
                                                solidHeader = T, status = "info",  width = 12, collapsible = TRUE,
                                                uiOutput(ns("IterDEResults"))
                                                
                            ),
                            uiOutput(ns("maindeplot"))
                        )               
                        ),
               tabPanel(title = "Pure (Homogeneous) Conditions",
                        fluidRow(
                            shinydashboard::box(title = "Differentially Expressed Genes",
                                                solidHeader = T, status = "info",  width = 12, collapsible = TRUE,
                                                uiOutput(ns("DEResults"))
                            ),
                            uiOutput(ns("mainiterdeplot"))
                        )
                        )
        )
    )
}

#' cutOffSelectionUI
#'
#' Gathers the cut off selection for DE analysis
#'
#' @param id, namespace id
#' @note \code{cutOffSelectionUI}
#' @return returns the left menu according to the selected tab;
#' @examples
#'     x <- cutOffSelectionUI("cutoff")
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
#' Gathers the cut off selection for Scoring 
#'
#' @param id, namespace id
#' @note \code{ScoreCutOffSelectionUI}
#' @return returns the left menu according to the selected tab;
#' @examples
#'     x <- ScoreCutOffSelectionUI("cutoff")
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
#' @return
#' @export
#'
#' @examples
getDEgenesDownloadButtons <- function(output, session,  DEgenes, IterDEgenes){
    
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

runDE <- function (data = NULL, columns = NULL, conds = NULL, params = NULL) 
{
    if (is.null(data)) 
        return(NULL)
    de_res <- NULL
    if (startsWith(params[1], "DESeq2")) 
        de_res <- runDESeq2(data, columns, conds, params)
    else if (startsWith(params[1], "EdgeR")) 
        de_res <- runEdgeR(data, columns, conds, params)
    else if (startsWith(params[1], "Limma")) 
        de_res <- runLimma(data, columns, conds, params)
    data.frame(de_res, gene = rownames(de_res))
}

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

runLimma <- function (data = NULL, columns = NULL, conds = NULL, params = NULL) 
{
    if (is.null(data)) 
        return(NULL)
    if (length(params) < 3) 
        params <- strsplit(params, ",")[[1]]
    normfact = if (!is.null(params[2])) 
        params[2]
    fitType = if (!is.null(params[3])) 
        params[3]
    normBet = if (!is.null(params[4])) 
        params[4]
    datatype = if(!is.null((params[5])))
        params[5]
    data <- data[, columns]
    conds <- factor(conds)
    cnum = summary(conds)[levels(conds)[1]]
    tnum = summary(conds)[levels(conds)[2]]
    filtd <- data
    des <- factor(c(rep(levels(conds)[1], cnum), rep(levels(conds)[2], 
                                                     tnum)))
    names(filtd) <- des
    design <- cbind(Grp1 = 1, Grp2vs1 = des)
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


