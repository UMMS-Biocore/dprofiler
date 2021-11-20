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
#' @examples
#'     x <- dprofilerdeanalysis()
#'     
#' @export
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
                                                column(12, 
                                                       p(strong("Note:")," Differentially expressed genes and ", strong("Membership Scores")," of samples are calculated to iteratively remove samples with low membership scores. Press ", 
                                                         strong("Start"), " to score lesional and non-lesional samples which results in", strong("P65_NL being removed."), " Here, P65_NL is a non-lesional Vitiligo sample with a low score, suggesting 
                                                         that its expression profile may ", strong("not match with the phenotypic profile of a non-lesional Vitiligo sample."))
                                                ),
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
                                                column(5,actionButtonDE("gotodeconvolute", "Go to Cellular Composition Analysis", 
                                                                        styleclass = "primary", style = 'margin-top:21px')),
                                                column(5,actionButtonDE("gotoprofile", "Go to Comparative Profiling", 
                                                                        styleclass = "primary", style = 'margin-top:21px')),
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
               tabPanel(title = "Results before Profiling",
                        fluidRow(
                            shinydashboard::box(title = "Differentially Expressed Genes",
                                                solidHeader = T, status = "info",  width = 12, collapsible = TRUE,
                                                uiOutput(ns("IterDEResults"))
                                                
                            ),
                            uiOutput(ns("maindeplot"))
                        )               
                        ),
               tabPanel(title = "Results before Profiling",
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
#' @param columns a vector that includes the columns that are going to be analyzed. These columns has to match with the given data.
#' @param conds experimental conditions. The order has to match with the column order
#' @param params all params for the DE methods
#'
#' @examples
#'      x <- runDE()
#'     
#' @export
#'     
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


