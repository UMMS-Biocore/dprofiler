#' getIntroText
#' 
#' Intro text
#'
#' @return the JS for tab updates
#'
#' @examples
#'     x<- getIntroText()
#'
#' @export
getIntroText<-function(){
  
  list(
    column(9,
           
      ####
      ## Introduction ####
      ####
           
      h2("Quick Start Guide"),
      p("Differential expression (DE) analysis has become an increasingly popular tool
  in determining and viewing up and/or down experssed genes between two sets of
  samples."),
      p("However, there might exist some confounding samples within each condition of a differential 
        expression analysis that hinders the significance of a function of a gene/protein"),
      p("The goal of Dprofiler is to detect this confounding elements in the data set by:"),
      tags$ul(
        tags$li(strong("Aim 1:"), "Detecting and Scoring samples that consistitute some heterogeneity within each condition"),
        tags$li(strong("Aim 2:"), "Validating the cellular composition of each scored sample with a reference scRNA data and estimate fractions of cellular compositions"),
        tags$li(strong("Aim 3:"), "Scoring third party data sets with condition-specific expression profiles from Aim 1"),
      ),
      
      ####
      ## Data input section ####
      ####
      
      h3("1. Data input"),
      p("There are three types input data in Dprofiler. These are: "),
      tags$ul(
        tags$li(strong("Reference Bulk Data Set"), ", used for detecting heterogeneous samples within and establish homogeneous conditions for scoring Profiling data sets"),
        tags$li(strong("scRNA Data Set"), ", used for validating scores of samples in Reference data set and for conducting cellular RNA deconvolution analysis"),
        tags$li(strong("Profiled Bulk Data Set"), "used for scoring its samples with the condition-specific expression profiles of Reference Bulk Data Set"),
      ),
      p("For both Reference and Profiling Bulk Data sets, there are two distinct inputs given in '.txt', '.csv' or '.tsv' formats.
        These are Count Data File and Metadata File where inputing metadata is optional"),
      p("However, for scRNA data set, the user should provide an ExpressionSet object"),
      p("If you do not have a dataset to upload, you can use the built in demo data file by 
    clicking on the 'Load Demo' button that loads a case study. To view the entire demo data file, you can 
    download"),
      tags$ol(
        tags$li(a("PRJNA554241",href="https://www.ncbi.nlm.nih.gov/bioproject/PRJNA554241"),": a reference bulk RNA-Seq count data of lesional and non-lesional vitiligo skin samples 
        processed by the Kallisto pipeline of ", a("DolphinNext", href = "https://github.com/UMMS-Biocore/dolphinnext")),
        tags$li(a("GSE65127",href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65127"),": a profiling bulk microarray data of vitiligo skin samples")
      ),
      h4("1.1 Count Data File"),            
      p("This input file should contain summarized count results of all samples in the experiment, 
    an example of the expected input data format is presented as below:"),
      tableOutput("countFile"),
      p("Where columns are samples, rows are the mapped genomic features 
    (e.g. genes, isoforms, miRNAs, ATAC or Chip regions etc.)."),
      h4("1.2 Metadata File"),
      p("In addition to the count data file; you may also upload metadata file to correct for 
  batch effects or any other normalizing conditions you might want to address that might be
  within your results. To handle for these conditions, simply create a metadata file by 
  using the example table at below "),
      tableOutput("metaFile"),
      p("To be able to upload the data please press ", strong("Upload"), ". "),
      p("After sucessfully uploading your data, you should see the summary of your data in the 'Upload Summary' section. 
        To move the filtering section please click ",  strong("Filter"), ". "),
      p("If you are ready to use Dprofiler, please click 'Upload' menu on the left to start using Dprofiler"),
    h4("1.3 scRNA Data File"),
    p("However, for scRNA data set, the user should provide an ExpressionSet object with metadata data slot 
        (see pData) that has three metadata variables available for each cell:", strong(" (i) "), "Samples",  
      strong(" (ii) "), "UMI counts",  strong(" (iii) "), "CellType")
    )
  )

}

#' getDataAssesmentText
#' 
#' DataAssesment text
#'
#' @return help text for data assesment 
#'
#' @examples
#'     x<- getDataAssesmentText()
#'
#' @export
getDataAssesmentText<-function(){
  list(
    column(9,
      h3("2. Data Assesment"),
      h4("2.1 Filtering"),    
      p("In this section, you can simultaneously visualise the changes of your dataset 
      while filtering out the low count genes. Choose your filtration criteria from ", 
      strong("Filtering Methods"), " box which is located just center of the screen. 
      Three methods are available to be used:"),
      tags$ul(
        tags$li(strong("Max:"), "Filters out genes where maximum count for each gene across all samples are less 
     than defined threshold."),
        tags$li(strong("Mean:"), "Filters out genes where mean count for each gene are less than defined threshold."),
        tags$li(strong("CPM:"), "First, counts per million (CPM) is calculated as the raw counts divided by the 
      library sizes and multiplied by one million. Then it filters out genes where at least 
      defined number of samples is less than defined CPM threshold.")
      ),
      withMathJax(),
      p("The expression cutoff value is determined according to the library size 
      and normalization factors with formula $$\\text{CPM} = \\frac{\\text{raw counts}}{\\text{library size} * \\text{normalization factors}} * 10^{-6}$$ 
      For example, if the cutoff CPM value is 10,
      the library size and normalization factors are estimated approximately equal to \\(\\ 3 \\text{ x} 10 ^ 6\\) and 1 for at least 4 samples, 
      then 10 CPM expression cutoff corresponds to about 30 read counts. 
      Therefore, in this example features in more than 4 samples have less than 
      30 read counts (10 CPM) is going to be low expression features
      and will be removed for batch effect correction and DE analysis."),
      p("To be able to filter out the low expression counts please press", strong("Filter"), ". "),
      h4("2.2 Quality control(QC)"),
      p("After filtering low count features, you may continue your analysis with Batch Effect 
      Detection & Correction or directly jump to differential expression analysis or view 
      quality control (QC) information of your dataset."),
      p("If specified metadata file containing your treatment and batch fields, by clicking ", strong("Batch effect correction"),
        " button, you have the option to conduct:"), 
      tags$ul(
        tags$li("Principal Components Analysis (PCA)"),
        tags$li("Interquartile range (IQR)"),
        tags$li("Density plots")
      ),
      p("to assess if the data requires batch effect correction or not."),
      p("If user wants to skip batch effect assesment and correction step, they can either click:"), 
      tags$ul(
        tags$li(strong("Go to DE Analysis"), " button to perform DE Analysis or "),
        tags$li(strong("Go to QC plots")," button for QC plots to draw PCA, all2all scatter, heatmaps, IQR and density plots.")
      ),
      h4("2.3. Data Preparation"),
      p("With metadata file containing your batch correction fields 
      then you have the option to conduct ",strong("batch effect correction"), " prior to 
      your analysis. By adjusting parameters of ", strong("Options") ," box, you can investigate 
      your dataset further. These parameters of the ", strong("Options"), " box are 
      explained as following:"),
      tags$ol(
        tags$li(strong("Normalization Method:"), "Dprofiler allows performing normalization 
      prior the batch effect correction. You may choose your normalization method 
      (among MRE, TMM, RLE, upperquartile), or if you don't want to normalize your 
      data you can select none for this item."),
        tags$li(
          strong("Correction Method:"), "Dprofiler uses ",
          a("ComBat", href = "https://bioconductor.org/packages/release/bioc/html/sva.html"),
          " (part of the SVA bioconductor package) or ",
          a("Harman", href = "https://www.bioconductor.org/packages/release/bioc/vignettes/Harman/inst/doc/IntroductionToHarman.html"),
          "to adjust for possible batch effect or conditional biases."),
        tags$li(strong("Treatment:"), "Please select the column that is specified in 
      metadata file for comparision, such as cancer vs control. It is named 
      condition for our sample metadata."),
        tags$li(strong("Batch:"), "Please select the column name in metadata file 
      which differentiate the batches. For example in our metadata, it is called batch.
      Upon clicking submit button, comparison tables and plots will be created on the right 
      part of the screen as shown below."),
      ),
      p("You can investigate the changes on the data by comparing following features:"),
      tags$ul(
        tags$li("Read counts for each sample."),
        tags$li("PCA, IQR and Density plot of the dataset."),
        tags$li("Gene/region vs samples data.")
      ),
      p("After batch effect correction, user can click ", strong("Go to DE Analysis"),
        " button to perform DE Analysis or ", strong("Go to QC plots"),
        " button for QC plots to draw PCA, all2all scatter, heatmaps, IQR and density plots.")
    )
  )
}

#' getDEAnalysisText
#' 
#' DEAnalysis text
#'
#' @return help text for DE Analysis 
#'
#' @examples
#'     x<- getDEAnalysisText()
#'
#' @export
getHeteroAnalysisText<-function(){
  list(
    column(9,
      h3("3. Differential Heterogeneity analysis"),
      p("Dprofiler provides methods for detecting heterogeneity within data sets that confound the detection of genes that are differentially expressed"),
      
      h4("3.1 An Iterative DE analysis and Membership Scores"),
      
      p("Dprofiler iteratively sweeps the data sets to assess the adequecy of each samples to associated conditions and on each iteration:"),
      tags$ul(
        tags$li("Conducts a DE analysis (DESeq2, EdgeR or Limma) given remaining samples within the data"),
        tags$li("Estimates the membership scores (Silhouette and NNLS) of all samples given expression profiles limited to differentially expressed genes"),
        tags$li("Removes samples with low membership scores"),
        tags$li("Repeats until no more samples have to removed from the data")
      ),
      p("Parameters that are used to conduct iterative DE analysis are:"),
      tags$ul(
        tags$li(strong("Score Method:"), " The algorithm to determine membership scores of each individual sample:",
                tags$ul(
                  tags$li(strong("Silhouette:"), "predicts the membership score via the normalized difference of 
                          average distances to one condition and average distance to the other condition"),
                  tags$li(strong("NNLS"), "estimates score with non-negative least squares regression where each 
                          samples modeled against mean expression profiles of each condition")
                )
        ),
        tags$li(strong("Min. Score:"), " A threshold for membership scores to determine dismembered samples of the dataset,
                samples whose scores are smaller than the threshold are eliminated on each iteration (However, their score are calculated
                n the last iteration)",
        ),
        tags$li(strong("DE Selection Method:"), " The protocol for selecting DE genes on each iteration",
                tags$ul(
                  tags$li(strong("log2FC+Padj:"), "selects DE genes whose log fold changes are higher and adjusted p-values are lower 
                          than a predetermined value"),
                  tags$li(strong("Top n Stat:"), "selects DE genes whose test statistics are among the top n genes with highest statistics"),
                )
        )
      ),
      p("Users who wish to analyze the data further should click ", strong("Go to Cellular Compositions"), " button to initiate a RNA deconvolution using 
        the reference scRNA data set"),
      
      h4("3.2 DE analysis"),
      p("The goal of differential gene expression analysis is to find genes
  or transcripts whose difference in expression, when accounting for the
  variance within condition, is higher than expected by chance."),
      tags$ul(
        tags$li(
          a("DESeq2", href="https://bioconductor.org/packages/release/bioc/html/DESeq2.html"),
          "is an R package available via Bioconductor and is designed to normalize count 
          data from high-throughput sequencing assays such as RNA-Seq and test for
          differential expression (Love et al. 2014).  With multiple parameters such as
          padjust values, log fold changes, plot styles, and so on, altering plots
          created with your DE data can be a hassle as well as time consuming. The
          Differential Expression Browser uses DESeq2 (Love et al., 2014)"
        ),
        tags$li(
          a("EdgeR",href="https://bioconductor.org/packages/release/bioc/html/edgeR.html"),
          "(Robinson et al., 2010), and ",
          a("Limma", href="https://bioconductor.org/packages/release/bioc/html/limma.html"),
          "(Ritchie et al., 2015) coupled with shiny (Chang, W. et al., 2016)
          to produce real-time changes within your
          plot queries and allows for interactive browsing of your DE results.
          In addition to DE analysis, Dprofiler also offers a variety of other plots
          and analysis tools to help visualize your data even further."),
      ),
      # p("If you are ready to discover and explore your data, please click ", strong("Go to Main Plot"), 
      #   " button in DE Results section."),
      p("Parameters that are used to conduct DE analysis are:"),
      tags$ul(
        tags$li(strong("DESeq2"),
                tags$ol(
                  tags$li(strong("fitType:"), " Either 'parametric', 'local', or 'mean' for the type of fitting of 
                          dispersions to the mean intensity. See estimateDispersions for description."),
                  tags$li(strong("betaPrior:"), " Whether or not to put a zero-mean normal prior on the non-intercept 
                          coefficients See nbinomWaldTest for description of the calculation of 
                          the beta prior. By default, the beta prior is used only for the Wald test, 
                          but can also be specified for the likelihood ratio test"),
                  tags$li(strong("testType:"), " Either 'Wald' or 'LRT', which will then use either Wald significance tests 
                          (defined by nbinomWaldTest), or the likelihood ratio test on the difference in 
                          deviance between a full and reduced model formula (defined by nbinomLRT)")
                )
          
        )
      ),
      tags$ul(
        tags$li(strong("EdgeR"),
                tags$ol(
                  tags$li(strong("Normalization:"), " Calculate normalization factors to scale the raw library sizes. Values 
                          can be 'TMM','RLE','upperquartile','none'."),
                  tags$li(strong("Dispersion:"), " Either a numeric vector of dispersions or a character string indicating 
                          that dispersions should be taken from the data object."),
                  tags$li(strong("testType:"), "ExactTest or glmLRT.",
                          tags$ul(
                            tags$li(strong("exactTest:")," Computes p-values for differential 
                                    abundance for each gene between two samples, conditioning on the total 
                                    count for each gene. The counts in each group are assumed to follow a 
                                    binomial distribution. "),
                            tags$li(strong("glmLRT:")," Fits a negative binomial generalized 
                                    log-linear model to the read counts for each gene and conducts 
                                    genewise statistical tests.")
                          )
                  ),
                )
                
        )
      ),
      tags$ul(
        tags$li(strong("Limma"),
                tags$ol(
                  tags$li(strong("Normalization:"), " Calculate normalization factors to scale the raw library sizes. Values 
                          can be 'TMM','RLE','upperquartile','none'."),
                  tags$li(strong("Fit Type:"), " fitting method; 'ls' for least squares or 'robust' for robust regression"),
                  tags$li(strong("Norm. Bet. Arrays:"), " Normalization Between Arrays; Normalizes expression intensities so that the 
                          intensities or log-ratios have similar distributions across a set of arrays.")
                )
                
        )
      )
    )
  )
}


#' getQAText
#' Some questions and answers
#'
#' @return help text for QA
#'
#' @examples
#'     x<- getQAText()
#'
#' @export
getQAText<-function(){
  list(
    h3("5. Frequently asked questions (FAQ)"),
    h4("5.1 Why un-normalized counts?"),
    p("DESeq2 requires count data as input obtained from 
          RNA-Seq or another high-thorughput sequencing experiment 
          in the form of matrix values. Here we convert un-integer 
          values to integer to be able to run DESeq2. The matrix values 
          should be un-normalized, since DESeq2 model internally corrects for 
          library size. So, transformed or normalized values such as counts 
          scaled by library size should not be used as input. Please use edgeR 
          or limma for normalized counts."),
    h4("5.2 Why am I getting error while uploading files?"),
    p("* Dprofiler supports tab, comma or semi-colon separated files. However spaces or characters in numeric regions not supported and causes an error while uploading files. It is crutial to remove these kind of instances from the files before uploading files."),
    p("* Another reason of getting an error is using same gene name multiple times. This may occurs after opening files in programs such as Excel, which tends to automatically convert some gene names to dates (eg. SEP9 to SEP.09.2018). This leads numerous problems therefore you need to disable these kind of automatic conversion before opening files in these kind of programs."),
    p("* Some files contain both tab and space as an delimiter which lead to error. It is required to be cleaned from these kind of files before loading."),
    h4("5.3 Why some columns not showed up after upload?"),
    p("If a character in numeric area or space is exist in one of your column, either column will be eliminated or you will get an error. Therefore it is crutial to remove for these kind of instances from your files before uploading."),
    h4("5.4 Why am I getting error while uploading CSV/TSV files exported from Excel?"),
    p("* You might getting an error, because of using same gene name multiple times. This may occurs after opening files in programs such as Excel, which tends to automatically convert some gene names to dates (eg. SEP9 to SEP.09.2018). Therefore you need to disable these kind of automatic conversion before opening files in these kind of programs."),
    h4("5.5 Why can't I see all the background data in Main Plots?"),
    p("In order to increase the performance, by default 10% of non-significant(NS) genes are used to generate plots. We strongly suggest you to use all of the NS genes in your plots while publishing your results. You can easily change this parameter by clicking **Main Options** button and change Background Data(%) to 100% on the left sidebar."),
    h4("5.6 Why am I getting error when I click on DE Genes in Go Term Analysis?"),
    p("To start ", strong("Go Term"), " analysis, it is important to select correct organism from ", strong("Choose an organism"), " field. After selecting other desired parameters, you can click ", strong("Submit")," button to run Go Term analysis. After this stage, you will able to see", strong(" categories")," regarding to your selected gene list in the ", strong("Table")," Tab. Once you select this category, you can click DE Genes button to see gene list regarding to selected category."),
    h4("5.7 How to download selected data from Main plots/QC Plots/Heatmaps?"),
    p("First, you need to choose ", strong("Choose dataset"), " field as ",strong("selected")," under ",strong("Data Options")," in the left sidebar. When you select this option, new field: ",strong("The plot used in selection")," will appear under ", strong("Choose dataset")," field. You need to specify the plot you are interested from following options: Main plot, Main Heatmap, QC Heatmap. Finally you can click ", strong("Download Data"), " button to download data, or if you wish to see the selected data, you can click ",strong("Tables")," tab.")
  )
}

