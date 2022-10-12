#' getReferenceSingleCellRNA
#' 
#' Prepares Seurat objects for Compositional Profiling 
#'
#' @param object Seurat Object
#' @param name name of the Single Cell Reference
#' @param marker_ident Ident of Seurat Object for marker analysis
#' @param subsetProp Proportion of single cell data to downsize, if NULL the entire single cell data is used
#' @param subset_seed The seed for random sampling
#' @param subset_ident Ident to subset the single cell data, if NULL cells are randomly selected
#' @param subset_ident_min Minimum number of cells in an ident, ignored if subset_ident is NULL
#' @param subset_assay Seurat assay used in the subset
#' @param subset_data Seurat data used in the subset
#' @param subset_embedding Seurat embedding used in the subset
#' 
#' @import Seurat
#' 
#' @export
#'
#' @examples
#' # library
#' library(Seurat)
#' 
#' # prepare scRNA reference data
#' referencePBMC <- getReferenceSingleCellRNA(pbmc_small)
#'                                            
#' # prepare scRNA reference data, select idents for markers analysis and subsetting                                          
#' referencePBMC <- getReferenceSingleCellRNA(pbmc_small, marker_ident = c("RNA_snn_res.0.8","RNA_snn_res.1"), 
#'                                            subsetProp = 0.8, subset_ident = "RNA_snn_res.0.8", subset_ident_min = 10)
#'
getReferenceSingleCellRNA <- function(object, name, marker_ident = NULL, subsetProp = NULL, 
                                      subset_seed = 1, subset_ident = NULL, subset_ident_min = 10, 
                                      subset_assay = "RNA", subset_data = "counts",
                                      subset_embedding = "tsne"){
  
  if(is.null(object) || class(object) != "Seurat")
    stop("Please provide a single cell data (Seurat object)")
  
  # create single cell reference folder
  dir.create(name)
  
  # marker analysis
  if(!is.null(marker_ident)){
    
    if(!any(colnames(object@meta.data) %in% marker_ident))
      stop("Please select an existing identification from the Seurat object")
    
    cat("Finding markers for idents:", paste(marker_ident, collapse = " "), "\n")
    markers_list <- NULL
    for(cur_ident in marker_ident){
      Idents(object) <- cur_ident
      markers <- FindAllMarkers(object)
      markers$Level <- cur_ident
      markers_list <- rbind(markers_list, markers)
    } 
    saveRDS(markers_list, file = paste0(name, "/markerdata.rds"))
  } else {
    warning("No identification is provided for single cell data, deconvolution cannot be conducted with markers")
    markers_list <- NULL
  }
  
  # subset of the data
  if(!is.null(subsetProp)){
    
    if(!is.numeric(subsetProp))
      stop("subsetProp should be between 0 and 1")
    
    if(subsetProp < 0 || subsetProp > 1)
      stop("subsetProp should be between 0 and 1")
    
    # random subset of cells
    cells_object <- Cells(object)
    set.seed(subset_seed)
    if(is.null(subset_ident)){
      subset_size <- round(length(cells_object)*subsetProp)
      cat("Random Subset of Seurat Object: ", subset_size, " cells are selected! \n")
      cells_subset <- sample(cells_object, subset_size, replace = FALSE) 
    } else {
      cat("Ident based subsampling: \n")
      identtypes <- unique(object@meta.data[[subset_ident]])
      cells_subset <- sapply(identtypes, function(x) {
        cells_object_identtypes <- cells_object[object@meta.data[[subset_ident]] == x]
        subset_size <- max(round(length(cells_object_identtypes)*subsetProp), 
                           min(subset_ident_min, length(cells_object_identtypes)))
        cat("Random Subset of group ", x, ": ", subset_size, " cells are selected! \n")
        cells_object_identtypes_subset <- sample(cells_object_identtypes, subset_size, replace = FALSE)
        return(cells_object_identtypes_subset)
      })
      cells_subset <- unlist(cells_subset)
    }
    object_subset <- object[,cells_subset]
  } else {
    object_subset <- object
  }
  
  # turn into ExpressionSet
  countdata <- as.matrix(GetAssayData(object_subset[[subset_assay]], slot = subset_data))
  metadata <- object_subset@meta.data
  metadata$x <- Embeddings(object_subset, reduction = subset_embedding)[,1]
  metadata$y <- Embeddings(object_subset, reduction = subset_embedding)[,2]
  metadata <- AnnotatedDataFrame(metadata)
  scdata <- ExpressionSet(countdata, metadata)
  saveRDS(scdata, file = paste0(name, "/countdata.rds"))
  
  # returning object
  return(list(scdata = scdata, markers = markers_list))
}

#' getReferenceProfileRNA
#' 
#' Prepares Seurat objects for Compositional Profiling 
#'
#' @param data profiling data
#' @param meta.data profiling metadata
#' @param name name of the Bulk Reference
#' @param sample_id the field in the metadata that identifies a sample/biosample, if null, the first field will be the sample_id
#' @param batch the batch term from the meta.data
#' @param secondary_batch specify is there exists a secondary batch within the primary batch, e.g. Patient, or Donor
#' @param subset_profiledata a list of subsets based on meta.data fields, the names of each element should correspond to a field in the meta.data
#' @param DEfield the meta.data field to indicate a comparison, this should be a list with one element that this element should be a vector of two characters indicating the conditions
#' @param DEparams DE Analysis and Computational Profiling parameters 
#' @param max_count remove genes whose maximum count is below a small count (e.g. 10)
#' 
#' @export
#'
getReferenceProfileRNA<- function(data, meta.data, name, sample_id = NULL, batch = NULL, secondary_batch = NULL, subset_profiledata = NULL, 
                                  DEfield = NULL,
                                  DEparams = c("DESeq2", "parametric", FALSE, "Wald", "None", "Silhouette", "auto", "Log2FC+Padj", 1, 0.05, ""),             
                                  max_count = 10){
  
  # create single cell reference folder
  dir.create(name)
  
  # subset data based on metadata
  set.seed(1)
  if(!is.null(subset_profiledata)){
    if(all(names(subset_profiledata) %in% colnames(meta.data))){
      for(i in 1:length(subset_profiledata)){
        meta.data <- meta.data[meta.data[,names(subset_profiledata)[i]] %in% subset_profiledata[[i]],]
      }
    } else {
      no_names <-  names(subset_profiledata)[!names(subset_profiledata) %in% colnames(meta.data)]
      stop(paste("Fields: ", paste(no_names, collapse = ", "), " are missing in the meta.data"))
    }
  } else {
    warning("No subsetting has been provided. The entire data will be used for profiling")
  }
  
  if(!is.null(DEfield) & is.list(DEfield)){
    if(length(DEfield) == 1 | length(DEfield[[1]]) == 2){
      meta.data <- meta.data[meta.data[,names(DEfield)[1]] %in% DEfield[[1]],]
    } else {
      stop("DE field should be a list of one element and that element should be a vector of two characters: two conditions")
    }
  } else {
    stop("DE field is empty or not a list")
  }
  
  # do all batches have comparison conditions
  if(batch %in% colnames(meta.data)){
    
    # select batch that has all secondary batch cases
    crosstable <- table(meta.data[,batch], meta.data[,names(DEfield)[1]])
    crosstable <- data.frame(rbind(crosstable))
    ind <- apply(crosstable,1,function(x) all(x > 0))
    selected_batches <- rownames(crosstable)[ind]
    meta.data <- meta.data[meta.data[,batch] %in% selected_batches,]
    
  }
  
  # sample ID
  if(is.null(sample_id)){
    sample_id <- colnames(meta.data)[1]
  }
  
  # overlap samples in profile data and metadata
  common_samples <- intersect(colnames(data), meta.data[,sample_id])
  datax_int <- data[,common_samples]
  meta.data <- meta.data[meta.data[,sample_id] %in% common_samples,]
  datax_int <- datax_int[,meta.data[,sample_id]]
  
  # remove max count
  max_count <- apply(datax_int,1,max)
  datax_int <- datax_int[which(max_count > 10),]
  
  # Batch Correct 
  datax_correct <- NULL
  if(batch %in% colnames(meta.data)){
    
    # is a secondary batch term specified 
    if(!is.null(secondary_batch)){
      
      # correct for secondary batch term for all batches
      unique_batch <- unique(meta.data[,batch])
      for(i in 1:length(unique_batch)){
        datax_temp <- datax_int[,meta.data[,batch] == unique_batch[i]]
        second_batch <- meta.data[meta.data[,sample_id] %in% colnames(datax_temp), secondary_batch]
        genes <- rownames(datax_temp)
        samples <- colnames(datax_temp)
        datax_temp <- apply(datax_temp, 2, function(x) as.integer(x) + 
                              runif(1, 0, 0.01))
        datax_temp = sva::ComBat(dat=as.matrix(datax_temp), batch=second_batch)
        datax_temp <- t(apply(datax_temp,1,as.integer))
        datax_temp <- apply(datax_temp, 2, function(x) ifelse(x < 0, 0, x))
        rownames(datax_temp) <- genes
        colnames(datax_temp) <- samples
        datax_correct <- cbind(datax_correct, datax_temp)
      }
    } else {
      datax_correct <- datax_int
    }
    
    # correct for the primary batch
    genes <- rownames(datax_correct)
    samples <- colnames(datax_correct)
    datax_correct <- apply(datax_correct, 2, function(x) as.integer(x) + 
                             runif(1, 0, 0.01))
    datax_correct = sva::ComBat(dat=as.matrix(datax_correct), batch=meta.data[,batch])
    datax_correct <- t(apply(datax_correct,1,as.integer))
    datax_correct <- apply(datax_correct, 2, function(x) ifelse(x < 0, 0, x))
    rownames(datax_correct) <- genes
    colnames(datax_correct) <- samples
    
  } else{
    datax_correct <- datax_int
    warning("Batch term in not specified, batch correction will not be applied")
  }
  
  # Computational Phenotypic Profiling 
  # ManualDEgenes <- NULL
  # TopStat <- NA
  columns <- colnames(datax_correct)
  conds <- meta.data[,names(DEfield)[1]]
  
  # configure conditions
  conds <- factor(conds)
  levels(conds) <- c("cond1","cond2")
  conds <- as.character(conds)
  columns <- columns[order(conds)]
  conds <- conds[order(conds)]
  
  # Profiling
  IterResults <- runIterDE(datax_correct, meta.data, columns, conds, DEparams, verbose = TRUE, visualize_prefix = names(DEfield)[1])
  # FinalScores <- getFinalScores(IterResults, datax_correct, columns, conds, DE, ManualDEgenes = ManualDEgenes, TopStat = TopStat)
  
  # make profiles
  datax_correct_cleaned <- datax_correct[,!colnames(datax_correct) %in% IterResults$cleaned_columns]
  saveRDS(datax_correct_cleaned, file = paste0(name, "countdata.rds"))
  
  # PCA embedding and metadata
  datax_correct_cleaned <- datax_correct_cleaned[IterResults$IterDEgenes,]
  datax_pr <- apply(datax_correct_cleaned, 1, scale)
  prtemp <- prcomp(datax_pr)
  datax_pr <- data.frame(prtemp$x)
  var_perc <- round(prtemp$sdev^2/sum(prtemp$sdev^2)*100)
  meta.data_cleaned <- meta.data[!meta.data[,sample_id] %in% IterResults$cleaned_columns,]
  meta.data_cleaned <- data.frame(meta.data_cleaned, x = datax_pr$PC1, y = datax_pr$PC2)
  colnames(meta.data_cleaned)[colnames(meta.data_cleaned) %in% c("x","y")] <- paste0("PC", 1:2, "_", var_perc[1:2])
  saveRDS(meta.data_cleaned, file = paste0(name,"/metadata.rds"))
  
  # marker table
  prof_marker_table <- IterResults$IterDEResults
  prof_marker_table$Level <- paste(names(DEfield)[1], paste(DEfield[[1]], collapse = "vs"), sep = "_")
  saveRDS(prof_marker_table, file = paste0(name, "/markerdata.rds"))
  
  # returning object
  return(list(profiledata = datax_correct_cleaned, meta.data = meta.data_cleaned, markers = prof_marker_table))
}