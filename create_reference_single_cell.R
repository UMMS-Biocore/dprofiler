# library
library(cluster)
library(xbioc)
library(MuSiC)
library(VennDiagram)
library(SignallingSingleCell)
library(Seurat)
library(dplyr)

# transform Seurat into ExpressionSet
load("demo/subset_integrated.rdata")

# meta data Samples
metadata <- srt_subset@meta.data
metadata$Sample <- NA
rownames_metadata <- rownames(metadata)
rownames_metadata_indrop <- rownames_metadata[metadata$Dataset=="inDrop"]
rownames_metadata_indrop <- sapply(rownames_metadata_indrop, function(x){
  temp <- strsplit(x, split = "_")[[1]]
  return(temp[1])
})
metadata$Sample[metadata$Dataset=="inDrop"] <- rownames_metadata_indrop

rownames_metadata_fibro <- rownames_metadata[metadata$Dataset=="Fibroblast"]
rownames_metadata_fibro <- sapply(rownames_metadata_fibro, function(x){
  temp <- strsplit(x, split = "_")[[1]]
  return(temp[2])
})
rownames_metadata_fibro <- factor(rownames_metadata_fibro)
levels(rownames_metadata_fibro) <- c("Old1","Old2","Old3","Young1","Young2")
rownames_metadata_fibro <- as.character(rownames_metadata_fibro)
metadata$Sample[metadata$Dataset=="Fibroblast"] <- rownames_metadata_fibro

# metadata Cell Types
metadata$CellType <- metadata$orig_CellType
# metadata$SubCellType <- metadata$orig_SubCellType
metadata$orig_CellType <- NULL
metadata$orig_SubCellType <- NULL
metadata$SubCellType <- NULL

# tailor celltype names
metadata$CellType[grepl("Smooth", metadata$CellType)] <- "SMC"
metadata$CellType <- gsub("Vascular SMC", "VascularSMC", metadata$CellType)
metadata$CellType <- gsub("Macrophages(Hofbauer Cells)", "HofbauerCellsMacrophages", metadata$CellType)
metadata$CellType <- gsub("Lymphatic EC", "LymphaticEC", metadata$CellType)
metadata$CellType <- gsub("Mast Cells", "MastCells", metadata$CellType)
metadata$CellType <- gsub("KRT", "Keratinocytes", metadata$CellType)
metadata$CellType <- gsub("MEL", "Melanocytes", metadata$CellType)
metadata$CellType <- gsub("FB", "Fibroblasts", metadata$CellType)
metadata$CellType <- gsub("TC", "T-cells", metadata$CellType)
metadata$CellType <- gsub("DC", "Dendritic-cells", metadata$CellType)
metadata$CellType <- gsub("MAC", "Macrophages", metadata$CellType)

# meta data Embeddings
metadata$x <- subset_embedding[,1]
metadata$y <- subset_embedding[,2]

# write to seurat
srt_subset@meta.data <- metadata

# count data
srt <- ExpressionSet(assayData=as.matrix(srt_subset@assays$RNA@data))
pData(srt) <- metadata
plot_tsne_metadata(srt, color_by = "CellType", title = "ident",
                   legend_dot_size = 5, text_sizes = c(20, 10, 5, 10, 15, 15))

# calculate markers and store in fData
srt <- id_markers(srt, id_by = "CellType")

# Find markers
Idents(srt_subset) <- "CellType"
markers <- FindAllMarkers(srt_subset, logfc.threshold = 0, return.thresh = 1)
saveRDS(markers, "demo/seurat_markers_celltype.rds")

# configure markers
# _marker_score_CellType
markers_pval <- markers[markers$avg_log2FC > 0 & markers$p_val_adj < 0.01, ]
celltypes <- unique(srt_subset@meta.data$CellType)
rownames_fdata <- rownames(fData(srt))
fdata <- fData(srt)
for(i in 1:length(celltypes)){
  cell_population <- markers_pval[markers_pval$cluster==celltypes[i],]
  cell_order <- order(order(as.numeric(cell_population$avg_log2FC), decreasing = TRUE))
  fdata[[paste0(celltypes[i],"_marker_score_CellType")]] <- cell_order[match(rownames_fdata, cell_population$gene, nomatch = NA)]
}
fData(srt) <- fdata

# save reference file
saveRDS(srt, "demo/demodata_trimsc_integrated.rds")

# normalized integrated library sizes
exprs_srt <- exprs(srt)
metadata <- pData(srt)
metadata$nCount_integratedRNA_norm <- colSums(exprs_srt)
metadata <- as_tibble(metadata) %>% group_by(CellType) %>% mutate(nCount_integratedRNA_normmean = mean(nCount_integratedRNA_norm))
pData(srt) <- metadata
temp <- apply(exprs_srt, 1,function(x){
  return((x/metadata$nCount_integratedRNA_normmean)*100)
})
srt <- ExpressionSet(assayData=t(temp))
pData(srt) <- data.frame(metadata)
rownames(srt) <- colnames(temp)
colnames(srt) <- rownames(temp)
pData(srt)$UMInorm <- colSums(exprs(srt))
plot_density_ridge(srt, color_by = "CellType", val = "UMInorm") + theme(legend.position = "none")
