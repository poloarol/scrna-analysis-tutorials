
# Identification of distict cell populations and key genetic mechanims through
# scRNA-seq

# Samples are collected from several environments

# Integrate data from several patients and correct for batch effects

# 1. When to use integration: -------------
## - Multiple RNA-seq datasets from different samples
## - Cell label transfer - transfer cell type classifications from a reference to a query dataset
## - Integration of multimodal single cell data (scRNA-seq and scATAC-seq), into a multi-omic dataset
##   signals collected from separate assays
## - Integration of scRNA-seq and spatial expression data
##   Integrate topological arrangement of cells in tissues with gene expression data


# 2. Types of integration: ----------------
## - Horizontal: Same modality measured across different types of cells
## - Vertical: Multiple modalities profilied simultaneously from same cells
## - Diagonal: Different modalities from different cells (scRNA-seq and scATAC-seq)

# 3. Batch correction methods: ------------
## - Methods: Seurat v3, Harmony

# 4. Data integration: --------------------

if (!require('pacman')) install.packages("pacman")

# Load contributed packages with pacman
pacman::p_load(Seurat, tidyverse, ggplot, gridExtra, archive)

## a. Extracting files -------

# Untar GSE___ file

main_path <-"C:\\Users\\user\\OneDrive\\Desktop\\learning-R-lang\\scRNA-seq\\datasets"
untar(paste0(main_path, "\\", "GSE180665_RAW.tar"), exdir = main_path)

tar_files <- list.files(main_path, pattern="*.tar.gz")

for (tar_file in tar_files)
{
  untar(paste0(main_path, "\\", tar_file), exdir=paste0(main_path, "\\GSE180655"))
}

## b. Load data and create Seurat objects ----------

data <- list()
conditions <- list.files(paste0(main_path, "\\", "GSE180655"))

for (condition in conditions)
{
  tmp <- paste0(main_path, '\\GSE180655\\', condition)
  cts <- ReadMtx(mtx=paste0(tmp, "\\matrix.mtx.gz"),
                              cells=paste0(tmp, '\\barcodes.tsv.gz'),
                              features=paste0(tmp, "\\features.tsv.gz"))
  
  data <- append(data, CreateSeuratObject(counts=cts))

}

# Merge datasets

merge_seurat <- merge(data[[1]], y=c(data[[2]], data[[3]], data[[4]], data[[5]], data[[6]], data[[7]]), 
               add.cell.ids=c("HB17_Background", "HB17_Tumor", "HB17_PDX",
                              "HB30_Tumor", "HB30_PDX", "HB53_Background",
                              "HB53_Tumor"), project="HB") # Need to find a way to loop and perform this operation
merge_seurat

## c. QC & filter -------------

View(merge_seurat@meta.data)

## Create a smaple column
merge_seurat$sample <- rownames(merge_seurat@meta.data)

# split sample column
merge_seurat@meta.data <- separate(merge_seurat@meta.data, 
                             col = "sample", into=c("Patient", "Type", "Barcode"), sep="_")
View(merge_seurat@meta.data)

# calculate mitochodrial percentage
merge_seurat$mitoPercent <- PercentageFeatureSet(merge_seurat, pattern="^MT-")

VlnPlot(merge_seurat, feature=c("nFeature_RNA", "nCount_RNA", "mitoPercent"), ncol=3)
FeatureScatter(merge_seurat, feature1="nCount_RNA", feature2 ="nFeature_RNA") + geom_smooth(method="lm")

## d. filtering -----------

merge_seurat_filter <- subset(merge_seurat, subset=nCount_RNA > 800 & nFeature_RNA > 500 & mitoPercent < 10)
merge_seurat_filter

# Perform standard workflow steps to figure out if we see any batch effects -----
merge_seurat_filter <- NormalizeData(object = merge_seurat_filter)
merge_seurat_filter <- FindVariableFeatures(object = merge_seurat_filter)
merge_seurat_filter <- ScaleData(object = merge_seurat_filter)
merge_seurat_filter <- RunPCA(object = merge_seurat_filter)

ElbowPlot(merge_seurat_filter)
merge_seurat_filter <- FindNeighbors(object = merge_seurat_filter, dims=1:20)
merge_seurat_filter <- FindClusters(object = merge_seurat_filter)
merge_seurat_filter <- RunUMAP(object = merge_seurat_filter, dims=1:20)

# Plot
p1 <- DimPlot(merge_seurat_filter, reduction = "umap", group.by = "Patient")
p2 <- DimPlot(merge_seurat_filter, reduction = "umap", group.by = "Type", cols = c("red", "green", "blue"))

grid.arrange(p1, p2, ncol=2)

# Perform integration to correct for bacth effects ----------
obj.list <- SplitObject(merge_seurat_filter, split.by="Patient")

for (i in 1:length(obj.list))
{
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures((object = obj.list[[i]]))
}

# select integration features
features <- SelectIntegrationFeatures(object.list = obj.list)

# Find integration anchors (CCA): Canonical Correlation Analysis
anchors <-  FindIntegrationAnchors(object.list = obj.list, anchor.features = features)


# integrate data
seurat.integrated <- IntegrateData(anchorset = anchors)

# Scale data, run PCA and UMAP and visualize integrated data
seurat.integrated <- ScaleData(object = seurat.integrated)
seurat.integrated <- RunPCA(object= seurat.integrated)
seurat.integrated <- RunUMAP(object =seurat.integrated, dims=1:50)

p3 <- DimPlot(seurat.integrated, reduction="umap", group.by = "Patient")
p4 <- DimPlot(seurat.integrated, reduction="umap", group.by = "Type", cols=c("red", "green", "blue"))


grid.arrange(p3, p4, ncol=2, ncol=2)

p <- grid.arrange(p1, p2, p3, p4, ncol=2, nrow=2)





















