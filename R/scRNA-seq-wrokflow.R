

if (!require('pacman')) install.packages("pacman")
 
# Load contributed packages with pacman
pacman::p_load(Seurat, tidyverse)


# Load datasets
nsclc.sparse.matrix <- Read10X_h5(filename="C:\\Users\\user\\OneDrive\\Desktop\\learning-R-lang\\scRNA-seq\\datasets\\40k_NSCLC_DTC_3p_HT_nextgem_Multiplex_count_raw_feature_bc_matrix.h5")
cts <- nsclc.sparse.matrix$`Gene Expression`

# Initialize the Seurat object with the raw (non-normalized data).
nsclc.seurat.obj <- CreateSeuratObject(cts, project = "NSCLC", min.cells = 3, min.features = 200)
nsclc.seurat.obj

# 1. QC ----------------------
# % MT reads
nsclc.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(nsclc.seurat.obj, pattern="^MT-")
View(nsclc.seurat.obj@meta.data)

VlnPlot(nsclc.seurat.obj, feature=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
FeatureScatter(nsclc.seurat.obj, feature1="nCount_RNA", feature2 ="nFeature_RNA") + geom_smooth(method="lm")

# 2. Filtering ---------------
nsclc.seurat.obj <- subset(nsclc.seurat.obj, subset=nFeature_RNA>200 & nFeature_RNA<2500 & percent.mt<5, idents=c("NSCLC"))
## Can remove doublets here, using DoubletFinder

# 3. Normalize data ----------
## NormalizeData(nsclc.seurat.obj, normalization.method="logNormalize", scale.factor=10000)
## scTransform has been shown to work best

nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj)

# 4. Identify highly variable features --------------
nsclc.seurat.obj <- FindVariableFeatures(nsclc.seurat.obj, selection.method = "vst", nfeatures = 2000)

## Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(nsclc.seurat.obj), 10)
top10

## plot variable features with and without labels
plot1 <- VariableFeaturePlot(nsclc.seurat.obj)
LabelPoints(plot=plot1, points = top10, repel = TRUE)

# 5. Scaling -------------
all.genes <- rownames(nsclc.seurat.obj)
nsclc.seurat.obj <- ScaleData(nsclc.seurat.obj, repel=TRUE)

# 6. Perform linear dimensionality reduction ---------------
nsclc.seurat.obj <- RunPCA(nsclc.seurat.obj, feature=VariableFeatures(object=nsclc.seurat.obj))
DimHeatmap(nsclc.seurat.obj, dims=1, cells=500, balanced=TRUE)

## Determine dimensionality of data i.e Find statistically significant PCs
ElbowPlot(nsclc.seurat.obj)

# 7. Clustering -----------
nsclc.seurat.obj <- FindNeighbors(nsclc.seurat.obj, dims=1:15)

## Visualize PCA results
print(nsclc.seurat.obj[["pca"]], dims=1, nfeatures=5)

## Understanding resolution
nsclc.seurat.obj <- FindClusters(nsclc.seurat.obj, resolution = c(0, 0.1, 0.3, 0.5, 0.7, 1))
View(nsclc.seurat.obj@meta.data)

DimPlot(nsclc.seurat.obj, group.by="RNA_snn_res.0.3", label=TRUE)

## Setting identity of clusters
Idents(nsclc.seurat.obj) <- "RNA_snn_res.0.1"
Idents(nsclc.seurat.obj)


# Non-linear dimensionality reduction
install.packages("umap")
library(umap)

nsclc.seurat.obj <- RunUMAP(nsclc.seurat.obj, dims=1:15)
DimPlot(nsclc.seurat.obj, reduction = "umap")











