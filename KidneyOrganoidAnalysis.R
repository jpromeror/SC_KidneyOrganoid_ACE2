############
## Script for Single-cell Kidney Organoid Analysis

#Inhibition of SARS-CoV-2 infections in engineered human tissues using clinical-grade soluble human ACE2

# Juan P. Romero
# jromeror@unav.es
############

## Load R packages

library(Seurat)

## Read .h5 file (Output from Cell Ranger)
KidneyOrganoid<-Read10X_h5("/Users/jpromero/Desktop/KidneyOrg/KidneyOrganoid_FilteredGeneBCMatrices.h5")

## Create Seurat Object with initial filtering
KidneyOrganoid <- CreateSeuratObject(counts = KidneyOrganoid, project = "KidneyOrganoid_ACE2", min.cells = 3, min.features = 400)

# Add %MT as QC value
KidneyOrganoid[["percent.mt"]] <- PercentageFeatureSet(KidneyOrganoid, pattern = "^MT-")

## QC Metrics Plots
VlnPlot(KidneyOrganoid, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.3)

## Get QC Thresholds
quantile(KidneyOrganoid@meta.data$nCount_RNA,c(0.025,0.975))
quantile(KidneyOrganoid@meta.data$nFeature_RNA,c(0.025,0.975))

## QC Plots
plot(KidneyOrganoid@meta.data$nCount_RNA,KidneyOrganoid@meta.data$nFeature_RNA,pch=16,cex=0.7,bty="n")
abline(h=c(488,5653),v=c(667,23108),lty=2,lwd=1,col="red")

## Filtering based on QC parameters
KidneyOrganoid <- subset(KidneyOrganoid, subset = nFeature_RNA > 488 & nFeature_RNA < 5653 & nCount_RNA > 667 & nCount_RNA < 23108 & percent.mt < 50)

## Log Normalization 
KidneyOrganoid<-NormalizeData(KidneyOrganoid)

## Scale Data
KidneyOrganoid <- ScaleData(KidneyOrganoid, features = rownames(KidneyOrganoid))

## Cell Cycle Effect
KidneyOrganoid<-CellCycleScoring(KidneyOrganoid,s.features = cc.genes$s.genes,g2m.features = cc.genes$g2m.genes)
KidneyOrganoid <- RunPCA(KidneyOrganoid, features = unlist(cc.genes))
DimPlot(KidneyOrganoid, reduction = "pca",dims = c(1,2),group.by = "Phase")

## SCTransform 
KidneyOrganoid<-SCTransform(KidneyOrganoid,vars.to.regress = c("S.Score","G2M.Score","percent.mt","nFeature_RNA"))

## PCA
KidneyOrganoid <- RunPCA(KidneyOrganoid, features = VariableFeatures(object = KidneyOrganoid))
DimPlot(KidneyOrganoid, reduction = "pca",dims = c(1,2))

## Some Plots
VizDimLoadings(KidneyOrganoid, dims = 1:2, reduction = "pca")
DimPlot(KidneyOrganoid, reduction = "pca",dims = c(1,2))

## Selecting PCA Components
ElbowPlot(KidneyOrganoid,ndims = 30)

## Clustering
KidneyOrganoid <- FindNeighbors(KidneyOrganoid, dims = 1:20)
KidneyOrganoid <- FindClusters(KidneyOrganoid, resolution = 0.4)

# Non Linear Dimensional Reduction
KidneyOrganoid <- RunUMAP(KidneyOrganoid, dims = 1:20)

# UMAP plot
colss<-c("#A6CEE3", "#1F78B4", "#08306B", "#B2DF8A", "#006D2C", "#8E0152",
         "#DE77AE", "#CAB2D6", "#6A3D9A", "#FB9A99", "#E31A1C", "#B15928",
         "#619CFF","#FF67A4","#00BCD8")

DimPlot(KidneyOrganoid, reduction = "umap",label = T,col=colss)

# Feature Plots on interesting genes
FeaturePlot(KidneyOrganoid,c("ACE2"),cols = c("lightgray","red"),order = T)
FeaturePlot(KidneyOrganoid,c("SLC3A1","SLC27A2","PODXL","NPHS2","NPHS1","CLDN4","MAL","CD93"),cols = c("lightgray","red"),order = T)

# Find Cluster Markers
KidneyOrganoid.markers <- FindAllMarkers(KidneyOrganoid, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
KidneyOrganoid.markers$Specific<-KidneyOrganoid$pct.1-KidneyOrganoid$pct.2


