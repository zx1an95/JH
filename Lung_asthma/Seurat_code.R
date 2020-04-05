rm(list = ls())   #remove all things in the environment
load("D:/Prof Campbell Project/Source data/GSE130148_raw_counts.RData/GSE130148_raw_counts.RData") #column=cell,row=gene
anno_data <- read.csv("D:/Prof Campbell Project/Source data/GSE130148_barcodes_cell_types.txt/GSE130148_barcodes_cell_types.txt",sep = '\t')
library(ggplot2)
library(sva)
samples <- anno_data$GEO_Sample
#new.dat <- as.matrix(raw_counts)
#new.dat <- ComBat(new.dat,batch = batch)
#new.dat <- Matrix(new.dat,sparse = TRUE)
#Seurat
library(dplyr)
library(Seurat)
lung_resection <- CreateSeuratObject(counts = raw_counts, project = "lung_resection", 
                                     min.cells = 3, min.features = 200)

lung_resection[["percent.mt"]] <- PercentageFeatureSet(lung_resection, pattern = "^MT-")
head(lung_resection@meta.data, 5)
VlnPlot(lung_resection, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(lung_resection, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(lung_resection, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
lung_resection <- subset(lung_resection, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
lung_resection <- NormalizeData(lung_resection, normalization.method = "LogNormalize", scale.factor = 10000)
lung_resection <- FindVariableFeatures(lung_resection, selection.method = "vst", nfeatures = 4500)
top10 <- head(VariableFeatures(lung_resection), 10)
plot3 <- VariableFeaturePlot(lung_resection)
plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE)
CombinePlots(plots = list(plot3, plot4))
all.genes <- rownames(lung_resection)
lung_resection <- ScaleData(lung_resection, features = all.genes)
#lung_resection@assays$RNA@scale.data <- ComBat(lung_resection@assays$RNA@scale.data,batch = batch)
lung_resection <- RunPCA(lung_resection, features = VariableFeatures(object = lung_resection))

check.batch.effect <- lung_resection@reductions$pca@cell.embeddings[,1:2]
ggplot(as.data.frame(check.batch.effect),aes(x = PC_1, y = PC_2))+geom_point(aes(col=batch),size=1.5)+
  labs(title="PCA plot")+xlab("PC1")+ylab("PC2")

print(lung_resection[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(lung_resection, dims = 1:2, reduction = "pca")
DimPlot(lung_resection, reduction = "pca")
DimHeatmap(lung_resection, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(lung_resection, dims = 50, cells = 500, balanced = TRUE)
lung_resection <- JackStraw(lung_resection, num.replicate = 100,dims = 30)
lung_resection <- ScoreJackStraw(lung_resection, dims = 1:25)
JackStrawPlot(lung_resection, dims = 1:30)
ElbowPlot(lung_resection,ndims = 30)
lung_resection <- FindNeighbors(lung_resection, dims = 1:43,k.param = 30)
lung_resection <- FindClusters(lung_resection,resolution = 0.8)
head(Idents(lung_resection), 5)

lung_resection <- RunTSNE(lung_resection, dims = 1:43)
lung_resection@meta.data$orig.ident <- samples
DimPlot(lung_resection, reduction = "tsne",group.by = "orig.ident")
DimPlot(lung_resection, reduction = "tsne",split.by = "orig.ident",label = TRUE)
lung_resection <- RunUMAP(lung_resection, dims = 1:43)
DimPlot(lung_resection, reduction = "umap",group.by = "orig.ident")
DimPlot(lung_resection, reduction = "umap",split.by = "orig.ident",label = TRUE)
#tsne.check.batch <- lung_resection@reductions$tsne@cell.embeddings[,1:2]
#ggplot(as.data.frame(tsne.check.batch),aes(x = tSNE_1, y = tSNE_2))+geom_point(aes(col=samples),size=1)+
  #labs(title="tSNE plot")+xlab("tSNE1")+ylab("tSNE2")


lung_section.markers <- FindAllMarkers(lung_resection, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
lung_section.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
cluster1.markers <- FindMarkers(lung_resection, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
cluster16.markers <- FindMarkers(lung_resection, ident.1 = 16, min.pct = 0.25)
head(cluster16.markers, n = 20)
FeaturePlot(lung_resection, features = c("APOC1", "MARCO", "C1QB", "C19orf59", "FABP4")) 
VlnPlot(lung_resection, features = c("APOC1", "MARCO", "C1QB", "C19orf59", "FABP4"), slot = "counts", log = TRUE)
FeaturePlot(lung_resection, features = c("SFTPB", "SFTPC", "SFTPA2", "SLPI", "NAPSA"))
VlnPlot(lung_resection, features = c("SFTPB", "SFTPC", "SFTPA2", "SLPI", "NAPSA"),slot = "counts", log = TRUE)
FeaturePlot(lung_resection, features = c("IGHA1", "IGHM", "IGHG3", "IGHG1", "IGKC"))
VlnPlot(lung_resection, features = c("IGHA1", "IGHM", "IGHG3", "IGHG1", "IGKC"),slot = "counts", log = TRUE)
FeaturePlot(lung_resection, features = c("CCL5", "TRBC2", "CD2", "IL32", "TRAC"))
VlnPlot(lung_resection, features = c("CCL5", "TRBC2", "CD2", "IL32", "TRAC"),slot = "counts", log = TRUE)
FeaturePlot(lung_resection, features = c("S100A9", "VCAN", "APOBEC3A", "CD300E", "THBS1"))
VlnPlot(lung_resection, features = c("S100A9", "VCAN", "APOBEC3A", "CD300E", "THBS1"),slot = "counts", log = TRUE)
FeaturePlot(lung_resection, features = c("SCGB1A1", "SCGB3A1", "SCGB3A2", "BPIFB1", "WFDC2"))
VlnPlot(lung_resection, features = c("SCGB1A1", "SCGB3A1", "SCGB3A2", "BPIFB1", "WFDC2"),slot = "counts", log = TRUE)
FeaturePlot(lung_resection, features = c("TPSAB1", "CPA3", "HPGDS", "KIT", "TPSD1"))
VlnPlot(lung_resection, features = c("TPSAB1", "CPA3", "HPGDS", "KIT", "TPSD1"),slot = "counts", log = TRUE)
FeaturePlot(lung_resection, features = c("PRG4", "MT2A", "MT1E", "KRT19", "SLC6A14"))
VlnPlot(lung_resection, features = c("PRG4", "MT2A", "MT1E", "KRT19", "SLC6A14"),slot = "counts", log = TRUE)
FeaturePlot(lung_resection, features = c("C9orf24", "TPPP3", "C20orf85", "RSPH1", "C11orf88"))
VlnPlot(lung_resection, features = c("C9orf24", "TPPP3", "C20orf85", "RSPH1", "C11orf88"),slot = "counts", log = TRUE)
FeaturePlot(lung_resection, features = c("EMP2", "AGER", "CAV1", "RTKN2", "CEACAM6"))
VlnPlot(lung_resection, features = c("EMP2", "AGER", "CAV1", "RTKN2", "CEACAM6"),slot = "counts", log = TRUE)
FeaturePlot(lung_resection, features = c("SPARCL1", "EPAS1", "EMP1", "VWF", "SERPINE1"))
VlnPlot(lung_resection, features = c("SPARCL1", "EPAS1", "EMP1", "VWF", "SERPINE1"),slot = "counts", log = TRUE)
FeaturePlot(lung_resection, features = c("GNLY", "NKG7", "GZMA", "FGFBP2"))
VlnPlot(lung_resection, features = c("GNLY", "NKG7", "GZMA", "FGFBP2"),slot = "counts", log = TRUE)
FeaturePlot(lung_resection, features = c("HLA-DPA1", "HLA-DQA1", "HLA-DQB1", "GPR183", "FCER1A"))
VlnPlot(lung_resection, features = c("HLA-DPA1", "HLA-DQA1", "HLA-DQB1", "GPR183", "FCER1A"),slot = "counts", log = TRUE)
FeaturePlot(lung_resection, features = c("DCN", "COL1A2", "LUM", "COL3A1", "A2M"))
VlnPlot(lung_resection, features = c("DCN", "COL1A2", "LUM", "COL3A1", "A2M"),slot = "counts", log = TRUE)
FeaturePlot(lung_resection, features = c("CCL21", "TFF3", "GNG11", "MMRN1", "LYVE1"))
VlnPlot(lung_resection, features = c("CCL21", "TFF3", "GNG11", "MMRN1", "LYVE1"),slot = "counts", log = TRUE)
FeaturePlot(lung_resection, features = c("KRT15", "BCAM", "KRT17", "DST", "KRT5"))
VlnPlot(lung_resection, features = c("KRT15", "BCAM", "KRT17", "DST", "KRT5"),slot = "counts", log = TRUE)
VlnPlot(lung_resection, features = c("S100A2", "CYR61", "NPPC", "SERPINB4", "SERPINB13"),slot = "counts", log = TRUE)
VlnPlot(lung_resection, features = c("FCGR3B","CMTM2","S100A9","CSF3R","S100A8","IFITM2"),slot = "counts", log = TRUE)
VlnPlot(lung_resection, features = c("MNDA","VNN2","CXCR2","S100P","G0S2","SELL") ,slot = "counts", log = TRUE)

top5genes <- lung_section.markers %>% group_by(cluster) %>% top_n(n = 5,wt = avg_logFC)
top20genes <- lung_section.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
Indexlist2 <- c("APOC1", "MARCO", "C1QB", "C19orf59", "FABP4",
                "SFTPB", "SFTPC", "SFTPA2", "SLPI", "NAPSA",
                "S100A9", "VCAN", "APOBEC3A", "CD300E", "THBS1",
                "IGHA1", "IGHM", "IGHG3", "IGHG1", "IGKC",
                "CCL5", "TRBC2", "CD2", "IL32", "TRAC",
                "EMP2", "AGER", "CAV1", "RTKN2", "CEACAM6",
                "C9orf24", "TPPP3", "C20orf85", "RSPH1", "C11orf88",
                "SCGB1A1", "SCGB3A1", "SCGB3A2", "BPIFB1", "WFDC2",
                "SPARCL1", "EPAS1", "EMP1", "VWF", "SERPINE1",
                "GNLY", "NKG7", "GZMA", "FGFBP2",
                "TPSAB1", "CPA3", "HPGDS", "KIT", "TPSD1",
                "PRG4", "MT2A", "MT1E", "KRT19", "SLC6A14",
                "HLA-DPA1", "HLA-DQA1", "HLA-DQB1", "GPR183", "FCER1A",
                "DCN", "COL1A2", "LUM", "COL3A1", "A2M",
                "CCL21", "TFF3", "GNG11", "MMRN1", "LYVE1")
DoHeatmap(lung_resection, features = Indexlist2) + NoLegend() 
DoHeatmap(lung_resection, features = top5genes$gene) + NoLegend()#
new.cluster.ids <- c("Macrophages", "Macrophages", "Type 2", "Neutrophils", "B cell", "T cell", 
                     "T cell", "Type 1", "Ciliated", "Type 2", "Macrophages", "Club", "Endothelium",
                     "NK cell", "Mast cell", "Mast cell", "Neutrophils", "Transformed epithelium",
                     "B cell", "Dendritic cell", "Type 2", "T cell", "Fibroblast", "Lymphatic", "Fibroblast")
top5genes$gene[top5genes$cluster=='24']
top5genes$gene[top5genes$cluster=='Basal cell']
names(new.cluster.ids) <- levels(lung_resection)
lung_resection <- RenameIdents(lung_resection, new.cluster.ids)
DimPlot(lung_resection, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()