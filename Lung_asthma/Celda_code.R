load("GSE130148_raw_counts.RData")
lung_resection <- CreateSeuratObject(counts = raw_counts, project = "lung_resection", 
                                     min.cells = 3, min.features = 250)
lung_resection[["percent.mt"]] <- PercentageFeatureSet(lung_resection, pattern = "^MT-")
lung_resection <- subset(lung_resection, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
lung_resection <- NormalizeData(lung_resection, normalization.method = "LogNormalize", scale.factor = 10000)
lung_resection <- FindVariableFeatures(lung_resection, selection.method = "vst", nfeatures = 4500)

topgenes <- lung_resection@assays$RNA@var.features
library(celda)
new_lung_resection <- as.matrix(lung_resection@assays$RNA@counts)
new_lung_resection <- new_lung_resection[topgenes,]
moduleSplit <- recursiveSplitModule(counts = new_lung_resection,
                                    initialL = 3, maxL = 150)   #L = 85
plotGridSearchPerplexity(celdaList = moduleSplit)
moduleSplitSelect <- subsetCeldaList(moduleSplit, params = list(L = 85))
cellSplit <- recursiveSplitCell(counts = new_lung_resection,
                                initialK = 3,
                                maxK = 45,
                                yInit = clusters(moduleSplitSelect)$y) #K = 30
plotGridSearchPerplexity(celdaList = cellSplit)
goodceldaModel <- subsetCeldaList(celdaList = cellSplit, params = list(K = 30, L = 85))
factorized <- factorizeMatrix(new_lung_resection, celdaMod = goodceldaModel)
names(factorized)
cellPop <- factorized$proportions$cellPopulation
topGenes <- topRank(matrix = factorized$proportions$module,
                    n = 5, threshold = NULL)
celdaHeatmap(counts = new_lung_resection, celdaMod = goodceldaModel, nfeatures = 5)
tsne <- celdaTsne(counts = new_lung_resection, celdaMod = goodceldaModel)
plotDimReduceCluster(dim1 = tsne[, 1],
                     dim2 = tsne[, 2],
                     cluster = clusters(goodceldaModel)$z)
plotDimReduceFeature(dim1 = tsne[, 1],
                     dim2 = tsne[, 2],
                     counts = new_lung_resection,
                     features = "DCN")
celdaProbabilityMap(counts = new_lung_resection, celdaMod = goodceldaModel)
moduleHeatmap(counts = new_lung_resection, celdaMod = goodceldaModel,
              featureModule = 55, topCells = 100)
genes <- c(topGenes$names$L76, topGenes$names$L78)
geneIx <- which(rownames(new_lung_resection) %in% genes)
normCounts <- normalizeCounts(counts = new_lung_resection, scaleFun = scale)
plotHeatmap(counts = normCounts,
            z = clusters(goodceldaModel)$z,
            y = clusters(goodceldaModel)$y,
            featureIx = geneIx,
            showNamesFeature = TRUE)

plotDimReduceCluster(dim1 = tsne[, 1], dim2 = tsne[, 2], cluster = clusters(goodceldaModel)$z, specificClusters = c(2,3))
plotDimReduceCluster(dim1 = tsne[, 1], dim2 = tsne[, 2], cluster = clusters(goodceldaModel)$z)
plotDimReduceFeature(dim1 = tsne[, 1],dim2 = tsne[, 2],counts = new_lung_resection,features = c("APOC1", "MARCO", "C1QB", "C19orf59", "FABP4"))
plotDimReduceFeature(dim1 = tsne[, 1],dim2 = tsne[, 2],counts = new_lung_resection,features = c("SFTPB", "SFTPC", "SFTPA2", "SLPI", "NAPSA"))
plotDimReduceFeature(dim1 = tsne[, 1],dim2 = tsne[, 2],counts = new_lung_resection,features = c("IGHA1", "IGHM", "IGHG3", "IGHG1", "IGKC"))
plotDimReduceFeature(dim1 = tsne[, 1],dim2 = tsne[, 2],counts = new_lung_resection,features = c("CCL5", "TRBC2", "CD2", "IL32", "TRAC"))
plotDimReduceFeature(dim1 = tsne[, 1],dim2 = tsne[, 2],counts = new_lung_resection,features = c("S100A9", "VCAN", "APOBEC3A", "CD300E", "THBS1"))
plotDimReduceFeature(dim1 = tsne[, 1],dim2 = tsne[, 2],counts = new_lung_resection,features = c("SCGB1A1", "SCGB3A1", "SCGB3A2", "BPIFB1", "WFDC2"))
plotDimReduceFeature(dim1 = tsne[, 1],dim2 = tsne[, 2],counts = new_lung_resection,features = c("TPSAB1", "CPA3", "HPGDS", "KIT", "TPSD1"))
plotDimReduceFeature(dim1 = tsne[, 1],dim2 = tsne[, 2],counts = new_lung_resection,features = c("PRG4", "MT2A", "MT1E", "KRT19", "SLC6A14"))
plotDimReduceFeature(dim1 = tsne[, 1],dim2 = tsne[, 2],counts = new_lung_resection,features = c("C9orf24", "TPPP3", "C20orf85", "RSPH1", "C11orf88"))
plotDimReduceFeature(dim1 = tsne[, 1],dim2 = tsne[, 2],counts = new_lung_resection,features = c("EMP2", "AGER", "CAV1", "RTKN2", "CEACAM6"))
plotDimReduceFeature(dim1 = tsne[, 1],dim2 = tsne[, 2],counts = new_lung_resection,features = c("SPARCL1", "EPAS1", "EMP1", "VWF", "SERPINE1"))
plotDimReduceFeature(dim1 = tsne[, 1],dim2 = tsne[, 2],counts = new_lung_resection,features = c("GNLY", "NKG7", "GZMA", "FGFBP2"))
plotDimReduceFeature(dim1 = tsne[, 1],dim2 = tsne[, 2],counts = new_lung_resection,features = c("HLA-DPA1", "HLA-DQA1", "HLA-DQB1", "GPR183", "FCER1A"))
plotDimReduceFeature(dim1 = tsne[, 1],dim2 = tsne[, 2],counts = new_lung_resection,features = c("DCN", "COL1A2", "LUM", "COL3A1", "A2M"))
plotDimReduceFeature(dim1 = tsne[, 1],dim2 = tsne[, 2],counts = new_lung_resection,features = c("CCL21", "TFF3", "GNG11", "MMRN1", "LYVE1"))
plotDimReduceFeature(dim1 = tsne[, 1],dim2 = tsne[, 2],counts = new_lung_resection,features = c("FCGR3B","S100A9","CSF3R","S100A8","IFITM2"))

plotDimReduceFeature(dim1 = tsne[, 1],dim2 = tsne[, 2],counts = new_lung_resection,features = c("CMTM2","MNDA","VNN2","CXCR2","S100P"))
plotDimReduceFeature(dim1 = tsne[, 1],dim2 = tsne[, 2],counts = new_lung_resection,features = c("S100A12","RGS2","FPR1","NAMPT","G0S2","IFITM3"))
DE <- array()

cluster1 <- which(rownames(new_lung_resection) %in% head(diffexpClust1$Gene,20))
Indlist <- c()
for (ct in cell_type){
  y <- differentialExpression(counts = new_lung_resection,
                              celdaMod = goodceldaModel,
                              c1 = ct,
                              c2 = NULL)
  y <- y[y$FDR < 0.05 & abs(y$Log2_FC) > 2, ]
  Indlist <- append(Indlist,which(rownames(new_lung_resection) %in% head(y$Gene,5)))
}

Indexlist <- which(rownames(new_lung_resection) %in% head(cluster1$Gene,5))
Indexlist <- append(Indexlist,which(rownames(new_lung_resection) %in% head(cluster35$Gene,5)))
Indexlist2 <- which((rownames(new_lung_resection) %in% c("APOC1", "MARCO", "C1QB", "C19orf59", "FABP4",
                                                         "SFTPB", "SFTPC", "SFTPA2", "SLPI", "NAPSA",
                                                         "IGHA1", "IGHM", "IGHG3", "IGHG1", "IGKC",
                                                         "CCL5", "TRBC2", "CD2", "IL32", "TRAC",
                                                         "S100A9", "VCAN", "APOBEC3A", "CD300E", "THBS1",
                                                         "SCGB1A1", "SCGB3A1", "SCGB3A2", "BPIFB1", "WFDC2",
                                                         "TPSAB1", "CPA3", "HPGDS", "KIT", "TPSD1",
                                                         "PRG4", "MT2A", "MT1E", "KRT19", "SLC6A14",
                                                         "C9orf24", "TPPP3", "C20orf85", "RSPH1", "C11orf88",
                                                         "EMP2", "AGER", "CAV1", "RTKN2", "CEACAM6",
                                                         "SPARCL1", "EPAS1", "EMP1", "VWF", "SERPINE1",
                                                         "GNLY", "NKG7", "GZMA", "FGFBP2",
                                                         "HLA-DPA1", "HLA-DQA1", "HLA-DQB1", "GPR183", "FCER1A",
                                                         "DCN", "COL1A2", "LUM", "COL3A1", "A2M",
                                                         "CCL21", "TFF3", "GNG11", "MMRN1", "LYVE1")))

normCounts <- normalizeCounts(counts = new_lung_resection, scaleFun = scale)
cell_type <- c("Macrophages", "Type 2", "Neutrophils", "B cell", "T cell", "Type 1", 
               "Ciliated", "Club", "Endothelium","NK cell", "Mast cell", "Transformed epithelium", 
               "Dendritic cell", "Fibroblast", "Lymphatic")
plotHeatmap(counts = normCounts[, clusters(goodceldaModel)$z %in% cell_type],
            z = clusters(goodceldaModel)$z[clusters(goodceldaModel)$z %in% cell_type],
            featureIx = Indlist,
            clusterFeature = FALSE,
            clusterCell = TRUE,
            showNamesFeature = TRUE)
plotHeatmap(counts = normCounts[, clusters(goodceldaModel)$z %in% cell_type],
            clusterCell = TRUE,
            clusterFeature = FALSE,
            z = clusters(goodceldaModel)$z[clusters(goodceldaModel)$z %in% cell_type],
            y = clusters(goodceldaModel)$y,
            featureIx = Indexlist2,
            showNamesFeature = TRUE)
goodceldaModel@clusters$z[goodceldaModel@clusters$z=='27'] <- "26"
goodceldaModel@clusters$z[goodceldaModel@clusters$z=='28'] <- "26"
goodceldaModel@clusters$z[goodceldaModel@clusters$z=='9'] <- "8"
goodceldaModel@clusters$z[goodceldaModel@clusters$z=='10'] <- "8"
goodceldaModel@clusters$z[goodceldaModel@clusters$z=='29'] <- "7"
goodceldaModel@clusters$z[goodceldaModel@clusters$z=='5'] <- "4"
goodceldaModel@clusters$z[goodceldaModel@clusters$z=='6'] <- "4"
goodceldaModel@clusters$z[goodceldaModel@clusters$z=='23'] <- "22"
goodceldaModel@clusters$z[goodceldaModel@clusters$z=='21'] <- "20"
goodceldaModel@clusters$z[goodceldaModel@clusters$z=='26'] <- "Macrophages"
goodceldaModel@clusters$z[goodceldaModel@clusters$z=='8'] <- "Type 2"
goodceldaModel@clusters$z[goodceldaModel@clusters$z=='7'] <- "B cell"
goodceldaModel@clusters$z[goodceldaModel@clusters$z=='18'] <- "T cell"
goodceldaModel@clusters$z[goodceldaModel@clusters$z=='24'] <- "Neutrophils"
goodceldaModel@clusters$z[goodceldaModel@clusters$z=='4'] <- "Club"
goodceldaModel@clusters$z[goodceldaModel@clusters$z=='22'] <- "Mast cell"
goodceldaModel@clusters$z[goodceldaModel@clusters$z=='11'] <- "Transformed epithelium"
goodceldaModel@clusters$z[goodceldaModel@clusters$z=='11'] <- "Ciliated"
goodceldaModel@clusters$z[goodceldaModel@clusters$z=='20'] <- "Type 1"
goodceldaModel@clusters$z[goodceldaModel@clusters$z=='15'] <- "Endothelium"
goodceldaModel@clusters$z[goodceldaModel@clusters$z=='19'] <- "NK cell"
goodceldaModel@clusters$z[goodceldaModel@clusters$z=='30'] <- "Dendritic cell"
goodceldaModel@clusters$z[goodceldaModel@clusters$z=='13'] <- "Lymphatic"
goodceldaModel@clusters$z[goodceldaModel@clusters$z=='1'] <- "B cell"
goodceldaModel@clusters$z[goodceldaModel@clusters$z=='12'] <- "Club"
goodceldaModel@clusters$z[goodceldaModel@clusters$z=='14'] <- "B cell"
goodceldaModel@clusters$z[goodceldaModel@clusters$z=='2'] <- "B cell"
goodceldaModel@clusters$z[goodceldaModel@clusters$z=='3'] <- "B cell"
goodceldaModel@clusters$z[goodceldaModel@clusters$z=='25'] <- "Neutrophils"

length(unique(goodceldaModel@clusters$z))
