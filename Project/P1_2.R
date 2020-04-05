#1
BiocManager::install(c("affy","affyPLM","sva","AnnotationDbi","hgu133plus2.db"))
#2
library(affy)
library(affyPLM)
library(sva)
library(AnnotationDbi)
library(hgu133plus2.db)
setwd("/projectnb/bf528/users/group2/project1/samples")
#3
dat <- ReadAffy()
norm.dat <- rma(dat)
#4
RLE.and.NUSE <- fitPLM(dat,normalize = TRUE,background = TRUE)
medianRLE.for.samples <- RLE(RLE.and.NUSE,type="stats")
medianNUSE.for.samples <- NUSE(RLE.and.NUSE,type="stats")
hist(medianRLE.for.samples[1,],col="blue",nclass=8,main="Histogram of median RLE",xlab="median RLE")
hist(medianNUSE.for.samples[1,],col="blue",nclass=8,main="Histogram of median NUSE",xlab="median NUSE")
#5
anno_data <- read.csv("/project/bf528/project_1/doc/proj_metadata.csv")
new.dat <- exprs(norm.dat)
new.dat <- ComBat(new.dat,batch = anno_data$normalizationcombatbatch,mod = anno_data$normalizationcombatmod)
write.csv(new.dat,file = "/projectnb/bf528/users/group2/project1/programmer_deliverables/Step5data.csv")
test.data <- read.csv("/projectnb/bf528/users/group2/project1/programmer_deliverables/Step5data.csv")
#6
new.dat <- t(scale(t(new.dat)))
PCA.data <- prcomp(new.dat,scale. = FALSE,center = FALSE)
summary.data <- summary(PCA.data)
PC.for.plot <- PCA.data$rotation[,1:2]
plot(PC.for.plot,col="red",xlab=paste0("PC1:",summary.data$importance[2,1]),ylab=paste0("PC2:",summary.data$importance[2,2]),
     main="PCA Plot")
ggplot(PC.for.plot,aes(x = PC1, y = PC2))+geom_point()+labs(title = "PCA plot")+
  xlab(paste0("PC1:",summary.data$importance[2,1]))+
  ylab(paste0("PC2:",summary.data$importance[2,2]))
#visualization for different batch
batch <- anno_data$normalizationcombatbatch
library(ggplot2)
PC.for.plot <- as.data.frame(PC.for.plot)
ggplot(PC.for.plot,aes(x = PC1, y = PC2))+geom_point(aes(col=batch))+
  labs(title="PCA plot")+xlab(paste0("PC1:",summary.data$importance[2,1]))+
  ylab(paste0("PC2:",summary.data$importance[2,2]))


#Analysts part1
P1.data <- read.csv("/projectnb/bf528/users/group2/project1/programmer_deliverables/programmer_deliverables/Step5data.csv")
rownames(P1.data) <- P1.data$X
P1.data <- P1.data[,-1]
n <- floor(134*0.2)
filter1 <- P1.data > log2(15)
filter1 <- apply(filter1,1,sum)>n
P1.data <- P1.data[filter1,]
gene.variance <- apply(P1.data,1,var)
median.variance <- median(gene.variance)
num_measurement <- ncol(P1.data)
test.stat <- (num_measurement-1)*(gene.variance/median.variance) 
#critical.values1 <- qchisq(0.005,num_measurement-1)#if two tail
critical.values2 <- qchisq(0.99,num_measurement-1)
filter2 <- test.stat > critical.values2 
Step2.data <- P1.data[filter2,]
gene.std <- apply(Step2.data,1,sd)
gene.mean <- apply(Step2.data,1,mean)
gene.cv <- gene.std / gene.mean
filter3 <- gene.cv > 0.186
Step3.data <- Step2.data[filter3,]
#write.csv(Step3.data,"Analyst_Step4.csv")
#write.csv(Step2.data,"Analyst_Step5.csv")

#Analyst part2
#Step3.data <- read.csv("Analyst_Step4.csv")
anno_data <- read.csv("/project/bf528/project_1/doc/proj_metadata.csv")
hc <- hclust(dist(t(Step3.data)))
clus = cutree(hc, 2)
group_list= as.factor(clus)
table(group_list)#1:57 2:77
color.for.heatmap <- ifelse(anno_data$cit.coloncancermolecularsubtype=="C3","red","blue")
#save heatmap
heatmap(as.matrix(Step3.data),ColSideColors = color.for.heatmap)
#t-test
ttest.results <- apply(Step3.data, 1, function(x) t.test(x = x[clus==1], y = x[clus==2]))
probeset.ID <- names(ttest.results)
t.stats <- sapply(ttest.results,function(x) x$statistic)
p_values <- sapply(ttest.results,function(x) x$p.value)
p_adjusted <- p.adjust(p_values,method = "fdr")
DE.genes <- data.frame(t_stats=t.stats,p_value=p_values,adjusted_p_values=p_adjusted)
rownames(DE.genes) <- probeset.ID
sig.DE.genes <- DE.genes[DE.genes$adjusted_p_values<0.05,]
num_of_sig_genes <- nrow(sig.DE.genes)#1236
sort.DE.genes <- DE.genes[sort(DE.genes$adjusted_p_values,index.return = TRUE)$ix,]
sort.sig.DE.genes <- sig.DE.genes[sort(sig.DE.genes$adjusted_p_values,index.return = TRUE)$ix,]
#write.csv(DE.genes,"Differentially_expressed_genes.csv")
Best.described <- head(sort.DE.genes,20)

#function to get pos and neg markers
Diff_exprs <- function(dataset,clusterlist,pos_or_neg,c1,c2){
  markers <- apply(dataset, 1, function(x) t.test(x = x[clusterlist==c1], y = x[clusterlist==c2],alternative = pos_or_neg))
  probeset.ID <- names(markers)
  t.stats <- sapply(markers,function(x) x$statistic)
  p_values <- sapply(markers,function(x) x$p.value)
  p_adjusted <- p.adjust(p_values,method = "fdr")
  marker.genes <- data.frame(t_stats=t.stats,p_value=p_values,adjusted_p_values=p_adjusted)
  rownames(marker.genes) <- probeset.ID
  return(marker.genes)
}
#positive and negative marker for 2 clusters
cluster1.pos.marker <- Diff_exprs(Step3.data,clus,"greater",1,2)
cluster1.neg.marker <- Diff_exprs(Step3.data,clus,"less",1,2)
cluster2.pos.marker <- Diff_exprs(Step3.data,clus,"greater",2,1)
cluster2.neg.marker <- Diff_exprs(Step3.data,clus,"less",2,1)
cluster1.pos.marker <- cluster1.pos.marker[sort(cluster1.pos.marker$adjusted_p_values,index.return = TRUE)$ix,]
cluster1.neg.marker <- cluster1.neg.marker[sort(cluster1.neg.marker$adjusted_p_values,index.return = TRUE)$ix,]
cluster2.pos.marker <- cluster2.pos.marker[sort(cluster2.pos.marker$adjusted_p_values,index.return = TRUE)$ix,]
cluster2.neg.marker <- cluster2.neg.marker[sort(cluster2.neg.marker$adjusted_p_values,index.return = TRUE)$ix,]
length(cluster1.pos.marker[cluster1.pos.marker$adjusted_p_values<0.05,]$t_stats)#919
length(cluster2.neg.marker[cluster1.neg.marker$adjusted_p_values<0.05,]$t_stats)#313
length(cluster2.pos.marker[cluster2.pos.marker$adjusted_p_values<0.05,]$t_stats)#313
length(cluster2.neg.marker[cluster2.neg.marker$adjusted_p_values<0.05,]$t_stats)#919
head(cluster1.pos.marker,10) #best described using pos markers
head(cluster2.pos.marker,10) #best described using pos markers


#For biologists files(not useful as two sided test)
hc.for.bio <- hclust(dist(t(Step2.data)))
clus.for.bio = cutree(hc.for.bio, 2)
group_list.for.bio = as.factor(clus.for.bio)
table(group_list.for.bio)
ttest.results.for.bio <- apply(Step2.data, 1, function(x) t.test(x = x[clus.for.bio==1], y = x[clus.for.bio==2]))
probeset.ID.for.bio <- names(ttest.results.for.bio)
t.stats.for.bio <- sapply(ttest.results.for.bio,function(x) x$statistic)
p_values.for.bio <- sapply(ttest.results.for.bio,function(x) x$p.value)
p_adjusted.for.bio <- p.adjust(p_values.for.bio,method = "fdr")
DE.genes.for.bio <- data.frame(t_stats=t.stats.for.bio,p_value=p_values.for.bio,adjusted_p_values=p_adjusted.for.bio)
rownames(DE.genes.for.bio) <- probeset.ID.for.bio
sig.DE.genes.for.bio <- DE.genes.for.bio[DE.genes.for.bio$adjusted_p_values<0.05,]
num_of_sig_genes.for.bio <- nrow(sig.DE.genes.for.bio)
sort.DE.genes.for.bio <- DE.genes.for.bio[sort(DE.genes.for.bio$adjusted_p_values,index.return = TRUE)$ix,]
sort.sig.DE.genes.for.bio <- sig.DE.genes.for.bio[sort(sig.DE.genes.for.bio$adjusted_p_values,index.return = TRUE)$ix,]
#write.csv(DE.genes.for.bio,"differentially_expressed_genes_for_biologists.csv")

#pos and neg marker gene for 2 clusters(useful)
hc.for.bio <- hclust(dist(t(Step2.data)))
clus.for.bio <- cutree(hc.for.bio, 2)
group_list.for.bio <- as.factor(clus.for.bio)
cluster1.pos.marker.for.bio <- Diff_exprs(Step2.data,clus.for.bio,"greater",1,2)
cluster2.pos.marker.for.bio <- Diff_exprs(Step2.data,clus.for.bio,"greater",2,1)
cluster1.neg.marker.for.bio <- Diff_exprs(Step2.data,clus.for.bio,"less",1,2)
cluster2.neg.marker.for.bio <- Diff_exprs(Step2.data,clus.for.bio,"less",2,1)
#cluster1.pos.marker.for.bio <- cluster1.pos.marker.for.bio[sort(cluster1.pos.marker.for.bio$adjusted_p_values,index.return = TRUE)$ix,]
#cluster2.pos.marker.for.bio <- cluster2.pos.marker.for.bio[sort(cluster2.pos.marker.for.bio$adjusted_p_values,index.return = TRUE)$ix,]
#length(cluster1.pos.marker.for.bio[cluster1.pos.marker.for.bio$adjusted_p_values<0.05,]$t_stats)
#length(cluster2.pos.marker.for.bio[cluster2.pos.marker.for.bio$adjusted_p_values<0.05,]$t_stats)
#head(cluster1.pos.marker,10) #best described using pos markers
#head(cluster2.pos.marker,10) #best described using pos markers

#Biologists part
#example.results <- read.csv("/project/bf528/project_1/data/differential_expression_results.csv")
gene.symbol <- select(hgu133plus2.db,keys = rownames(cluster1.pos.marker.for.bio),column = "SYMBOL")
#sum(rownames(DE.genes.for.bio) %in% gene.symbol$PROBEID) See if all mapped
#gene.symbol[duplicated(gene.symbol$PROBEID),] #one probe ID to mutiple symbol or mutiple symbol to one gene
gene.symbol.rm.dup <- gene.symbol[!duplicated(gene.symbol$PROBEID),]
#length(unique(gene.symbol.rm.dup[!is.na(gene.symbol.rm.dup$SYMBOL),]$SYMBOL)) check unique symbol
maxp <- by(Step2.data,gene.symbol.rm.dup$SYMBOL,function(x) rownames(x)[which.max(rowMeans(x))])
uniprobes <- as.character(maxp)
gene.symbol.selected <- gene.symbol.rm.dup[gene.symbol.rm.dup$PROBEID %in% uniprobes,]
cluster1.pos.marker.for.bio$symbol <- gene.symbol.rm.dup$SYMBOL
cluster2.pos.marker.for.bio$symbol <- gene.symbol.rm.dup$SYMBOL
cluster1.neg.marker.for.bio$symbol <- gene.symbol.rm.dup$SYMBOL
cluster2.neg.marker.for.bio$symbol <- gene.symbol.rm.dup$SYMBOL
cluster1.pos.with.eff.symbol <- cluster1.pos.marker.for.bio[uniprobes,]
cluster2.pos.with.eff.symbol <- cluster2.pos.marker.for.bio[uniprobes,]
cluster1.neg.with.eff.symbol <- cluster1.neg.marker.for.bio[uniprobes,]
cluster2.neg.with.eff.symbol <- cluster2.neg.marker.for.bio[uniprobes,]
sorted.pos.cluster1 <- cluster1.pos.with.eff.symbol[sort(cluster1.pos.with.eff.symbol$adjusted_p_values,index.return = TRUE)$ix,]
sorted.pos.cluster2 <- cluster2.pos.with.eff.symbol[sort(cluster2.pos.with.eff.symbol$adjusted_p_values,index.return = TRUE)$ix,]
sorted.neg.cluster1 <- cluster1.neg.with.eff.symbol[sort(cluster1.neg.with.eff.symbol$adjusted_p_values,index.return = TRUE)$ix,]
sorted.neg.cluster2 <- cluster2.neg.with.eff.symbol[sort(cluster2.neg.with.eff.symbol$adjusted_p_values,index.return = TRUE)$ix,]
c1top1000up <- head(sorted.pos.cluster1,1000)
c1top1000down <- head(sorted.neg.cluster1,1000)
c2top1000up <- head(sorted.pos.cluster2,1000)
c2top1000down <- head(sorted.neg.cluster2,1000)
c1top10up <- head(sorted.pos.cluster1,10)
c1top10down <- head(sorted.neg.cluster1,10)
c2top10up <- head(sorted.pos.cluster2,10)
c2top10down <- head(sorted.neg.cluster2,10)
#write.csv(c1top10up,"top10up.csv")
#write.csv(c1top10down,"top10down.csv")

#GSEABase
library(GSEABase)
KEGG <- getGmt("c2.cp.kegg.v7.0.symbols.gmt")#186
GO <- getGmt("c5.all.v7.0.symbols.gmt")#9996
Hall <- getGmt("h.all.v7.0.symbols.gmt")#50
filter.for.sig.c1 <- sorted.pos.cluster1$adjusted_p_values < 0.05 
filter.for.sig.c2 <- sorted.pos.cluster2$adjusted_p_values < 0.05 
#sig.DE.pos.genes.c1 <- sorted.pos.cluster1[filter.for.sig.c1,]#use for enrichment analysis
#sig.DE.pos.genes.c2 <- sorted.pos.cluster2[filter.for.sig.c2,]
geneID.c1 <- sorted.pos.cluster1$symbol
geneID.c2 <- sorted.pos.cluster2$symbol
terms.and.genes.KEGG <- sapply(KEGG@.Data,function(x) x@geneIds)
names.KEGG <- sapply(KEGG@.Data,function(x) x@setName) 
names(terms.and.genes.KEGG) <- names.KEGG
terms.and.genes.Hall <- sapply(Hall@.Data,function(x) x@geneIds)
names.Hall <- sapply(Hall@.Data,function(x) x@setName) 
names(terms.and.genes.Hall) <- names.Hall
terms.and.genes.GO <- sapply(GO@.Data,function(x) x@geneIds)
names.GO <- sapply(GO@.Data,function(x) x@setName) 
names(terms.and.genes.GO) <- names.GO
#Fisher exact test
enrichment.analysis <- function(genesets){
  DE.in.genesets <- sum(geneID.c1[filter.for.sig.c1] %in% genesets)
  DE.notin.genesets <- sum(!(geneID.c1[filter.for.sig.c1] %in% genesets))
  NDE.in.genesets <- sum(geneID.c1[!filter.for.sig.c1] %in% genesets)
  NDE.notin.genesets <- sum(!(geneID.c1[!filter.for.sig.c1] %in% genesets))
  results <- c(DE.in.genesets,DE.notin.genesets,NDE.in.genesets,NDE.notin.genesets)
  test.result <- fisher.test(matrix(results,nrow=2))
  return(c(test.result$p.value,test.result$estimate))
}
Hall.result.c1 <- sapply(terms.and.genes.Hall,enrichment.analysis)
KEGG.result.c1 <- sapply(terms.and.genes.KEGG,enrichment.analysis)
GO.result.c1 <- sapply(terms.and.genes.GO,enrichment.analysis)

#for c2
enrichment.analysis2 <- function(genesets){
  DE.in.genesets <- sum(geneID.c2[filter.for.sig.c2] %in% genesets)
  DE.notin.genesets <- sum(!(geneID.c2[filter.for.sig.c2] %in% genesets))
  NDE.in.genesets <- sum(geneID.c2[!filter.for.sig.c2] %in% genesets)
  NDE.notin.genesets <- sum(!(geneID.c2[!filter.for.sig.c2] %in% genesets))
  results <- c(DE.in.genesets,DE.notin.genesets,NDE.in.genesets,NDE.notin.genesets)
  test.result <- fisher.test(matrix(results,nrow=2))
  return(c(test.result$p.value,test.result$estimate))
}
Hall.result.c2 <- sapply(terms.and.genes.Hall,enrichment.analysis2)
KEGG.result.c2 <- sapply(terms.and.genes.KEGG,enrichment.analysis2)
GO.result.c2 <- sapply(terms.and.genes.GO,enrichment.analysis2)
#creat data frame for each cluster and geneset type
get.enrich.df <- function(results){
  return(data.frame(stats_estimate=results[2,],
                    p_value=results[1,],
                    p_adjusted=p.adjust(results[1,],method = "fdr")))
}
Hall.enrichment.c1 <- get.enrich.df(Hall.result.c1)
KEGG.enrichment.c1 <- get.enrich.df(KEGG.result.c1)
GO.enrichment.c1 <- get.enrich.df(GO.result.c1)
Hall.enrichment.c2 <- get.enrich.df(Hall.result.c2)
KEGG.enrichment.c2 <- get.enrich.df(KEGG.result.c2)
GO.enrichment.c2 <- get.enrich.df(GO.result.c2)
sorted.Hall.c1 <- Hall.enrichment.c1[sort(Hall.enrichment.c1$p_adjusted,index.return = TRUE)$ix,]
sorted.Hall.c2 <- Hall.enrichment.c2[sort(Hall.enrichment.c2$p_adjusted,index.return = TRUE)$ix,]
sorted.KEGG.c1 <- KEGG.enrichment.c1[sort(KEGG.enrichment.c1$p_adjusted,index.return = TRUE)$ix,]
sorted.KEGG.c2 <- KEGG.enrichment.c2[sort(KEGG.enrichment.c2$p_adjusted,index.return = TRUE)$ix,]
sorted.GO.c1 <- GO.enrichment.c1[sort(GO.enrichment.c1$p_adjusted,index.return = TRUE)$ix,]
sorted.GO.c2 <- GO.enrichment.c2[sort(GO.enrichment.c2$p_adjusted,index.return = TRUE)$ix,]
sig.Hall.c1 <- sorted.Hall.c1[sorted.Hall.c1$p_adjusted<0.05,]#30
sig.Hall.c2 <- sorted.Hall.c2[sorted.Hall.c2$p_adjusted<0.05,]#33
sig.KEGG.c1 <- sorted.KEGG.c1[sorted.KEGG.c1$p_adjusted<0.05,]#68
sig.KEGG.c2 <- sorted.KEGG.c2[sorted.KEGG.c2$p_adjusted<0.05,]#85
sig.GO.c1 <- sorted.GO.c1[sorted.GO.c1$p_adjusted<0.05,]#1619
sig.GO.c2 <- sorted.GO.c2[sorted.GO.c2$p_adjusted<0.05,]#1001

top3Hallc1 <- head(sig.Hall.c1,3)
top3Hallc2 <- head(sig.Hall.c2,3)
top3KEGGc1 <- head(sig.KEGG.c1,3)
top3KEGGc2 <- head(sig.KEGG.c2,3)
top3GOc1 <- head(sig.GO.c1,3)
top3GOc2 <- head(sig.GO.c2,3)

#write.csv(sorted.Hall.c1,"Hall_c1.csv")
#write.csv(sorted.Hall.c2,"Hall_c2.csv")
#write.csv(sorted.KEGG.c1,"KEGG_c1.csv")
#write.csv(sorted.KEGG.c2,"KEGG_c2.csv")
#write.csv(sorted.GO.c1,"GO_c1.csv")
#write.csv(sorted.GO.c2,"GO_c2.csv")