library(ggplot2)
library(reshape2)
library(DESeq2)
#4
sample1 <- read.table("count_file/SRR1177981.txt",header = TRUE)[,c(1,7)]
sample2 <- read.table("count_file/SRR1177982.txt",header = TRUE)[,c(1,7)]
sample3 <- read.table("count_file/SRR1177983.txt",header = TRUE)[,c(1,7)]
sample4 <- read.table("count_file/SRR1178008.txt",header = TRUE)[,c(1,7)]
sample5 <- read.table("count_file/SRR1178009.txt",header = TRUE)[,c(1,7)]
sample6 <- read.table("count_file/SRR1178010.txt",header = TRUE)[,c(1,7)]
sample7 <- read.table("count_file/SRR1178014.txt",header = TRUE)[,c(1,7)]
sample8 <- read.table("count_file/SRR1178021.txt",header = TRUE)[,c(1,7)]
sample9 <- read.table("count_file/SRR1178047.txt",header = TRUE)[,c(1,7)]

clean <- function(x,sp){
  rown <- x[,1]
  coln <- sp
  dataset <- as.data.frame(x[,-1])
  rownames(dataset) <- rown
  colnames(dataset) <- sp
  return(dataset)
}

sample1 <- clean(sample1,"SRR1177981")
sample2 <- clean(sample2,"SRR1177982")
sample3 <- clean(sample3,"SRR1177983")
sample4 <- clean(sample4,"SRR1178008")
sample5 <- clean(sample5,"SRR1178009")
sample6 <- clean(sample6,"SRR1178010")
sample7 <- clean(sample7,"SRR1178014")
sample8 <- clean(sample8,"SRR1178021")
sample9 <- clean(sample9,"SRR1178047")
dat <- cbind(sample1,sample2,sample3,sample4,sample5,sample6,sample7,sample8,
             sample9)
#write.csv(dat,"concat_counts.csv")
gene.id <- rownames(dat)
dat.for.plot <- cbind(gene.id,dat)
rownames(dat.for.plot) <- c(1:length(dat.for.plot[,1]))
dat.for.plot <- melt(dat.for.plot,value.name = "counts")
dat.for.plot <- subset(dat.for.plot,counts>300&counts<10000)
dat.for.plot$counts <- log(dat.for.plot$counts)
ggplot(data = dat.for.plot,aes(x=variable,y=counts)) + 
  geom_boxplot(aes(fill=variable))

#4.1
meta.info <- read.csv("toxgroup_3_rna_info.csv",sep = ",")
controls <- read.csv("/project/bf528/project_3/samples/control_counts.csv",
                     sep = ",")
rown.ct <- controls[,1]
controls <- controls[,c("SRR1178050","SRR1178061","SRR1178063","SRR1178004",
                        "SRR1178006","SRR1178013")]
rownames(controls) <- rown.ct
dat <- cbind(dat,controls)

#4.2 4.3
dat1 <- subset(dat,rowSums(dat)!=0)
##LEF
LEF <- dat1[,c("SRR1178008","SRR1178009","SRR1178010",
               "SRR1178050","SRR1178061","SRR1178063")]
LEF.info <- meta.info[c(1,2,3,10,11,12),] 
LEF.dds <- DESeqDataSetFromMatrix(
  countData = LEF,
  colData = LEF.info,
  design= ~ mode_of_action
)
LEF.dds$mode_of_action <- relevel(LEF.dds$mode_of_action, ref='Control')
LEF.dds <- DESeq(LEF.dds)
LEF.res <- results(LEF.dds, contrast=c('mode_of_action','AhR','Control'))
LEF.res <- lfcShrink(LEF.dds, coef=2)
LEF.res <- LEF.res[order(LEF.res$padj),]
LEF.res.sig <- subset(LEF.res,padj<0.05)
LEF.res.top10 <- head(LEF.res,10)
LEF.for.plot <- as.data.frame(LEF.res)
ggplot(data = LEF.for.plot,aes(x = log2FoldChange))+geom_histogram()
ggplot(data = LEF.for.plot,aes(x = log2FoldChange,y=pvalue))+ geom_point()
ggplot(data = LEF.for.plot,aes(x = pvalue,y=log2FoldChange))+ geom_point()
#write.csv(LEF.res,"LEF.csv")
##FLU
FLU <- dat1[,c("SRR1178014","SRR1178021","SRR1178047",
               "SRR1178050","SRR1178061","SRR1178063")]
FLU.info <- meta.info[c(4,5,6,10,11,12),] 
FLU.dds <- DESeqDataSetFromMatrix(
  countData = FLU,
  colData = FLU.info,
  design= ~ mode_of_action
)
FLU.dds$mode_of_action <- relevel(FLU.dds$mode_of_action, ref='Control')
FLU.dds <- DESeq(FLU.dds)
FLU.res <- results(FLU.dds, contrast=c('mode_of_action','CAR/PXR','Control'))
FLU.res <- lfcShrink(FLU.dds, coef=2)
FLU.res <- FLU.res[order(FLU.res$padj),]
FLU.res.sig <- subset(FLU.res,padj<0.05)
FLU.res.top10 <- head(FLU.res,10)
FLU.for.plot <- as.data.frame(FLU.res)
ggplot(data = FLU.for.plot,aes(x = log2FoldChange))+geom_histogram()
ggplot(data = FLU.for.plot,aes(x = log2FoldChange,y=pvalue))+ geom_point()
#ggplot(data = FLU.for.plot,aes(x = pvalue,y=log2FoldChange))+ geom_point()
#write.csv(FLU.res,"FLU.csv")

##IFO
IFO <- dat1[,c("SRR1177981","SRR1177982","SRR1177983",
               "SRR1178004","SRR1178006","SRR1178013")]
IFO.info <- meta.info[c(7,8,9,13,14,15),] 
IFO.dds <- DESeqDataSetFromMatrix(
  countData = IFO,
  colData = IFO.info,
  design= ~ mode_of_action
)
IFO.dds$mode_of_action <- relevel(IFO.dds$mode_of_action, ref='Control')
IFO.dds <- DESeq(IFO.dds)
IFO.res <- results(IFO.dds, contrast=c('mode_of_action','DNA_Damage','Control'))
IFO.res <- lfcShrink(IFO.dds, coef=2)
IFO.res <- IFO.res[order(IFO.res$padj),]
IFO.res.sig <- subset(IFO.res,padj<0.05)
IFO.res.top10 <- head(IFO.res,10)
IFO.for.plot <- as.data.frame(IFO.res)
ggplot(data = IFO.for.plot,aes(x = log2FoldChange))+geom_histogram()
ggplot(data = IFO.for.plot,aes(x = log2FoldChange,y=pvalue))+ geom_point()
#ggplot(data = FLU.for.plot,aes(x = pvalue,y=log2FoldChange))+ geom_point()
#write.csv(IFO.res,"IFO.csv")

##normalized counts
all.dds <- DESeqDataSetFromMatrix(
  countData = dat1,
  colData = meta.info,
  design= ~ mode_of_action
)
all.dds$mode_of_action <- relevel(all.dds$mode_of_action, ref='Control')
all.dds <- estimateSizeFactors(all.dds)
normed.dat1 <- counts(all.dds,normalized=TRUE)
#write.csv(normed.dat1,"normed_counts.csv")

##heatmap
norm.dat <- read.csv("normed_counts.csv",row.names = 1)
LEF <- read.csv("LEF.csv",row.names = 1)
FLU <- read.csv("FLU.csv",row.names = 1)
IFO <- read.csv("IFO.csv",row.names = 1)
genelist <- c(head(rownames(LEF),1000),head(rownames(FLU),1000),
              head(rownames(IFO),1000))
#variance <- apply(norm.dat, 1, sd)
#norm.dat <- norm.dat[order(variance),]
#top1000 <- tail(norm.dat,1000)
norm.dat1 <- norm.dat[genelist,]
heatmap(as.matrix(norm.dat1))

LEF_control <- norm.dat[head(rownames(LEF),1000),c(4,5,6,10,11,12)]
heatmap(as.matrix(LEF_control))
FLU_control <- norm.dat[head(rownames(FLU),1000),c(7,8,9,10,11,12)]
heatmap(as.matrix(FLU_control))
IFO_control <- norm.dat[head(rownames(IFO),1000),c(1,2,3,13,14,15)]
heatmap(as.matrix(IFO_control))
