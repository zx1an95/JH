source("https://bioconductor.org/biocLite.R")
options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
BiocManager::install(c("ALL","CLL", "pasilla", "airway"))
BiocManager::install(c("limma","DESeq2", "clusterProfiler"))
install.packages("reshape2")
install.packages("ggplot2")
library(CLL)
library(limma)
library(DESeq2)
library(clusterProfiler)
library(reshape2)
library(dplyr)
data(sCLLex)
sCLLex
e <- exprs(sCLLex)
str(e)
head(e)
dim(e)
colnames(e)
rownames(e)
sampleNames(sCLLex)
varMetadata(sCLLex)
featureNames(sCLLex)[1:100]
featureNames(sCLLex) %>% unique() %>% length()
pdata <- pData(sCLLex)
group_list <- as.character(pdata[,2])
table(group_list)
par(cex = 0.7)
n.sample <- ncol(e)
if(n.sample>40) par(cex = 0.5)
cols <- rainbow(n.sample*1.2)
boxplot(e, col = cols,main="expression value",las=2)
BiocManager::install("hgu95av2.db")
library(hgu95av2.db)
ls("package:hgu95av2.db")
capture.output(hgu95av2())
as.list(hgu95av2SYMBOL)[1]
C <- hgu95av2CHR
Llength(C)
Rlength(C)
Rkeys(C)
Rkeys(C) <- c("6","8")
Rkeys(C)
as.list(C)
table(toTable(C)[,2])
s <- hgu95av2SYMBOL
summary(s)
mapped_probe <- mappedkeys(s)
Lkeys(s)[!Lkeys(s) %in% mapped_probe]
ss <- as.list(s[mapped_probe]) #as.list(s)[mapped_probe]
colnames(sst)
sst[grep("^TP53$",sst$symbol),]
unique(sst$symbol) %>% length()
table(sst$symbol) %>% sort() %>% tail()
table(sst$symbol) %>%table()
e1 <- e[rownames(e)%in%mapped_probe,]
maxp <- by(e1,sst$symbol,function(x) rownames(x)[which.max(rowMeans(x))])
uniprobes <- as.character(maxp)
efilt <- e[rownames(e)%in%uniprobes,]
rownames(efilt) <- sst[match(rownames(efilt),sst$probe_id),2]
library(reshape2)
dim(efilt)
m_efilt <- melt(efilt)
colnames(m_efilt) <- c('symbol','sample','value')
m_efilt$group <- rep(group_list,each=nrow(efilt))
head(efilt)
head(m_efilt)
e_mean <- tail(sort(apply(efilt,1,mean)),30)
e_median <- tail(sort(apply(efilt,1,median)), 30)
e_max <- tail(sort(apply(efilt,1,max)),30)
e_min <- tail(sort(apply(efilt,1,min)),30)
e_sd <- tail(sort(apply(efilt,1,sd)),30)
e_var <- tail(sort(apply(efilt,1,var)),30)
e_mad <- tail(sort(apply(efilt,1,mad)),30)
install.packages("UpSetR")
library("UpSetR")
e_all <- c(names(e_mean),names(e_median),names(e_max),names(e_min),
           names(e_sd),names(e_var),names(e_mad)) %>% unique()
edat <- data.frame(e_all,
                e_mean=ifelse(e_all %in% names(e_mean) ,1,0),
                e_median=ifelse(e_all %in% names(e_median) ,1,0),
                e_max=ifelse(e_all %in% names(e_max) ,1,0),
                e_min=ifelse(e_all %in% names(e_min) ,1,0),
                e_sd=ifelse(e_all %in% names(e_sd) ,1,0),
                e_var=ifelse(e_all %in% names(e_var) ,1,0),
                e_mad=ifelse(e_all %in% names(e_mad) ,1,0)
)
upset(edat,nsets = 7,sets.bar.color = "#56B4E9")
gl <- as.factor(group_list)
group1 = which(group_list == levels(gl)[1]) 
group2 = which(group_list == levels(gl)[2])
et1 <- e[, group1]
et2 <- e[, group2]
et <- cbind(et1,et2)
pvals <- apply(e, 1, function(x){
  t.test(as.numeric(x)~group_list)$p.value
})
p.adj <- p.adjust(pvals, method = "BH")
eavg_1 <- rowMeans(et1)
eavg_2 <- rowMeans(et2)
log2FC <- eavg_2-eavg_1
DEG_t.test <- cbind(eavg_1, eavg_2, log2FC, pvals, p.adj)
DEG_t.test <- DEG_t.test[order(DEG_t.test[,4]),]
DEG_t.test <- as.data.frame(DEG_t.test)
head(DEG_t.test)
top30_gene=names(e_mad)
top30_matrix=efilt[top30_gene,] 
top30_matrix=t(scale(t(top30_matrix)))
test <- list(countsdata=top30_matrix,group_list=group_list)
test$group_list <- as.data.frame(test$group_list)
rownames(test$group_list) <- colnames(test$countsdata)
pheatmap(test$countsdata, annotation_col = test$group_list)
suppressMessages(library(limma))
design1 <- model.matrix(~factor(group_list))
colnames(design1) <- levels(factor(group_list))
rownames(design1) <- colnames(e)
fit1 <- lmFit(e,design1)
fit1 <- eBayes(fit1)
options(digits = 3)
mtx1 <- topTable(fit1,coef=2,adjust='BH',n=Inf)
DEG_mtx1 <- na.omit(mtx1)
head(DEG_mtx1)

design2=model.matrix(~0+factor(group_list))
colnames(design2)=levels(factor(group_list))
rownames(design2)=colnames(e)
fit2=lmFit(e,design2)

fit2=contrasts.fit(fit2, contrast.matrix) 
fit2=eBayes(fit2) 

mtx2=topTable(fit2, coef=1, n=Inf) 
DEG_mtx2= na.omit(mtx2)
#write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
head(DEG_mtx2)

logFC_cutoff <- with(DEG,mean(abs(logFC)) + 2*sd(abs(logFC)) ) 
DEG$result = as.factor(ifelse(DEG$P.Value < 0.05 & abs(DEG$logFC) > logFC_cutoff,
                              ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
this_tile <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                    '\nThe number of up gene is ',nrow(DEG[DEG$result =='UP',]) ,
                    '\nThe number of down gene is ',nrow(DEG[DEG$result =='DOWN',])
)
library(ggplot2)
ggplot(data=DEG, aes(x=logFC, y=-log10(adj.P.Val), color=result)) +
  geom_point(alpha=0.4, size=1.75) +
  theme_set(theme_set(theme_bw(base_size=20)))+
  xlab("log2 fold change") + ylab("-log10 p-value") +
  ggtitle( this_tile ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
  scale_colour_manual(values = c('blue','black','red')) 
