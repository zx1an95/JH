library(limma)
rma <- read.table('/projectnb/bf528/project_3/samples/liver-normalization-rma.txt',
                  sep='\t',
                  as.is=TRUE,
                  header=TRUE,
                  row.names=1,
)
meta.info <- read.csv("group_3_mic_info.csv")
dat <- rma[,paste0("X",meta.info$array_id)]
##LEF
LEF.id <- as.character(meta.info[meta.info$chemical=="LEFLUNOMIDE",]$array_id)
LEF.control.id <- as.character(meta.info[meta.info$vehicle=="CORN_OIL_100_%"&meta.info$chemical=="Control",]$array_id)
LEF <- dat[,paste0("X",c(LEF.id,LEF.control.id))]
LEF.meta <- meta.info[(meta.info$vehicle=="CORN_OIL_100_%"&
                        meta.info$chemical=="Control")|
                        (meta.info$chemical=="LEFLUNOMIDE"),]
LEF.design <- model.matrix(
  ~factor(
    LEF.meta$chemical,
    levels=c('Control','LEFLUNOMIDE')
  )
)
colnames(LEF.design) <- c('Intercept','LEFLUNOMIDE')
LEF.fit <- lmFit(LEF, LEF.design)
LEF.fit <- eBayes(LEF.fit)
LEF.t <- topTable(LEF.fit, coef=2, n=nrow(LEF), adjust='BH')
#write.csv(LEF.t,'LEF_microarray.csv')
LEF.t.sig <- subset(LEF.t,adj.P.Val<0.05)
LEF.t.top10 <- head(LEF.t,10)
LEF.for.plot <- as.data.frame(LEF.t)
ggplot(data = LEF.for.plot,aes(x = logFC))+geom_histogram()
ggplot(data = LEF.for.plot,aes(x = logFC,y=P.Value))+ geom_point()

##FLU
##LEF
FLU.id <- as.character(meta.info[meta.info$chemical=="FLUCONAZOLE",]$array_id)
FLU.control.id <- as.character(meta.info[meta.info$vehicle=="CORN_OIL_100_%"&meta.info$chemical=="Control",]$array_id)
FLU <- dat[,paste0("X",c(FLU.id,FLU.control.id))]
FLU.meta <- meta.info[(meta.info$vehicle=="CORN_OIL_100_%"&
                         meta.info$chemical=="Control")|
                        (meta.info$chemical=="FLUCONAZOLE"),]
FLU.design <- model.matrix(
  ~factor(
    FLU.meta$chemical,
    levels=c('Control','FLUCONAZOLE')
  )
)
colnames(FLU.design) <- c('Intercept','FLUCONAZOLE')
FLU.fit <- lmFit(FLU, FLU.design)
FLU.fit <- eBayes(FLU.fit)
FLU.t <- topTable(FLU.fit, coef=2, n=nrow(FLU), adjust='BH')
#write.csv(FLU.t,'FLU_microarray.csv')
FLU.t.sig <- subset(FLU.t,adj.P.Val<0.05)
FLU.t.top10 <- head(FLU.t,10)
FLU.for.plot <- as.data.frame(FLU.t)
ggplot(data = FLU.for.plot,aes(x = logFC))+geom_histogram()
ggplot(data = FLU.for.plot,aes(x = logFC,y=P.Value))+ geom_point()

##IFO
IFO.id <- as.character(meta.info[meta.info$chemical=="IFOSFAMIDE",]$array_id)
IFO.control.id <- as.character(meta.info[meta.info$vehicle=="SALINE_100_%"&meta.info$chemical=="Control",]$array_id)
IFO <- dat[,paste0("X",c(IFO.id,IFO.control.id))]
IFO.meta <- meta.info[(meta.info$vehicle=="SALINE_100_%"&
                         meta.info$chemical=="Control")|
                        (meta.info$chemical=="IFOSFAMIDE"),]
IFO.design <- model.matrix(
  ~factor(
    IFO.meta$chemical,
    levels=c('Control','IFOSFAMIDE')
  )
)
colnames(IFO.design) <- c('Intercept','IFOSFAMIDE')
IFO.fit <- lmFit(IFO, IFO.design)
IFO.fit <- eBayes(IFO.fit)
IFO.t <- topTable(IFO.fit, coef=2, n=nrow(IFO), adjust='BH')
#write.csv(IFO.t,'IFO_microarray.csv')
IFO.t.sig <- subset(IFO.t,adj.P.Val<0.05)
IFO.t.top10 <- head(IFO.t,10)
IFO.for.plot <- as.data.frame(IFO.t)
ggplot(data = IFO.for.plot,aes(x = logFC))+geom_histogram()
ggplot(data = IFO.for.plot,aes(x = logFC,y=P.Value))+ geom_point()

##heatmap
norm.dat <- read.csv("../pg/normed_counts.csv",row.names = 1)
variance <- apply(norm.dat, 1, sd)
norm.dat <- norm.dat[order(variance),]
top1000 <- tail(norm.dat,1000)
heatmap(as.matrix(top1000))
