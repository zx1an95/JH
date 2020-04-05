#6.1

dat <- read.table("/projectnb/bf528/users/group2/project2/analyst/gene_exp.diff",header = TRUE)
sort.dat <- dat[order(dat$q_value),]
top10.DE <- head(sort.dat[,c(3,8,9,10,12,13)],10)
write.csv(top10.DE,"top10_DE.csv")
#6.2
hist(dat$log2.fold_change.,col = "blue",xlab = "Log2 fold change",
     nclass = 12,xlim = c(-6,6),main = "Histogram of Log2 fold change for all genes")
#The number of the negative values is more than the number of positive values. 

#6.3
sig.dat <- subset(dat,significant=="yes")#5188
sig.q001.genes <- subset(dat,q_value<0.01) #3788
sig.p001.genes <- subset(dat,p_value<0.01) #4686
up.regulated.q001 <- subset(sig.q001.genes,log2.fold_change.>0) #2101
down.regulated.q001 <- subset(sig.q001.genes,log2.fold_change.<0) #1687
up.regulated.p001 <- subset(sig.p001.genes,log2.fold_change.>0) #2532
down.regulated.p001 <- subset(sig.p001.genes,log2.fold_change.<0) #2154

#6.4
hist(sig.dat$log2.fold_change.,col = "blue",xlab = "Log2 fold change",
     xlim = c(-10,10),main = "Histogram of Log2 fold change from significant DE genes")


#6.5
up.regulated <- subset(sig.dat,log2.fold_change.>0)#2757
down.regulated <- subset(sig.dat,log2.fold_change.<0)#2431

#6.6
#write.csv(up.regulated$gene,"up-regulated_gene.csv")
#write.csv(down.regulated$gene,"down-regulated_gene.csv")

#6.7
  