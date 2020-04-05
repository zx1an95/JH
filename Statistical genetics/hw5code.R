setwd("D:/BU/Statistical genetics")
dat <- read.table("cc_model.model",header = TRUE)
Allele.test <- dat[dat$TEST=="ALLELIC",]
Allele.sort.by.p <- Allele.test[sort(Allele.test$P,index.return=TRUE)$ix,]
Geno.test <- dat[dat$TEST=="GENO",]
Geno.sort.by.p <- Geno.test[sort(Geno.test$P,index.return=TRUE,na.last = TRUE)$ix,]

#2
dat2 <- read.table("logistic.add.beta.assoc.logistic",header = TRUE)
d2.sort.by.p <- dat2[sort(dat2$P,index.return=TRUE,na.last = TRUE)$ix,]
dat2[dat2$SNP=="rs4866676",]

#3
dat3 <- read.table("logistic.hethom.assoc.logistic",header = TRUE)
test.dat3 <- dat3[dat3$TEST=="GENO_2DF",]
d3.sort.by.p <- test.dat3[sort(test.dat3$P,index.return=TRUE,na.last = TRUE)$ix,]
dat3[dat3$SNP=="rs4866676",]

#4
dat4 <- read.table("logistic.adjust.assoc.logistic",header = TRUE)
dat4[dat4$SNP=="rs4866676",]

#5
dat5 <- read.table("linear.general.assoc.linear",header = TRUE)
d5.test <- dat5[dat5$TEST=="GENO_2DF",]
d5.sort.by.p <- d5.test[sort(d5.test$P,index.return=TRUE,na.last = TRUE)$ix,]
dat5[dat5$SNP=="rs4866676",]
