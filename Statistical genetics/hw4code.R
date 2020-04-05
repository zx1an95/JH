setwd("D:/BU/Statistical genetics")
dat <- read.table("bs858.spring2020_freq.frq",header = TRUE)
num.maf.0.05 <- sum(dat$MAF<0.05)
dat[dat$MAF<0.05,]
sort.by.maf <- dat[sort(dat$MAF,index.return=TRUE)$ix,]
answer1 <- head(sort.by.maf,10)[,c(2,3,4,5)]
write.csv(answer1,"hw4p4q1.csv")

dat2 <- read.table("bs858.spring2020_freqx.frqx",header = TRUE,sep = "\t")
num.count.15.homA1 <- sum(dat2$C.HOM.A1.<15)
num.count.15.het <- sum(dat2$C.HET.<15)
num.count.15.homA2 <- sum(dat2$C.HOM.A2.<15)
dat2.10 <- dat2[dat2$C.HOM.A1.<10,]
#sort.by.count <- dat2[sort(dat2$C.HOM.A1.,index.return=TRUE)$ix,]
answer2 <- dat2.10[,c(2,3,4,5,6,7)]
write.csv(answer2,"hw4p4q2.csv")

dat3 <- read.table("bs858.spring2020_hwe.hwe",header = TRUE)
num.failed.0.001 <- sum(dat3$P<0.001 & dat3$TEST=="ALL")
failed.0.001 <- dat3[dat3$P<0.001 & dat3$TEST=="ALL",]
failed.0.001.alltest <- dat3[dat3$SNP%in%failed.0.001$SNP,] 
failed.0.001.controls <- failed.0.001.alltest[failed.0.001.alltest$P<0.001 & failed.0.001.alltest$TEST=="UNAFF",]
#sort.by.count <- dat2[sort(dat2$C.HOM.A1.,index.return=TRUE)$ix,]

dat4 <- read.table("bs858.spring2020_ld.ld",header = TRUE)
r2.0.5 <- dat4[dat4$R2>0.5,]
r2.0.5 <- r2.0.5[,c(3,6,7,8)]
r2.0.5$SNP_A_MAF <- c(0.2473,0.2428,0.3430,0.4808)
r2.0.5$SNP_B_MAF <- c(0.2428,0.2308,0.3168,0.4398)
write.csv(r2.0.5,"hw4p4q6.csv")
