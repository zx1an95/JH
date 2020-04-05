dat <- read.table("lm_adj.assoc.logistic.adjusted",header = TRUE)
sum(dat$BONF>dat$FDR_BH)
sum(dat$SIDAK_SS>dat$FDR_BH)
head(dat[order(dat$UNADJ),c(2,3)],5)
head(dat[order(dat$BONF),c(2,5)],5)
head(dat[order(dat$SIDAK_SS),c(2,7)],5)
head(dat[order(dat$FDR_BH),c(2,9)],5)
tests <- dat$UNADJ * 112
eff.m.bonf <- data.frame(SNP=dat$SNP,BONF=dat$UNADJ * 55)
