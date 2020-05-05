#Question 1
Q1data <- read.csv("project_sib_pheno_and_RV_data(1).csv")

#Heritability based on sibling pair
sibdata <- subset(Q1data,is.na(GRIP2)==FALSE)
kidgrip = c(sibdata$GRIP1,sibdata$GRIP2)
meankidgrip = mean(kidgrip,na.rm=T)
sdkidgrip = sd(kidgrip,na.rm=T)
N = nrow(sibdata)
adjkid1 = sibdata$GRIP1 - meankidgrip
adjkid2 = sibdata$GRIP2 - meankidgrip
icc = sum(adjkid1*adjkid2)/sdkidgrip^2/(N-1)
h2 = 2*icc
print(h2) #0.2575


rm(list = ls())

#Question 2
#..\plink\plink --file project --out project                 #get .bed and .bim file
#..\plink\plink --bfile project --freq --out project_freq    #get minor allele frequency data
#..\plink\plink --bfile project --freqx --out project_freqx  #get genotype frequency data
#..\plink\plink --bfile project --hardy --out project_hwe    #get hardy-weinburg report
Q2data <- read.table("project_hwe.hwe",header = TRUE)
num.failed <- sum(Q2data$P<0.05 & Q2data$TEST=="ALL")
failed.SNP <- Q2data[Q2data$P<0.05 & Q2data$TEST=="ALL",]
Alltest <- subset(Q2data,TEST=="ALL")
write.csv(Alltest,"HWE_for_16SNPs.csv")

rm(list = ls())

#Question 3
#..\plink\plink --file project --out project                 #get minor allele frequency data
Q3data <- read.table("project_freq.frq",header = TRUE)
write.csv(Q3data,"MAF_for_16SNPs.csv")


rm(list = ls())

#Question 4
#a: 
0.05/16 #0.003125

#b: 
0.05/16 #0.003125


#Question 5
#a csv file available.The effect size in standard deviation units is calculated as Effect/11.3 where Effect is 
#the effect size in combined part in Table 1 of the paper.
#b
Q5data <- read.csv("Adjusted_effect_size.csv",header = TRUE)
Q5data <- Q5data[order(Q5data$SNP),]
Q5af <- read.table("project_freq.frq",header = TRUE)
Q5af <- Q5af[order(Q5af$SNP),]
Q5af$Major <- 1-Q5af$MAF
Q5data$minor <- Q5af$MAF
Q5data$major <- Q5af$Major
Q5data$h2 <- 2*Q5af$MAF*Q5af$Major*(Q5data$Adjusted.Effect.size)^2
write.csv(Q5data,"locus_specific_h2.csv")

#c
h2 <- max(Q5data$h2)
power.h2<-function(n,h2,df,sig.level=0.05){
  lambda<-n*h2
  print(paste("non-centrality parameter=",lambda,sep=""))
  cv<-qf(1-sig.level,df,n-df)  
  pow<-pf(cv,df,n-df,lambda,lower.tail=FALSE)
  return(pow)}
power.h2(22000,h2,1,0.05) #0.3734289


rm(list = ls())

#Question 6
#a:
Q6data <- read.csv("project_sib_pheno_and_RV_data(1).csv")
K <- sum(Q6data$FRAC1==2) / nrow(Q6data)
n_case <- sum(Q6data$FRAC1==2)
n_control <- nrow(Q6data)-n_case
Q6raf <- read.table("project_freq.frq",header = TRUE)
Q6raf$Major <- 1-Q6raf$MAF
Q6raf$RAF <- c(0.70400,0.16390,0.40580,0.20860,0.55090,0.65740,0.01109,0.57650,0.38030,0.92689,0.04420,0.60400,
               0.20830,1.00000,0.53220,0.23000)
Q6raf$Risk_allele <- c("G","C","C","A","A","T","A","A","T","A","A","C","A","A","A","A")
##Get f0,f1,f2
Q6genotype <- read.table("project.ped")
Q6SNP <- read.table("project.map")
Q6SNP$Risk_allele <- c("A","A","C","A","C","G","A","T","A","T","A","A","C","A","A","A")
Q6SNP$Other_allele <- c("C","G","T","G","T","A","G","C","G","C","G","T","T","T","G","G")
penetrance <- c()
for(i in seq(1,16)){
  affect_status <- Q6genotype$V6
  allele1 <- Q6genotype[,i*2+5]
  allele2 <- Q6genotype[,i*2+6]
  Genotype <- paste0(allele1,allele2)
  dat <- as.data.frame(cbind(affect_status,Genotype))
  risk_al <- Q6SNP$Risk_allele[i]
  other_al <- Q6SNP$Other_allele[i]
  ###get f2
  two <- paste0(risk_al,risk_al)
  n2 <- subset(dat,Genotype==two)
  f2 <- sum(n2$affect_status==2)/nrow(n2)
  ###get f0
  zero <- paste0(other_al,other_al)
  n0 <- subset(dat,Genotype==zero)
  f0 <- sum(n0$affect_status==2)/nrow(n0)
  ###get f1
  n1 <- dat[dat$Genotype!=zero&dat$Genotype!=two,]
  f1 <- sum(n1$affect_status==2)/nrow(n1)
  penetrance <- rbind(penetrance,c(f0,f1,f2))
}
penetrance <- as.data.frame(penetrance)
colnames(penetrance) <- c("f0","f1","f2")
GRR <- cbind(Q6SNP,penetrance)
GRR$ratio1 <- GRR$f1/GRR$f0
GRR$ratio2 <- GRR$f2/GRR$f0
colnames(GRR)[1:4] <- c("chr","SNP","Genetic_dist","position")
rownames(GRR) <- c(1:16)
GRR <- GRR[order(GRR$SNP),]
Q6raf <- Q6raf[order(Q6raf$SNP),]
data_for_power <- cbind(Q6raf[,c(2,8)],GRR[,c(10,11)],replicate(16,K),replicate(16,n_case),replicate(16,n_control))
#write.csv(data_for_power,"Q6Power.csv")

rm(list=ls())

#Question 7
#a
#..\plink\plink --bfile project --linear --covar project.covar --pheno project.pheno --out linear_add
p <- 0.05/16
Q7data <- read.table("linear_add.assoc.linear",header = TRUE)
test <- subset(Q7data,TEST=="ADD")
write.csv(test,"Assoc_grip.csv")
example_snp <- subset(Q7data,SNP=="rs958685")
write.csv(example_snp,"Assoc_grip_example.csv")


rm(list=ls())

#Question 8
#..\plink\plink --bfile project --logistic beta --covar project.covar  --pheno project.pheno --pheno-name FRAC 
#--out logistic_add_all_covar
p <- 0.05/16
Q8data <- read.table("logistic_add_all_covar.assoc.logistic",header = TRUE)
test <- subset(Q8data,TEST=="ADD")
write.csv(test,"Assoc_fracture.csv")
example_snp <- subset(Q8data,SNP=="rs72762373")
write.csv(example_snp,"Assoc_fracture_example.csv")

rm(list=ls())

#Question 9
Q9data <- read.csv("project_sib_pheno_and_RV_data(1).csv",as.is = T)
Q9RV <- Q9data[,c(1,3,4,5,6,7,8,14,15,16,17,18,19,20,21,22)]
Q9RV <- subset(Q9RV,is.na(RV1)==FALSE)
##MB all
maf<-apply(Q9RV[,8:16],2,function(i)sum(i)/(2*length(i)))
weights<-1/sqrt(1000*maf*(1-maf))
xxx<-sapply(1:9,function(i)Q9RV[,i+7]*weights[i])
Q9RV$MB<-apply(xxx,1,sum)
mb.all<-glm(GRIP1~sex1+age1+pop+bmi1+hgt1+MB,data=Q9RV)
summary(mb.all)
##CMC all
Q9RV$CMC<-apply(Q9RV[,8:16],1,sum)
CMC.all<-glm(GRIP1~sex1+age1+pop+bmi1+hgt1+CMC,data=Q9RV)
summary(CMC.all)
