#Change directory first

famdata <- read.csv("HW1_data.csv",as.is=T,na=".")

print(famdata[1:20,])
print(summary(famdata))

Nind <- nrow(famdata)
print("Number of indivdiduals")
print(Nind)

Nfam <- length(table(famdata$famid))
print("Number of families")
print(Nfam)

Nchildren <- sum(famdata$mo !=0)
print("Number of children")
print(Nchildren)

N.nonmiss.pheno <- sum(famdata$aff!=0)
print("Number with non missing phenotype")
print(N.nonmiss.pheno)

prop.miss.pheno <- sum(famdata$aff==0)/Nind
print("Proportion with missing phenotype")
print(prop.miss.pheno)

N.nonmiss.pheno.kid <- sum(famdata$aff!=0 & famdata$mo!=0)
print("Number of children with non missing phenotype")
print(N.nonmiss.pheno.kid)

prop.miss.pheno.kid <- sum(famdata$aff==0 & famdata$mo!=0)/Nchildren
print("Proportion of children with missing phenotype")
print(prop.miss.pheno.kid)

N.DNA <- sum(!is.na(famdata$DNAdate))
print("Number with DNA")
print(N.DNA)

prop.noDNA <- sum(is.na(famdata$DNAdate))/Nind
print("Proportion without DNA")
print(prop.noDNA)

prop.kid.noDNA <- sum(is.na(famdata$DNAdate) & famdata$mo!=0)/sum(famdata$mo!=0)
print("Proportion children without DNA")
print(prop.kid.noDNA)

N.both.DNA.pheno <- sum(!is.na(famdata$DNAdate) & famdata$aff!=0)
print("Number with both DNA and phenotype")
print(N.both.DNA.pheno)

Nsibs <- tapply(famdata$mo!=0,famdata$famid,sum)
mean.sibship <- mean(Nsibs)
print("Average sibship size")
print(mean.sibship)

Nsibs.DNA.pheno <- tapply(famdata$mo!=0 & famdata$aff!=0 & !is.na(famdata$DNAdate),
famdata$famid,sum)
mean.sibship.DNA.pheno <- mean(Nsibs.DNA.pheno)
print("Average sibship size (with DNA and phenotype available)")
print(mean.sibship.DNA.pheno)

Nkid.aff <- tapply(famdata$mo!=0 & famdata$aff==2,famdata$famid,sum)
print("Number of affected children per families")
print(table(Nkid.aff))

N.genotype.available <- sum(!is.na(famdata$rs12075_Al1) & !is.na(famdata$rs12075_Al2))
print(N.genotype.available)
prop.genotype.available <- N.genotype.available / Nind

Nfounders <- sum(famdata$fa==0)
print(Nfounders)
prop.miss.geno.parents <- sum(is.na(famdata$rs12075_Al1) & is.na(famdata$rs12075_Al2) & famdata$fa == 0) / Nfounders
print("Proportion of founders with missing genotype")
print(prop.miss.geno.parents)

Nfoundermothers <- sum(famdata$fa==0 & famdata$sex==2)
print(Nfoundermothers)
prop.miss.geno.mothers <- sum(is.na(famdata$rs12075_Al1) & is.na(famdata$rs12075_Al2) & famdata$fa == 0 & famdata$sex==2)/Nfoundermothers
print(prop.miss.geno.mothers)
prop.miss.geno.children <- sum(is.na(famdata$rs12075_Al1) & is.na(famdata$rs12075_Al2) & famdata$fa != 0) / Nchildren
print(prop.miss.geno.children)
Nwith.geno.and.pheo <- sum(famdata$aff!=0 & !is.na(famdata$rs12075_Al1) & !is.na(famdata$rs12075_Al2))
print(Nwith.geno.and.pheo)

Nsibs.only.pheogeno <- tapply(famdata$mo!=0 & famdata$aff!=0 & !is.na(famdata$rs12075_Al1) & !is.na(famdata$rs12075_Al2),famdata$famid,sum)
mean.sibship.only.pheogeno <- mean(Nsibs.only.pheogeno)
print("Average sibship size")
print(mean.sibship.only.pheogeno)

Nparent.known.pheo <- sum(famdata$mo==0 & famdata$aff!=0)
prop.Nparents.aff <- sum(famdata$mo==0 & famdata$aff==2) / Nparent.known.pheo
print(prop.Nparents.aff)

snp.genotypes<-ifelse(famdata$rs12075_Al1<famdata$rs12075_Al2,
                      paste(famdata$rs12075_Al1,famdata$rs12075_Al2,sep="/"),
                      paste(famdata$rs12075_Al2,famdata$rs12075_Al1,sep="/"))

print(table(snp.genotypes))

snp.alleles <- c(famdata$rs12075_Al1,famdata$rs12075_Al2)
table(snp.alleles)