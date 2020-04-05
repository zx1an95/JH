one <- read.table("medelian2020.txt",header=T,as.is=T)

dim(one)
names(one)

table(one$FORMAT)

table(one$IND01)

##Missingness rate
length(one[which(one$IND01 == "./."),"IND01"])/length(one[,"IND01"])
length(one[which(one$IND02 == "./."),"IND02"])/length(one[,"IND02"])
length(one[which(one$IND03 == "./."),"IND03"])/length(one[,"IND03"])
length(one[which(one$IND04 == "./."),"IND04"])/length(one[,"IND04"])
length(one[which(one$IND05 == "./."),"IND05"])/length(one[,"IND05"])
##Number of variant sites
vs1 <- one[which(one$IND01 %in% c("1/1","0/1","1/0")),]
length(vs1[,"IND01"])
vs2 <- one[which(one$IND02 %in% c("1/1","0/1","1/0")),]
length(vs2[,"IND02"])
vs3 <- one[which(one$IND03 %in% c("1/1","0/1","1/0")),]
length(vs3[,"IND03"])
vs4 <- one[which(one$IND04 %in% c("1/1","0/1","1/0")),]
length(vs4[,"IND04"])
vs5 <- one[which(one$IND05 %in% c("1/1","0/1","1/0")),]
length(vs5[,"IND05"])
##Number of singletons
AC.count <- function(x){
	y <- ifelse(x == "1/1",2,NA)
	y <- ifelse(x %in% c("0/1","1/0"),1,y)
	y <- ifelse(x %in% c("0/0"),0,y)
	sum(y,na.rm=T)
}

one$AC <- apply(one[,paste("IND0",seq(1,5,1),sep="")],1,AC.count)

length(one[which(one$IND01 %in% c("0/1","1/0") & one$AC == 1),"IND01"])
length(one[which(one$IND02 %in% c("0/1","1/0") & one$AC == 1),"IND02"])
length(one[which(one$IND03 %in% c("0/1","1/0") & one$AC == 1),"IND03"])
length(one[which(one$IND04 %in% c("0/1","1/0") & one$AC == 1),"IND04"])
length(one[which(one$IND05 %in% c("0/1","1/0") & one$AC == 1),"IND05"])

##Ti/Tv ratio
vsnps1 <- vs1[which(vs1$REF %in% c("A","C","G","T") & vs1$ALT %in% c("A","C","G","T")),]
vsnps2 <- vs2[which(vs2$REF %in% c("A","C","G","T") & vs2$ALT %in% c("A","C","G","T")),]
vsnps3 <- vs3[which(vs3$REF %in% c("A","C","G","T") & vs3$ALT %in% c("A","C","G","T")),]
vsnps4 <- vs4[which(vs4$REF %in% c("A","C","G","T") & vs4$ALT %in% c("A","C","G","T")),]
vsnps5 <- vs5[which(vs5$REF %in% c("A","C","G","T") & vs5$ALT %in% c("A","C","G","T")),]

vsnps1$TiTv <- ifelse((vsnps1$REF == 'A' & vsnps1$ALT == 'G') | 
			   (vsnps1$REF == 'G' & vsnps1$ALT == 'A') | 
			   (vsnps1$REF == 'T' & vsnps1$ALT == 'C') | 
			   (vsnps1$REF == 'C' & vsnps1$ALT == 'T'),"Ti","Tv")
vsnps2$TiTv <- ifelse((vsnps2$REF == 'A' & vsnps2$ALT == 'G') | 
                        (vsnps2$REF == 'G' & vsnps2$ALT == 'A') | 
                        (vsnps2$REF == 'T' & vsnps2$ALT == 'C') | 
                        (vsnps2$REF == 'C' & vsnps2$ALT == 'T'),"Ti","Tv")
vsnps3$TiTv <- ifelse((vsnps3$REF == 'A' & vsnps3$ALT == 'G') | 
                        (vsnps3$REF == 'G' & vsnps3$ALT == 'A') | 
                        (vsnps3$REF == 'T' & vsnps3$ALT == 'C') | 
                        (vsnps3$REF == 'C' & vsnps3$ALT == 'T'),"Ti","Tv")
vsnps4$TiTv <- ifelse((vsnps4$REF == 'A' & vsnps4$ALT == 'G') | 
                        (vsnps4$REF == 'G' & vsnps4$ALT == 'A') | 
                        (vsnps4$REF == 'T' & vsnps4$ALT == 'C') | 
                        (vsnps4$REF == 'C' & vsnps4$ALT == 'T'),"Ti","Tv")
vsnps5$TiTv <- ifelse((vsnps5$REF == 'A' & vsnps5$ALT == 'G') | 
                        (vsnps5$REF == 'G' & vsnps5$ALT == 'A') | 
                        (vsnps5$REF == 'T' & vsnps5$ALT == 'C') | 
                        (vsnps5$REF == 'C' & vsnps5$ALT == 'T'),"Ti","Tv")

length(vsnps1[which(vsnps1$TiTv == "Ti"),"TiTv"])/length(vsnps1[which(vsnps1$TiTv == "Tv"),"TiTv"])
length(vsnps2[which(vsnps2$TiTv == "Ti"),"TiTv"])/length(vsnps2[which(vsnps2$TiTv == "Tv"),"TiTv"])
length(vsnps3[which(vsnps3$TiTv == "Ti"),"TiTv"])/length(vsnps3[which(vsnps3$TiTv == "Tv"),"TiTv"])
length(vsnps4[which(vsnps4$TiTv == "Ti"),"TiTv"])/length(vsnps4[which(vsnps4$TiTv == "Tv"),"TiTv"])
length(vsnps5[which(vsnps5$TiTv == "Ti"),"TiTv"])/length(vsnps5[which(vsnps5$TiTv == "Tv"),"TiTv"])

#check <- one[,c(10,11,12,13,14)]
check <- one[one$IND01=="0/0"&one$IND05=="0/0"&
                 one$IND02 %in% c("1/1","0/1","1/0")&
                 one$IND03 %in% c("1/1","0/1","1/0")&
                 one$IND04 %in% c("1/1","0/1","1/0"),]
check_info <- strsplit(check$INFO,split = "[|]")
protein_coding <- sapply(check_info,function(x) x[8]=="protein_coding")
check_pc <- check[protein_coding,]
check[c(15,21),]
