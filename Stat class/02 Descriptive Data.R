library(rockchalk)
data(religioncrime)
head(religioncrime, n=5)
mean(religioncrime)
religioncrime
mean(religioncrime[ ,"heaven"])
mean(religioncrime[ ,"hell"])
mean(religioncrime[ ,"crime"])
religioncrime$halfheaven <- ifelse(religioncrime$heaven>=0.5, "yes", "no")
religioncrime$halfhell <- ifelse(religioncrime$hell>=0.5, "yes", "no")
religioncrime.table <-table(religioncrime$halfheaven)
religioncrime.table
sum(religioncrime.table)
religioncrime.table/sum(religioncrime.table)
religioncrime.table <- table(Heaven = religioncrime$heaven, Hell = religioncrime$hell)
religioncrime.table
addmargins(religioncrime.table)
prop.table(religioncrime.table)
table(religioncrime$heaven, religioncrime$hell)
half.table <- table(religioncrime$halfheaven, religioncrime$halfhell)
addmargins(table)
half.table
addmargins(half.table)
prop.table(half.table)
religioncrime.table <-table(BelieveInHeaven=religioncrime$halfheaven, BelieveInHell=religioncrime$halfhell)
religioncrime.table
pnorm(3, mean=2, sd=5, lower.tail=TRUE, log.p=FALSE)
dnorm(3, mean=2, sd=5, log=FALSE)
a<-1/sqrt(2*pi*5^2) 
b<-exp(-((3-2)^2)/(2*5^2))
?pnorm
pnorm(4,mean=2,sd=5,lower.tail=TRUE, log.p=FALSE)-pnorm(2,mean=2,sd=5,lower.tail=TRUE,log.p=FALSE)
