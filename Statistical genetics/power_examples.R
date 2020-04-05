 
#power.prop.test syntax:
#power.prop.test(n = NULL, p1 = NULL, p2 = NULL, sig.level = 0.05, power = NULL, 
#alternative = c("two.sided", "one.sided"), strict = FALSE) 
##note: assumes equal numbers of cases and controls



##computint power using power.prop.test:  note, I am using the proportions we determined
##using the grr=2, dominant model, risk allele frequency 0.2, population prevalence 0.10.
##(the genotype proportions for cases and controls are on slide 22 of the notes).
##For other genetic models, you'd need to determine the proportions using the formulas in 
##the class notes.

##allele based test:
#power for class example, 200 alleles (100 people) per group:
power.prop.test(n=200,p1=0.2941,p2=0.1895,sig.level=0.05,alternative=c("two.sided"))
#sample size required for 80% power, sig level 0.05:  class example
##note: assumes equal numbers of cases and controls
power.prop.test(p1=0.2941,p2=0.1895,sig.level=0.05,power=.80,alternative=c("two.sided"))

#dominant model, 100 genotypes per group
power.prop.test(n=100,p1=0.3412,p2=0.5294,sig.level=0.05,alternative=c("two.sided"))
#sample size required for 80% power, sig level=0.05, dominant model (class example)
##note: assumes equal numbers of cases and controls
power.prop.test(p1=0.3412,p2=0.5294,sig.level=0.05,power=0.80,alternative=c("two.sided"))


###################
##In class quantitative trait example.
#N=200:  what is the critical value of an F distribution with 2 and 198 df?
cv<-qf(1-0.05,1,4999)
cv
#power for the in class example:
pow<-pf(cv,1,4999,0,lower.tail=FALSE)
pow
##pf(q, df1, df2, ncp, lower.tail = TRUE, log.p = FALSE) ##syntax

##more general:  function to determine power for given sample size and
##locus-specific heritability, using test with df degrees of freedom:

power.h2<-function(n,h2,df,sig.level=0.05){
 #power for a specified locus-specific heritability
 #note: locus h2 is a function of allele frequency and genotype means
 lambda<-n*h2
 print(paste("non-centrality parameter=",lambda,sep=""))
 cv<-qf(1-sig.level,df,n-df)  #critical value for given sample size and sig level
 pow<-pf(cv,df,n-df,lambda,lower.tail=FALSE)
 return(pow)}

##testing out the function:  
power.h2(200,.01,.05,df=2)  ## gives power == 0.223
power.h2(200,.01,.05,df=1)  ## gives power == 0.291 (1df test rather than 2df test above)
power.h2(5000,0,df=1,sig.level=0.05) ## gives 80% power --> we need h2=0.049 to have 80% power for n=200

