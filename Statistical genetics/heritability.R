famdata <- read.csv("HW3_ldl2020.csv")

print(summary(famdata))

print(famdata[1:20,])

famdata$paveldl = (famdata$momldl + famdata$dadldl)/2

#Heritability based on parent offspring

print(lm(kid1ldl ~ momldl,data=famdata))
print(lm(kid1ldl ~ dadldl,data=famdata))
print(lm(kid1ldl ~ paveldl,data=famdata))

print(summary(lm(kid2ldl ~ momldl,data=famdata)))
print(summary(lm(kid2ldl ~ dadldl,data=famdata)))
print(summary(lm(kid2ldl ~ paveldl,data=famdata)))

#Heritability based on sibling pair

kidldl = c(famdata$kid1ldl,famdata$kid2ldl)
meankidldl = mean(kidldl,na.rm=T)
sdkidldl = sd(kidldl,na.rm=T)
N = nrow(famdata)

adjkid1 = famdata$kid1ldl - meankidldl
adjkid2 = famdata$kid2ldl - meankidldl
icc = sum(adjkid1*adjkid2)/sdkidldl^2/(N-1)

h2 = 2*icc

print(h2)

#Heritability based on siblings (ANOVA)
kid1 = famdata[,c("famid","kid1id","kid1ldl")]
kid2 = famdata[,c("famid","kid2id","kid2ldl")]
names(kid1) <- c("famid","kidid","ldl")
names(kid2) <- c("famid","kidid","ldl")
kids <- rbind(kid1,kid2)

print(summary(aov(ldl ~ as.factor(famid),data=kids)))
print(2*(1106-654)/(1106+(2-1)*654))
