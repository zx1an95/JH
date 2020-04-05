#Set your working directory

setwd("/Users/adminuser/Desktop/Work/CPSC_Courses/CPSC440/Spring_2017/Homeworks/") #Type in the path to your working directory here


#################### Homework 5, Problem 1 ################
#Read in the data
the.data <- read.csv("HW5.Q1.csv", head = TRUE)

#Separate the data into two subsets: one for LactNumb = 1, and one for LactNumb = 3
Subaru <- the.data[which(the.data[,1] == "Subaru"),]
Volkswagen <- the.data[which(the.data[,1] == "Volkswagen"),]


#Conduct the Shapiro-Wilk test for normality
shapiro.test(t(Subaru[,2]))
shapiro.test(t(Volkswagen[,2]))

#Calculate the means so that we can include them in the box plot
means <- by(the.data[,2], the.data[,1], mean) 


t.test(t(Subaru[,2]), t(Volkswagen[,2]), alternative="less", var.equal=TRUE)


t.test(t(Subaru[,2]), t(Volkswagen[,2]), alternative="less", var.equal=FALSE)

#Make a box plot. This will be saved into your working directory
pdf("Boxplot_HW5_Q1.pdf")
  boxplot(the.data[,2]~the.data[,1], col = "Red", 
          ylab = "Quarter mile (seconds)", xlab = "1 = Subaru, 2 = Volkswagen )")
  points(c(1:length(means)), means, pch = 3, cex = 0.75) #Adds the means, which are indicated by the "+" symbol
dev.off()


#Make a QQ-plot
pdf("QQplot_HW5_Q1.pdf", width = 12)
  par(mfrow = c(1,2))
  qqnorm(t(Subaru[,2]), col = "Blue", main =  "Subaru")
  qqnorm(t(Volkswagen[,2]), col = "Blue", main =  "Volkswagen")
dev.off()




#################### Homework 5, Problem 2 ################
#Clear your R workspace of the objects created from problem 1
# Rationale: we don't want to accidentally use the data from
# problem 1 to answer problem 2
rm(list = ls())
setwd("/Users/adminuser/Desktop/Work/CPSC_Courses/CPSC440/Spring_2017/Homeworks/") #Type in the path to your working directory here

#Read in the data
the.data <- read.csv("HW5.Q2.csv", head = TRUE)

#Separate the data into two subsets: 
Two.times <- the.data[,2]
Four.times <- the.data[,3]


#Conduct the Shapiro-Wilk test for normality
shapiro.test(t(Two.times))
shapiro.test(t(Four.times))

#Calculate the means so that we can include them in the box plot
means <- c(mean(Two.times), mean(Four.times))

#Reorganize the data so that they can be read into the boxplot() function
the.data.for.boxplot <- data.frame(rbind(cbind(rep(1, length(Two.times)),  Two.times) ,   
                                         cbind(rep(2, length(Four.times)),  Four.times)) )


#Make a box plot. This will be saved into your working directory
pdf("Boxplot_HW5_Q2.pdf")
  par(mar = c(3,6,5,1), mfrow = c(1,1))
  boxplot(the.data.for.boxplot[,2]~the.data.for.boxplot[,1], col = "Red", 
          axes = FALSE, ylim = c(100,400), ylab = "Milk Production (kg)", cex.lab = 2.0)
  points(c(1:length(means)), means, pch = 3, cex = 0.75) #Adds the means, which are indicated by the "+" symbol
  axis(1, at=c(1,2),cex.axis=1.5,labels=c("Twice a day", "Four times a day"),tick=FALSE)
  axis(2, at=c(100,200,300,400), cex.axis=1.5, labels = c("100","200","300","400"), tick=FALSE)
  box()
dev.off()


#Make a QQ-plot
pdf("QQplot_HW5_Q2.pdf", width = 12)
  par(mfrow = c(1,2))
  qqnorm(t(the.data[,2]), col = "Blue", main =  "Twice a day")
  qqnorm(t(the.data[,3]), col = "Blue", main =  "Four times a day")
dev.off()

t.test(t(the.data[,3]), t(the.data[,2]), paired = TRUE, alternative = "greater")
