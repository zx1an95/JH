setwd("/Users/adminuser/Desktop/Work/CPSC_Courses/CPSC440/Fall_2018/Homeworks/") 

#Read in the data
the.data <- read.csv("HW4.Q3.csv", head = FALSE)

#Conduct the Shapiro-Wilk test for normality
shapiro.test(t(the.data))

#Make a box plot. This will be saved into your working directory
pdf("Boxplot_HW4_Q3.pdf")
  boxplot(the.data, col = "Red", 
          ylab = "Distortion Level")
  points(mean(t(the.data)), pch = 3, cex = 0.75) #Adds the mean, which is indicated by the "+" symbol
dev.off()


#Make a QQ-plot
pdf("QQplot_HW4_Q3.pdf")
  qqnorm(t(the.data), col = "Blue")
dev.off()



