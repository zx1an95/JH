
#Read in the sick workers example
setwd("/Users/adminuser/Desktop/Work/CPSC_Courses/CPSC440/Fall_2016/Lecture_Notes/")
sick.worker.data <- read.csv("Chapter_11_Sick_Worker_Data.csv", head = TRUE)


#Make a scatter plot with absences on the y-axis and experience on the x-axis
plot(sick.worker.data$abscences~sick.worker.data$experience)

#Fit a simple linear regression model absences = bo + b1*experience + epsilon
SLR.model <- lm(sick.worker.data$abscences ~ sick.worker.data$experience)


#Look at the output from the model fit
anova(SLR.model)
summary(SLR.model)

#Intermediate step: calcualte the resuldiuals and predicted values (for the observations without missing absences data)
residuals <- SLR.model$residuals
predicted.values <- SLR.model$fitted.values

#Make a residual plot
plot(residuals~predicted.values)
abline(h = 0, lty = 2)

#Make a QQ-plot of the residuals
qqnorm(residuals)


#Make a histogram of the residuals
hist(residuals)

#Conduct a Shapiro-Wilk test for normality of the residuals, which tests Ho: Residuals are normally distributed
shapiro.test(residuals)


#Calculate a 95% CI for the expected value of y at 22.4 and 30.4
predict(object = SLR.model, newdata = data.frame(c(22.4, 30.4)), interval = "confidence", level = 0.95)


#Calculate a 95% prediction interval of y at xxx and xxxx
predict(object = SLR.model, newdata = data.frame(c(22.4, 30.4)), interval = "prediction", level = 0.95)



















##############################Old code, which can be discarded################################################################
##############################################################################################
##### Read in the P content in apple tree leaves data 
setwd("/Users/adminuser/Desktop/Work/CPSC_Courses/CPSC440/Fall_2015/Lecture_Notes/")

apple.data <- read.table("Chapter_8_Apple_Data.txt", head = TRUE)


########Fit a one-way ANOVA model Y_ij = mu + tau_i + epsilon_ij
ANOVA.model <- lm(apple.data$value ~ as.factor(apple.data$trt))

#Look at the ANOVA table
anova(ANOVA.model)

#Obtain the predicted values and residuals
residuals <- ANOVA.model$residuals
predicted.values <- ANOVA.model$fitted.values

####Assess the assumption of constant variance

#Plot residuals against predicted values
plot(residuals~predicted.values)
abline(h = 0, lty = 2)

#Conduct the Brown-Forsyth test to test Ho: Variance within the three groups are equal
install.packages("lawstat") #Note you only need to run this command the first time you run this

library(lawstat)

levene.test(apple.data$value, apple.data$trt)

#######Assess the assumption of normally distributed error terms

#Make a QQ-plot of the residuals
qqnorm(residuals)


#Make a histogram of the residuals
hist(residuals)

#Conduct a Shapiro-Wilk test for normality of the residuals, which tests Ho: Residuals are normally distributed
shapiro.test(residuals)


##############################################################################################
##############################################################################################
##############################################################################################
#######Conduct the Kruskal-Wallis Test

#Read in the cleric data
cleric.data <- read.table("Chapter_8_Cleric_Data.txt", head = TRUE)

##Conduct the Kruskal-Wallis test
kruskal.test(testscore ~ denomination, data = cleric.data)







