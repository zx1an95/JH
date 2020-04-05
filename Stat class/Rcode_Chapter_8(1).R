
setwd("/Users/adminuser/Desktop/Work/CPSC_Courses/CPSC440/Fall_2018/Lecture_Notes/")

apple.data <- read.csv("Chapter_8_Apple_Data.csv", head = TRUE)


########Fit a one-way ANOVA model Y_ij = mu + tau_i + epsilon_ij
ANOVA.model <- lm(apple.data$value ~ as.factor(apple.data$trt))

residuals <- ANOVA.model$residuals
predicted.values <- ANOVA.model$fitted.values

plot(residuals~predicted.values)
abline(h = 0, lty = 2)

library(car)
leveneTest(value~as.factor(trt), data = apple.data)
qqnorm(residuals)


#Make a histogram of the residuals
hist(residuals)
shapiro.test(residuals)
psychology.data <- read.csv("Chapter_8_Psychology_Data.csv", head = TRUE)
kruskal.test(testscore ~ university, data = psychology.data)







