Lab9 <- read.csv("Lab9.csv",header = TRUE)
my_lm<-lm(Lab9$Prod~Lab9$Destination)

shapiro.test(residuals(my_lm))

install.packages("car")
library(car)
leveneTest(Lab9$Prod,Lab9$Destination)

install.packages("MASS")
library(MASS)
boxcox(my_lm)

new.variable<-(Lab9$Prod)^0.5

Lab9 <- cbind(Lab9,new.variable)

transformed_lm<-lm(Lab9$new.variable~Lab9$Destination)
anova(transformed_lm)

shapiro.test(residuals(transformed_lm))

library(car)
leveneTest(Lab9$new.variable,Lab9$Destination)

kruskal.test(Prod~Destination, data = Lab9)
