HW11 <- read.csv("HW11.csv", header = TRUE)
my_lm <- lm(Time ~ Timeofday*Route, data=HW11)
shapiro.test(residuals(my_lm))
library(car)
leveneTest(my_lm)
anova(my_lm)
#View interaction
#For A by B
interaction.plot(HW11$Timeofday,HW11$Route,HW11$Time)
#for B by A
interaction.plot(HW11$Route,HW11$Timeofday,HW11$Time)