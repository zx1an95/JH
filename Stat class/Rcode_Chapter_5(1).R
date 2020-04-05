
setwd(class.dir)
#Calculate a P-value of a t-distribution with a given degrees of freedom
pt(1.782, 12, lower.tail = FALSE)

#Calculate a critical value of a t-distribution with given degrees of freedom
qt(0.05, 12, lower.tail = FALSE)

x <- 33:41
beta <- c(.9500, .8266, .5935, .3200, .1206,
          .0301, .0049, .0005, .0001)
power <- 1-beta

#Graph the type II error rate at the value of the alternative hypothesis given by x
pdf("Type_II_and_Power_Curves.pdf", width = 6)
 
  #Set up the output graph so that there are two graphs on one sheet and the margins looks ok
  par(mfrow = c(2,1), mar = c(5,5,5,5))
  
  #Plot the Type II error rates
  plot(beta~x, main = "Type II Error Rate")
  lines(spline(x, beta, n = 201), col = 2)
  
  #Plot the power levels
  plot(power~x, main = "Power")
  lines(spline(x, power, n = 201), col = 2)  

dev.off()
Ecoli <- c(0.593, 0.142, 0.329, 0.691, 0.231, 0.793, 0.519, 0.392, 0.418)

shapiro.test(Ecoli)