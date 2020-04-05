
# Use the non-central F distribution to calculate power

#Input values
diff <- 15
sigma <- 7.5
alpha <- 0.05
number.of.treatments <- 4
rep.number.vector <- c(6,8,9,10)

power.vector <- NULL
beta.vector <- NULL

for(i in rep.number.vector){
  num.df <- number.of.treatments - 1
  den.df <- number.of.treatments*(i-1)
  critical.F.value <- qf(p = 0.05, df1 = num.df, df2 = den.df, lower.tail = FALSE) 
  lambda <- (i*(diff^2))/(2*(sigma^2))
  power <- pf(q = critical.F.value, df1 = num.df, df2 = den.df, ncp = lambda, lower.tail = FALSE)
  beta <- 1-power
  power.vector <- c(power.vector, power)
  beta.vector <- c(beta.vector, beta)
}#end for(i in rep.number.vector)

final.table <- cbind(rep.number.vector, power.vector, beta.vector)
colnames(final.table) <- c("Reps", "Power", "Beta")
