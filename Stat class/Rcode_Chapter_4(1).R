#####Simulate 10,000,000 flips of two coins

#Initialze the variables
number.of.flips <- 10000000
nTT <- 0
nTH <- 0
nHT <- 0
nHH <- 0

for(i in 1:number.of.flips){
  #Strategy: flip two coins, recorded as a random number between 0 and 1.
  # If the random number is > 0.5, call it "heads"; if it is < 0.5, call it "tails"
  tosscoin1 <- runif(n = 1, min = 0, max = 1)
  tosscoin2 <- runif(n = 1, min = 0, max = 1)
  
  #Note that "&" means "and"
  if((tosscoin1 < 0.5) & (tosscoin2 < 0.5)) nTT <- nTT + 1
  if((tosscoin1 < 0.5) & (tosscoin2 > 0.5)) nTH <- nTH + 1
  if((tosscoin1 > 0.5) & (tosscoin2 < 0.5)) nHT <- nHT + 1
  if((tosscoin1 > 0.5) & (tosscoin2 > 0.5)) nHH <- nHH + 1
  
}# end for(i in 1:number.of.flips)

#Calculate the relative frequencies
Rel.freq.TT <- nTT/number.of.flips
Rel.freq.TH <- nTH/number.of.flips
Rel.freq.HT <- nHT/number.of.flips
Rel.freq.HH <- nHH/number.of.flips

#View the results
Rel.freq.TT
Rel.freq.TH
Rel.freq.HT 
Rel.freq.HH



################################################################################
# Using functions in R to calcualte probabilities of a binomial random variable

#Let Y = number of seeds that germinate, 
# Supppose that Y follows a binomial distribution with n = 20 and pi = 0.85


#Calculate the probability that Y = 18
dbinom(18, size = 20, p = 0.85)

#Calculate the probability that Y >= 18
dbinom(18, size = 20, p = 0.85)+dbinom(19, size = 20, p = 0.85)+dbinom(20, size = 20, p = 0.85)

sum(dbinom(c(18,19,20), size = 20, p = 0.85))

#Calculate the probability that Y < 18
Prob.Y.geq.18 <- dbinom(18, size = 20, p = 0.85)+dbinom(19, size = 20, p = 0.85)+
                 dbinom(20, size = 20, p = 0.85)
1-Prob.Y.geq.18



################################################################################
# Using functions in R to calcualte probabilities for normal random variables

#Let X be a normal random variable with mean 20 and sigma 2. Then P{Z <= 1.5}
# is calculated as follows:
pnorm(23, mean = 20, sd = 2, lower.tail = TRUE)

#Let Z be a normal random variable with mean 0 and sigma 1. Then P{Z <= 1.5}
# is calculated as follows:
pnorm(1.5, mean = 0, sd = 1, lower.tail = TRUE)

#Let Z be a normal random variable with mean 0 and sigma 1. What is the value of Z
# in which 80% of all values are less that it?
qnorm(0.8, mean = 0, sd = 1)


###############################################################################
# Selecting a random number in R

#Select 5 random numbers from 1 to 40
sample(1:40, 5)

################################################################################
# Let Y be a binomial random variable with n = 1000 and pi = 0.5. 
# Then P{Y <= 460 } is calculated as follows:

sum(dbinom(0:460, size = 1000, p = 0.50))





