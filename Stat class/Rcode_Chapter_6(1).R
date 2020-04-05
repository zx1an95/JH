############################################
#####Conduct a two-sample t-test

Animal <- 1:12
Drug.treatment <- c(18, 43, 28, 50, 16, 32, 13, 35, 38, 33, 6, 7)
Untreated <- c(40, 54, 26, 63, 21, 37, 39, 23, 48, 58, 28, 39)


####One-side lower test
t.test(Drug.treatment, Untreated, alternative="less", var.equal=TRUE) #pooled t-test
t.test(Drug.treatment, Untreated, alternative="less", var.equal=FALSE) #approximate t-test for unequal variances

####One-side upper test
t.test(Drug.treatment, Untreated, alternative="greater", var.equal=TRUE) #pooled t-test
t.test(Drug.treatment, Untreated, alternative="greater", var.equal=FALSE) #approximate t-test for unequal variances

####Two-sided test
t.test(Drug.treatment, Untreated, alternative="two.sided", var.equal=TRUE) #pooled t-test
t.test(Drug.treatment, Untreated, alternative="two.sided", var.equal=FALSE) #approximate t-test for unequal variances

############################################################
#####Conduct a Wilcoxson (i.e., Mann-Whitney) Rank Sum test

placebo <- c(0.90, 0.37, 1.63, 0.83, 0.95, 0.78, 0.86, 0.61, 0.38, 1.97)
alcohol <- c(1.46, 1.45, 1.76, 1.44, 1.11, 3.07, 0.98, 1.27, 2.56, 1.32)

wilcox.test(placebo, alcohol, alternative = "less")

############################################
#####Conduct a paired t-test

garage1 <- c(17.6, 20.2, 19.5, 11.3, 13, 16.3, 15.3, 16.2, 12.2, 14.8, 21.3, 22.1, 16.9, 17.6, 18.4)
garage2 <- c(17.3, 19.1, 18.4, 11.5, 12.7, 15.8, 14.9, 15.3, 12, 14.2, 21, 21, 16.1, 16.7, 17.5)
t.test(garage1, garage2, paired = TRUE, alternative = "two.sided")


############################################################
#####Conduct a Wilcoxson Signed-Rank test
brand.A <- c(211.4, 204.4, 202, 201.9, 202.4, 202, 202.4, 207.1, 203.6, 216, 208.9, 208.7, 213.8, 201.6, 201.8, 200.3, 201.8, 201.5, 212.1, 203.4)
brand.B <- c(186.3, 205.7, 184.4, 203.6, 180.4, 202, 181.5, 186.6, 205.7, 189.1, 183.6, 188.7, 188.6, 204.2, 181.6, 208.7, 181.5, 208.7, 186.8, 182.9)
wilcox.test(brand.A, brand.B, paired = TRUE, alternative = "greater")
