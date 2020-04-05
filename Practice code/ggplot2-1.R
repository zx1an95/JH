
install.packages("ggplot2")
install.packages("RColorBrewer")

library(ggplot2)
library(RColorBrewer)

salaries <- read.csv( "Salaries.csv")

#Explore the dataset
head(salaries)
str(salaries)
summary(salaries)

ranks <- table(salaries$rank)

pie(ranks)
mycol <- c("darkorange","forestgreen","steelblue4")
pie(ranks,
    labels = names(ranks),
    col = mycol,
    main = "Teaching Faculty")

percent = round(100 * ranks / sum(ranks), 1)
rank.names <- c("Associate Professor", "Assistant Professor", "Professor")
labels <- paste( rank.names, " ",percent, "%",     sep="")
pie( ranks,
     labels = labels,
     col = mycol ,
     init.angle = 90,
     main = "Teaching Faculty")

text(x=0.25,y=-1.0, 
     labels = "Source: University Administration", 
     pos = 4, 
     cex=0.75)
dev.off()

#Check opened devices
dev.list()
pie( ranks, labels = labels, col = mycol , init.angle = 90, main = "Teaching Faculty")
dev.list()
dev.cur()
dev.copy(png,"piechart.png")
dev.off() 
png("piechart.png")
pie( ranks, labels = labels, col = mycol , init.angle = 90, main = "Teaching Faculty")
dev.off() 

b<-barplot(ranks,
           col = mycol,
           main = "Teaching Faculty",
           sub = "Source: University Administration",
           names.arg = rank.names, 
           ylim = c(0, max(ranks) * 1.1),
           axes = FALSE)

text(b,ranks,
     labels = ranks, 
     pos = 3)

boxplot(salaries$salary)
boxplot( formula = salary ~ discipline,
         data = salaries,
         names = c("A","B"),
         col = c("gold","forestgreen"),
         horizontal = TRUE,
         main = "Salary range for Disciplines A and B")

levels(salaries$rank)
col3 <- c("steelblue4","firebrick","forestgreen")
cols <- col3 [ as.numeric(salaries$rank) ]
plot(x = salaries$yrs.service,
     y = salaries$salary, 
     pch = 19, 
     main = "Analysis of dependency of 2 variables",
     xlab = "Years of Service",
     ylab = "Salary", 
     col = cols,
     cex=1.5)

lm.fit <- lm(salary ~ yrs.service, data = salaries)
abline(lm.fit, lty = 4, lwd = 2, col="darkgray")
legend("topleft",levels(salaries$rank),fill=col3 )

eq <- substitute(italic(salary) == a + b %.% italic(yrs.service)*","~~italic(r)^2~"="~r2, 
                 list(a = format(coef(lm.fit)[1], digits = 2), 
                      b = format(coef(lm.fit)[2], digits = 2), 
                      r2 = format(summary(lm.fit)$r.squared, digits = 3)))
text(15,170000,labels=as.expression(eq), pos=4, cex=.8 )
hist(salaries$salary)
hist(salaries$salary, 
     main="Histogram of Salaries of professors in Disciplines A and B",
     col="steelblue3",
     xlab="Salary, dollars")

hist(salaries$salary, 
     main="Histogram of Salaries of professors in Disciplines A and B",
     col="steelblue3",
     xlab="Salary, dollars", 
     prob=TRUE)

par( mfrow = c(1,2) )  # 1 row and 2 columns

pie( ranks, labels = labels, col = mycol , init.angle = 90, main = "Teaching Faculty")

b<-barplot(table(salaries$rank),
           col=c("steelblue4","firebrick","forestgreen"),
           main = "Professor Ranks",
           names.arg = c("Associate Professor", "Assistant Professor", "Professor"),
           ylim = c(0, max(table(salaries$rank)) * 1.25), 
           axes=FALSE)
text(b,table(salaries$rank),table(salaries$rank), pos=3)
par ( mfrow = c(1,1) )
par(par.save)

library(ggplot2)
ggplot(salaries, aes(x=factor(1), fill=rank)) + 
  geom_bar(width=1)

ggplot(salaries, aes(x=factor(1), fill=rank)) + 
  geom_bar(width=1) + 
  coord_polar("y") + 
  ggtitle("Teaching Faculty")

ggplot(salaries, aes(x=factor(1), fill=rank)) + 
  geom_bar(width=1) + 
  coord_polar("y") +
  ggtitle("Teaching Faculty") + 
  theme_minimal() + 
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    plot.title=element_text(size=14,face="bold")
  )

ggplot(salaries, aes( yrs.service, salary)) + geom_point() 

ggplot(salaries, aes( yrs.service, salary)) + 
  geom_point() +    # type of geom objects
  theme_bw() +      # theme
  ggtitle(" Salary analysis")   # plot title

ggplot(salaries, aes( yrs.service, salary)) + 
  geom_point( aes(color = factor(rank)  ), size=5) +    # type of geom objects
  ggtitle(" Salary analysis") +  # plot title
  theme_bw( ) +      # theme
  theme ( plot.title=element_text(size=14,face="bold" ) )

ggplot(salaries, aes( yrs.service, salary)) + 
  geom_point( aes(color = factor(rank)  ), size=5) +    # type of geom objects
  ggtitle(" Salary analysis") +  # plot title
  geom_smooth() + 
  theme_bw( ) +      # theme
  theme ( plot.title=element_text(size=14,face="bold" ) )

ggplot(salaries, aes( yrs.service, salary)) + 
  geom_point( aes(color = factor(rank)  ), size=5) +    # type of geom objects
  ggtitle(" Salary analysis") +  # plot title
  geom_smooth( method='lm' , col = "black") + 
  theme_bw( ) +      # theme
  theme ( plot.title=element_text(size=14,face="bold" ) )

ggplot(salaries, aes( x=discipline, y=salary ) ) + geom_boxplot() 

ggplot(salaries, aes( x=discipline, y=salary ) ) + 
  geom_boxplot() +
  geom_jitter(width = 0.2)

ggplot(salaries, aes( x=rank, y=salary) ) + 
  geom_boxplot( varwidth = TRUE ) + #width is proportional to the square roots of the number of observations
  geom_jitter(width = 0.2) + 
  ggtitle(" Salary analysis") + 
  coord_flip() + 
  theme ( plot.title=element_text(size=14,face="bold" ) )

ggplot(salaries, aes(x=salary)) + geom_histogram()
ggplot(salaries, aes(x=salary)) + geom_histogram(binwidth=10000)

ggplot(salaries, aes(x=salary)) + 
  geom_histogram(binwidth=10000, 
                 color="black", fill="white")

ggplot(salaries, aes(x=salary)) + 
  geom_histogram(binwidth=10000, 
                 color="black", fill="white", aes(y=..density..)) + 
  geom_density(alpha=.2, fill="darkorange")

