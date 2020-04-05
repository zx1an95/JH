library(MASS)

plot( Petal.Width ~ Petal.Length, data=iris )
sunflowerplot(Petal.Width ~ Petal.Length, data=iris)
pairs( iris[, 1:4] )
library(corrplot)
numvals <- iris[, 1:4]
corrm <- cor(numvals)
corrplot(corrm, method="ellipse")
plot( Bwt ~ Hwt, data = cats)

ln.res <- lm(Bwt ~ Hwt, data = cats)

colors <- c("brown", "steelblue")
plot( Bwt ~ Hwt, data = cats, 
      main = " Relationship between body weight and heart weight in cats",
      xlab = "Heart weight (g)",
      ylab = "Body weight (kg)",
      cex = 1.1,
      pch = 16,
      col=colors[ cats$Sex ] )

legend("topleft",col=colors[ unique(cats$Sex) ], pch=16, legend=levels(cats$Sex))

abline(ln.res, lwd=2, lty=2)

dev.off()
dev.list()

colors <- c("brown", "steelblue")
plot( Bwt ~ Hwt, data = cats, 
      main = " Relationship between body weight and heart weight in cats",
      xlab = "Heart weight (g)",
      ylab = "Body weight (kg)",
      cex = 1.1,
      pch = 16,
      col=colors[ cats$Sex ] )

legend("topleft",col=colors[ unique(cats$Sex) ], pch=16, legend=levels(cats$Sex))

abline(ln.res, lwd=2, lty=2)

dev.list()
dev.cur()


# Copy the existing graph to the file
dev.copy(png,"cats.png")

#Close the png file
dev.off() 

dev.off()

dev.list()



png("cats.png")

colors <- c("brown", "steelblue")
plot( Bwt ~ Hwt, data = cats, 
      main = " Relationship between body weight and heart weight in cats",
      xlab = "Heart weight (g)",
      ylab = "Body weight (kg)",
      cex = 1.1,
      pch = 16,
      col=colors[ cats$Sex ] )

legend("topleft",col=colors[ unique(cats$Sex) ], pch=16, legend=c("Female","Male"))

abline(ln.res, lwd=2, lty=2)

dev.off() 

par()


par.save <- par( no.readonly=TRUE)


par( bg = "lightgray", no.readonly=TRUE)

pie( table(cats$Sex), 
     labels = c("Female","Male"), 
     col = c("brown", "steelblue") , 
     init.angle = 90, 
     main = "Cats")

#restore background to original value 
par(par.save)  


x <- 1:10
y1 <- log(x)
y2 <- sqrt(x)


plot(x,y1)
lines(x,y2)


plot(x, y1,                                          
     type = "b",                                     
     ylim = c(min(c(y1,y2)), max(c(y1,y2))),         
     main = "Logarithm and Square root functions",   
     ylab = "y",     
     col = "darkgreen",                             
     lwd = 2,                                       
     cex.main = 2)                                   

lines(x,y2,                       
      col = "steelblue",           
      lwd =2,                      
      type="b",                    
      pch=19,                      
      cex.lab=1.5)               

#Add grid to the plot
abline(h = seq(0,3,by=0.5), v = c(0:10), col='lightgray',lty="dotted")


legend("topleft",                             
       col=c("darkgreen","steelblue"),       
       lty=1,                                 
       pch=c(15,19),                          
       cex=1.5,                               
       legend=c("logarithm","square root"))   # text


text(x = 6, y = sqrt(6), pos = 2, labels=expression(paste("y = ",sqrt(x))), cex=2, offset=0.75)
text(x = 6, y = log(6),  pos = 2, labels=expression(paste("y = ",ln(x))),   cex=2, offset=0.75)



cats.sex <- table(cats$Sex)
barplot(cats.sex)

colors <- c("brown", "steelblue")
b<-barplot(cats.sex,
           col = colors,
           main = "Cats",
           sub = "R. A. Fisher (1947)",
           names.arg = c("Female","Male"), 
           ylim = c(0, max(cats.sex) * 1.1),
           axes = FALSE)

text(b,cats.sex,
     labels = cats.sex, 
     pos = 3)
title("Anatomical Data from Domestic Cats", line=0, cex.main=.8)

data <- matrix( c(45,78,90, 30, 50, 60), nrow=2, byrow=TRUE)
data

# By default barplot will show stacked bars
barplot(data)

colors <- c("brown", "steelblue")
b <- barplot(data, 
             main="Hours worked by 2 teams during 2015-2017 ", 
             col=colors, 
             ylim = c(0, max(colSums(data)) * 1.2))
text(b, c(45,78,90, 45+30, 78+50, 90+60), c(45,78,90, 30, 50, 60), pos=1)

# Side by side plot
colors <- c("brown", "steelblue")
b <- barplot(data, 
             main="Hours worked by 2 teams during 2015-2017 ", 
             col=colors, 
             ylim = c(0, max(data) * 1.2),
             beside=T)
b
text(b, data, data, pos=3)

legend("topleft",col=colors, pch=15, legend=c("team1", "team2"))

barplot(c(1:10), col=c(1:10))

palette()

pie(rep(1,length(palette())),
    col=palette())

num.col=15
barplot(c(1:num.col), col=rainbow(num.col))
barplot(c(1:num.col), col=heat.colors(num.col))
barplot(c(1:num.col), col=terrain.colors(num.col))
barplot(c(1:num.col), col=topo.colors(num.col))
barplot(c(1:num.col), col=cm.colors(num.col))
barplot(c(1:num.col), col=gray.colors(num.col))

#install.packages("RColorBrewer")
library(RColorBrewer)
display.brewer.all()
barplot(c(1:8), col=brewer.pal(8,"Dark2") )
barplot(c(1:9), col=brewer.pal(9,"Set1") )

colors() 
barplot(c(21,19,16,14,12,10,8,5,4,3,2,1), col=c("steelblue4",
                                                "darkred",
                                                "dodgerblue",
                                                "forestgreen",
                                                "gold",
                                                "steelblue",
                                                "firebrick",
                                                "darkorange",
                                                "darkgrey",
                                                "tomato",
                                                "olivedrab",
                                                "dimgray"))

barplot(c(21,19,17,14,12,8,5,4,3,2,1), col=c("springgreen",
                                             "skyblue",
                                             "violetred",
                                             "slateblue",
                                             "sienna",
                                             "seagreen",
                                             "sandybrown",
                                             "salmon",
                                             "saddlebrown",
                                             "maroon",
                                             "limegreen"))



col1 <- rgb(0.1, 0.5, 0.1)
col2 <- rgb(0.9, 0.7, 0.0)
barplot(c(7,11), col=c(col1, col2))

salaries <- read.csv( "Salaries.csv")
boxplot(salaries$salary)

boxplot( formula = salary ~ discipline,
         data = salaries,
         names = c("A","B"),
         col = c("gold","forestgreen"),
         horizontal = TRUE,
         main = "Salary range for Disciplines A and B")

plot(c(-1, 26), 0:1, type = "n", axes = FALSE, xlab="", ylab="")
text(0:25, 0.6, 0:25, cex=0.75)
points(0:25, rep(0.5, 26), pch = 0:25, bg = "grey")


plot(x = salaries$yrs.service,
     y = salaries$salary, 
     pch = 19)

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

rug(salaries$salary)
lines(density(salaries$salary), col="darkblue", lwd=2)
pairs( iris[, 1:4] )
pairs( ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width, data=iris)

pairs( ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width, 
       data = iris,
       lower.panel = panel.smooth, 
       upper.panel = NULL) 
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

pairs( ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width, 
       data = iris,
       lower.panel = panel.smooth, 
       upper.panel = panel.cor) 


panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}
pairs(~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width, 
      data = iris,
      panel = panel.smooth,
      cex = 1.5, pch = 16, bg = "light blue",
      diag.panel = panel.hist, cex.labels = 2, font.labels = 2)




par( mfrow = c(1,2) ) 


ranks <- table(salaries$rank)
labels <-c ("Associate Professor", "Assistant Professor", "Professor")
pie( ranks, labels = labels, col = 2:4 , init.angle = 90, main = "Teaching Faculty")

b<-barplot(table(salaries$rank),
           col=c("steelblue4","firebrick","forestgreen"),
           main = "Professor Ranks",
           names.arg = c("Associate Professor", "Assistant Professor", "Professor"),
           ylim = c(0, max(table(salaries$rank)) * 1.25), 
           axes=FALSE)
text(b,table(salaries$rank),table(salaries$rank), pos=3)


par ( mfrow = c(1,1) )
par(par.save)
