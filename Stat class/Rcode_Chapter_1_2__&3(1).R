#####Chicken example

#Read in the data
data <- c(3.7, 4.2, 4.4, 4.4, 4.3, 4.2, 4.4, 4.8, 4.9, 4.4,
          4.2, 3.8, 4.2, 4.4, 4.6, 3.9, 4.3, 4.5, 4.8, 3.9,
          4.7, 4.2, 4.2, 4.8, 4.5, 3.6, 4.1, 4.3, 3.9, 4.2,
          4.0, 4.2, 4.0, 4.5, 4.4, 4.1, 4.0, 4.0, 3.8, 4.6,
          4.9, 3.8, 4.3, 4.3, 3.9, 3.8, 4.7, 3.9, 4.0, 4.2,
          4.3, 4.7, 4.1, 4.0, 4.6, 4.4, 4.6, 4.4, 4.9, 4.4,
          4.0, 3.9, 4.5, 4.3, 3.8, 4.1, 4.3, 4.2, 4.5, 4.4,
          4.2, 4.7, 3.8, 4.5, 4.0, 4.2, 4.1, 4.0, 4.7, 4.1,
          4.7, 4.1, 4.8, 4.1, 4.3, 4.7, 4.2, 4.1, 4.4, 4.8,
          4.1, 4.9, 4.3, 4.4, 4.4, 4.3, 4.6, 4.5, 4.6, 4.0)

mean(data)
sd(data)
length(data)

hist(data, xlim = c(3.5,5.0), breaks = 20, col = "blue") 

hist(data, xlim = c(3.5,5.0), breaks = 10, col = "blue") 

hist(data, xlim = c(3.5,5.0), breaks = 5, col = "blue") 

hist(data, xlim = c(3.5,5.0), breaks = 2, col = "blue") 


####################Flight delay example

#Set the working directory:
setwd("/Users/adminuser/Desktop/Work/CPSC_Courses/CPSC440/Fall_2016/Lecture_Notes/")


flight.delay.data <- read.csv("Lecture_1_Flight_Delays_1288.csv", head = TRUE)

means <- by(flight.delay.data$Delays, flight.delay.data$Date, mean) 

pdf("Flight_Delay_Boxplot.pdf")
  boxplot(flight.delay.data$Delays ~ flight.delay.data$Date, col = "Red", 
          xlab = "Date", ylab = "Delay (in minutes)")
  points(c(1:length(means)), means, pch = 3, cex = 0.75) #Adds the means, which are indicated by the "+" symbol
dev.off()

