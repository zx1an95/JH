##############################################################################################
# Read in the required packages
install.packages("SimComp") #Run this only if you do not have "SimComp" already installed on your computer
library(SimComp)

##### Read in the weed data
setwd("/Users/adminuser/Desktop/Work/CPSC_Courses/CPSC440/Fall_2016/Lecture_Notes/")

weed.data <- read.csv("Chapter_9_Weed_Data.csv", head = TRUE)
attach(weed.data)

#Fit an ANOVA model Y_ij = mu + method_i + epsilon_ij
# and look at the ANOVA table
ANOVA.model <- lm(y ~ method)
anova(ANOVA.model)


#####Make contrasts
contrast.matrix <- rbind(c(-1,-1,-1,-1,4),c(0.5,0.5,-0.5,-0.5,0),c(1,-1,0,0,0),c(0,0,1,-1,0) )
rownames(contrast.matrix) <- c("control vs agents", "bios vs chems", "bio1 vs bio2", 
                               "chem1 vs chem2")

SimTestDiff(data =weed.data, grp = "method", resp = "y", ContrastMat = contrast.matrix, covar.equal = TRUE)

#Conduct the Tukey and Dunnett Adjustements for the mutliple testing problem.
SimTestDiff(data =weed.data, grp = "method", resp = "y", type = "Dunnett", base = 5, covar.equal = TRUE)

SimTestDiff(data =weed.data, grp = "method", resp = "y", type = "Tukey", covar.equal = TRUE)

