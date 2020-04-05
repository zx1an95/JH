library(dplyr)
library(data.table)
library(stringr)
library(microbenchmark)


dt <- data.table ( ID = c(rep(1034,3), rep(5621,2),rep(3468,4)), 
                   admission = as.POSIXct(c("2018-03-14","2018-05-28", "2018-11-03", "2018-04-15","2018-10-11","2018-02-09","2018-03-17", "2018-08-31","2018-10-12")),
                   discharge = as.POSIXct(c("2018-03-16","2018-05-29", "2018-11-06", "2018-04-19","2018-10-18","2018-02-15","2018-03-28", "2018-09-15","2018-10-19")),
                   age = c(rep(38,3), rep(45,2), rep(69,4)), 
                   diastBP = c(80, 75, 78, 60, 65, 90, 100, 105, 90), 
                   systBP = c(120,110, 115, 90, 110, 130, 135, 140, 130),
                   weight = c(180, 182, 180, 130, 135, 210, 205, 212, 207),
                   height = c(5.9, 5.9, 5.9, 5.1, 5.1, 5.4, 5.4, 5.4, 5.4),
                   symptoms=c("cough, headache, rash","cough,sore throat","diarrhea", "chest pain, cough","headache", "chest pain", "shortness of breath, chest pain", "nausea, vomiting, chest pain","cough"))

class(dt) 
data(iris)
iris.dt <- as.data.table(iris)
iris.dt
hhc_df <- as.data.frame(hhc)

class(iris)

setDT(iris)  # make iris to be a data.table
class(iris)

setDF(iris)  # make iris to be back a data.frame
class(iris)
dt[ systBP %between% c(100,120) ]

set.seed(12345)
dt_e6 <- data.table( V1 = sample(10000, 10e6, T), 
                     V2 = sample(letters, 10e6, T),
                     V3 = sample(c(T,F), 10e6, T))

set.seed(12345)
df_e6 <- data.frame( V1 = sample(10000, 10e6, T), 
                     V2 = sample(letters, 10e6, T),
                     V3 = sample(c(T,F), 10e6, T), 
                     stringsAsFactors=F)
indices(dt_e6)

dt[symptoms %like% "cough"]

dt[diastBP %between% c(60,70) ]

# Traditional base R approach
dt[diastBP >= 60 & diastBP <= 70 ]

iris.dt[ , 2]
iris.dt[ , c(2, 4) ]
dt[ ,  ID ]        # output is a vector 
dt[ , "ID" ]       # output is a data.table 
dt[ ,.(ID) ]       # output is a data.table

dt[ ,.(ID, symptoms) ]   # multiple columns

dt[ ,!c("symptoms","age") ]  

dt[, plot(diastBP)]

# We can define a multi-line function operating on columns of a data.table
dt[, {
  map = (systBP + 2*diastBP)/3 
  boxplot(map, horizontal = T)
  NULL    # Return value
}]

dt[ , .( mean_weight= mean(weight), Total=.N) ]
dt[ , map := systBP + 2*diastBP ]

# To create more than one column at a time:
dt[ , `:=` (X = rnorm(.N), Y = rbinom(.N, 100, 0.5)) ]
dt[ , c("X","Y") := NULL]   
dt

X <- data.table( a=1:5, b=rnorm(5))
setnames( X, 1:2, c("A", "B") )
setcolorder( X, c("B", "A") )

dt <- data.table ( ID = c(rep(1034,3), rep(5621,2),rep(3468,4)), 
                   admission = as.POSIXct(c("2018-03-14","2018-05-28", "2018-11-03", "2018-04-15","2018-10-11","2018-02-09","2018-03-17", "2018-08-31","2018-10-12")),
                   discharge = as.POSIXct(c("2018-03-16","2018-05-29", "2018-11-06", "2018-04-19","2018-10-18","2018-02-15","2018-03-28", "2018-09-15","2018-10-19")),
                   age = c(rep(38,3), rep(45,2), rep(69,4)), 
                   diastBP = c(80, 75, 78, 60, 65, 90, 100, 105, 90), 
                   systBP = c(120,110, 115, 90, 110, 130, 135, 140, 130),
                   weight = c(180, 182, 180, 130, 135, 210, 205, 212, 207),
                   height = c(5.9, 5.9, 5.9, 5.1, 5.1, 5.4, 5.4, 5.4, 5.4),
                   symptoms=c("cough, headache, rash","cough,sore throat","diarrhea", "chest pain, cough","headache", "chest pain", "shortness of breath, chest pain", "nausea, vomiting, chest pain","cough"))
dt [, .( weight = mean(weight), height = mean(height)), by = ID]

dt [, .( weight = mean(weight), height = mean(height)), keyby = ID]
dt [, .( .N ), by = .(admission < as.POSIXct("2018-07-01"), age > 50)]
uniqueN( c( 1,2,4,1,2,1) )
dt[ , uniqueN(ID), by = .(month(admission) )]

microbenchmark(
  aggregate(Grade ~ State, data=hhc_df, FUN=function(x){ c(mean=mean(x,na.rm=TRUE), sum=sum(x,na.rm=TRUE)) }),
  
  hhc_df %>% group_by(State) %>% summarise( mean=mean(Grade,na.rm=TRUE), sum=sum(Grade,na.rm=TRUE) ),
  
  hhc[ , .( mean=mean(Grade,na.rm=TRUE), sum=sum(Grade,na.rm=TRUE) ), by=State]
)

visits.dt <- data.table (id = c(11425, 11425 , 10873, 14562 ,19112, 19112, 19112, 18567, 14475, 15940, 15940, 15940, 15940),
                         admission=c("20-Feb-17","25-Mar-18","7-Aug-17","17-Dec-17","14-Mar-17","3-Jul-18","23-Jan-18","9-Nov-17","18-Aug-18","5-Feb-18","11-Mar-17","21-Oct-17","3-Nov-18"),
                         discharge=c("22-Feb-17","26-Mar-18","11-Aug-17","19-Dec-17", "19-Mar-17", "8-Jul-18","27-Jan-18","11-Nov-17","19-Aug-18","7-Feb-18","13-Mar-17","24-Oct-17", "4-Nov-18"),
                         temp=c(101.5,98.3,98.6,98.8,98.2,98.4,103.3,102.7,98.8,102.5,103.4,101.2,98.6),
                         DBP=c(55, 60,100,70,56,60,60,80,75,65,70,NA,70),
                         SBP=c(95,90,150,100,80,85,90,120, 110, 115,120,110,110),
                         heartrate=c(110,80,105,70,65,70,80,75,80,77,84,80, 75))
patient.dt <- data.table(ID = c(11425,10873,14567,19112,18567,14475,15940,15786),
                         Gender = c("Male","Female","Male","Female","Female","Male","Male","Female"))

merge( x = visits.dt, y = patient.dt, by.x = "id", by.y="ID")
merge( x = visits.dt, y = patient.dt, by.x = "id", by.y="ID", all = T)
merge( x = visits.dt, y = patient.dt, by.x = "id", by.y="ID", all.x = TRUE)
merge( x = visits.dt, y = patient.dt, by.x = "id", by.y="ID", all.y = TRUE)
visits.dt[patient.dt, on=.(id=ID)]   # all rows in Y will be preserved
patient.dt[visits.dt, on=.(ID=id)] 

visits.dt[patient.dt, on=.(id=ID), nomatch=0]
visits.dt[!patient.dt, on=.(id=ID)]

setkey(visits.dt, id)
haskey(visits.dt)
key(visits.dt)
tables()
