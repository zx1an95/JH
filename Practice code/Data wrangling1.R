
library(readxl)
library(dplyr)


visit <- read_excel("MedData.xlsx", sheet = "Visits")

visit <- read_excel("MedData.xlsx", 
                    sheet = "Visits", 
                    skip = 3,  # how many lines to skip at the top of the sheet
                    n_max = 13)  # how many observations to read
str(visit)
glimpse(visit)

visit <- read_excel("MedData.xlsx", 
                    sheet = "Visits", 
                    skip = 3,
                    n_max = 13,
                    na=c("","NA","##N/A"))
glimpse(visit)

visit %>% glimpse()
visit %>% summary()

sort(unique(visit$`Patient ID`), decreasing=TRUE)

visit$`Patient ID` %>% 
  unique()  %>%            # find unique values
  sort(decreasing=TRUE)    # sort  in descending order



grepl ("pain", visit.clean$Symptoms,  ignore.case=TRUE)
grepl ("fever|pain", visit.clean$Symptoms,  ignore.case=TRUE)

grep ("fever", visit.clean$Symptoms,  ignore.case=TRUE, value=FALSE)
grep ("fever", visit.clean$Symptoms,  ignore.case=TRUE, value=TRUE)
sub (",", ";", visit.clean$Symptoms)
gsub (",", ";", visit.clean$Symptoms)

visit.clean %>% filter( !is.na(DBP) )
visit.clean %>% select(id, Temperature:pulse)
visit.clean %>% select( id, ends_with("BP") )
visit.clean %>% arrange(id, admission)
visit.clean %>% arrange(id, desc(admission))
pinfo <- read_excel("MedData.xlsx", sheet = "Patient Data")
pinfo %>%head()
pinfo %>%glimpse()
pinfo %>%summary()
full_join(visit.clean, pinfo, by = c("id"="ID"))
result <- visit.clean  %>% left_join(pinfo, by = c("id"="ID"))

library(tidyr)
city.temps <- data.frame( time = as.Date('2018-09-03') + 0:4,
                          Boston = rnorm( 5, 75, 3),
                          Denver = rnorm( 5, 85, 5),
                          SanFrancisco = rnorm(5, 80, 5),
                          Austin = rnorm( 5, 90, 5) )
city.temps2 <- gather( city.temps, 
                       key = "City",
                       value = "Temperature",
                       -time,  
                       factor_key = TRUE) 

city.temps2 <- gather( city.temps, 
                       key = "City",
                       value = "Temperature",
                       Boston: Austin,  
                       factor_key = TRUE) 

city.temps2 <- gather( city.temps, 
                       key = "City",
                       value = "Temperature",
                       c("Boston","Denver","SanFrancisco","Austin"),  
                       factor_key = TRUE) 


city.temps2
glimpse(city.temps2)
city.temps3 <- spread( city.temps2, City, Temperature)
city.temps3