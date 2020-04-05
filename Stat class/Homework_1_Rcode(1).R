#####Simulate a predetermined number of flips of two coins.
#### To adjust the number of flips of two coins, please just
### change the value of “number of flips”

#Instructions:
#1.) Highlight and run lines 17-60
#2.) To answer the homework question, run line 10


multiple.flips.of.coin(1000)






#Don't change the code below
multiple.flips.of.coin <- function(number.of.flips = NA){

  nTT <- 0
  nTH <- 0
  nHT <- 0
  nHH <- 0
  
  for(i in 1:number.of.flips){
    #Strategy: flip two coins, recorded as a random number between 0 and 1.
    # If the random number is > 0.25, call it "heads"; if it is < 0.25, call it "tails"
    tosscoin1 <- runif(n = 1, min = 0, max = 1)
    tosscoin2 <- runif(n = 1, min = 0, max = 1)
    
    #Note that "&" means "and"
    if((tosscoin1 >  0.90) & (tosscoin2 >  0.90)) nTT <- nTT + 1
    if((tosscoin1 >  0.90) & (tosscoin2 <=  0.90)) nTH <- nTH + 1
    if((tosscoin1 <=  0.90) & (tosscoin2 >  0.90)) nHT <- nHT + 1
    if((tosscoin1 <=  0.90) & (tosscoin2 <=  0.90)) nHH <- nHH + 1
    
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
  
  
  #Explicitly print out the proporiton of heads that were obtained
  print(paste("The proportion of ", number.of.flips, " flips in which exactly zero heads was obtained is: ",
              Rel.freq.TT, sep = ""))
  
  print(paste("The proportion of ", number.of.flips, " flips in which exactly one head was obtained is: ",
              (Rel.freq.TH + Rel.freq.HT), sep = ""))
  
  print(paste("The proportion of ", number.of.flips, " flips in which exactly two heads was obtained is: ",
            Rel.freq.HH, sep = ""))
}#end flip.the.coin


