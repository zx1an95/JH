the.mean<-NULL
a<-NULL
b<-NULL
dd<-NULL
for(i in unique(hw13q2$Half_Marathon_Time_Indianapolis)){
  + if(sum(hw13q2$Half_Marathon_Time_Indianapolis==i)>1){
    + these.data<-hw13q2[which(hw13q2$Half_Marathon_Time_Indianapolis==i),]
    + this.mean <- mean(these.data$Half_Marathon_Time_Champaign_Urbana)
    + the.mean<-c(the.mean,this.mean)
      a<-c(a,i)
      b<-c(b,sum(hw13q2$Half_Marathon_Time_Indianapolis==i))
      ss<-sum((these.data$Half_Marathon_Time_Champaign_Urbana-this.mean)^2)
      dd<-c(dd,ss)
    + }
  + }


