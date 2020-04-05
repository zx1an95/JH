x1<-c(1,2,3)
dim(x1)<-c(1,3)
x<-t(x1)
y1<-c(4,5,6)
dim(y1)<-c(1,3)
y<-t(y1)
e<-array(1,dim=c(3,1))
z<-2*x+y+e
print(z)
crossprod(x,y)
tcrossprod(x,y)

A<-matrix(1:20,nrow=4,ncol=5,byrow=FALSE)
B<-matrix(1:20,nrow=4,ncol=5,byrow=TRUE)
C=A+B
D=A*B
E1=array(D,dim=c(20,1))
tcrossprod(E1)
F=matrix(A[1:3,1:3],ncol=3)
G=B[,c(1,2,4,5)]
G<-B[,-3]
X=c(rep(1,5),rep(2,3),rep(3,4),rep(4,2))
X<-as.vector(X)
n<-5;H<-array(0,dim=c(5,5))
for(i in 1:n){
  for(j in 1:n){
    H[i,j]<-1/(i+j-1)
  }
}
det(H)
solve(H)
eigen(H)
svd(H)