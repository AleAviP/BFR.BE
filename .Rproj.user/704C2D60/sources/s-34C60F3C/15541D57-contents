n<-100
p<-1000
q<-10

grid<-seq(-1,1,length.out=p)
M<-matrix(grid,ncol=q,nrow=p)
rate<-trunc(p/(q*2))
for(i in 2:q){
  M[,i]<-grid[c((i*rate):p,1:(i*rate-1))]
}

E <- rmvnorm(n,numeric(p),diag(p))
Z <- rmvnorm(n,numeric(q),diag(q))
X <- Z%*%t(M) + E

result=BFR.BE.EM.CV(x=X, q = 10,varianceBE=FALSE)
