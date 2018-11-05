n<-200
p<-250
q<-10
Pv = 1
Pb = 2

grid<-seq(-1,1,length.out=p)
M<-matrix(grid,ncol=q,nrow=p)
rate<-trunc(p/(q*2))
for(i in 2:q){
  M[,i]<-grid[c((i*rate):p,1:(i*rate-1))]
}

Z <- rmvnorm(n,numeric(q),diag(q))
Psi <- diag(c(rep(.5,p)),p)
V <- matrix(ncol = Pv, nrow = n, runif(n*Pv,0,3))
Theta = matrix(ncol=Pv,c(rep(-2,round(p/2)),rep(2,p-round(p/2))))
batch = round(runif(n,0,1),0)
B = matrix(ncol=Pb, c(batch,ifelse(batch==1,0,1)))
Beta = matrix(ncol=Pb, c(rep(2,p),rep(0,p)))
Tau_inv = cbind(diag(Psi),diag(Psi)*1.5)
Er = matrix(ncol=p, nrow=n)
for (i in 1:n){
  if(i%%10==0){print(i)}
  Er[i,] = mvrnorm(1, rep(0,p), diag(c(Tau_inv %*% B[i,]),p), tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
}

X <- Z%*%t(M) + V%*%t(Theta) + B%*%t(Beta) + Er

result=BFR.BE.EM.CV(x=X,v=V,b=B,q = 10)
