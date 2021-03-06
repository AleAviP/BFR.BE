#####E-step####
###Z Factors
E.Z <- function(x,M,sigma,A){
  if(is.vector(x)){
    x<-t(as.matrix(x))
  }
  n<-nrow(x)
  q<-ncol(M)
  aux.W<-t(M/sigma)%*%M+solve(A)
  W<-solve(aux.W)
  Z=W%*%t(M/sigma)%*%(t(x))
  Ezz.x <- Z%*%t(Z)+n*W
  return(list(z=t(Z),W=W,zz=Ezz.x))
}

###Gamma
##MoM
E.gamma.i.pMOM <- function(m.i,l0,l1,theta.i){
  dens0 = dnorm(m.i,0,sqrt(l0))
  dens1pMOM =dmom(m.i,l1)
  p.i =dens1pMOM*theta.i/(dens0*(1-theta.i)+dens1pMOM*theta.i)
  d.i = (1-p.i)/l0 +  p.i/l1
  return(list(p=p.i,d=d.i))
}

E.gamma.pMOM <- function(M,l0,l1,theta){
  p = nrow(M)
  Gamma.aux = llply(.data = 1:p, .fun = function(y){
    E.gamma.i.pMOM(M[y,],l0,l1,theta)
  }, .parallel = FALSE)
  return(Reduce(function(x,y) Map(rbind, x, y),Gamma.aux))
}

##Normal
E.gamma.i <- function(m.i,l0,l1,theta.i){
  dens0 = dnorm(m.i,0,sqrt(l0))
  dens1 = dnorm(m.i,0,sqrt(l1))
  p.i =dens1*theta.i/(dens0*(1-theta.i)+dens1*theta.i)
  d.i = (1-p.i)/l0 +  p.i/l1
  return(list(p=p.i,d=d.i))
}

E.gamma.2 <- function(M,l0,l1,theta){
  p = nrow(M)
  Gamma.aux = llply(.data = 1:p, .fun = function(y){
    #E.gamma.i(M[y,],diag(sigma)[y],l0,l1,theta)
    E.gamma.i(M[y,],l0,l1,theta)
  }, .parallel = FALSE)
  return(Reduce(function(x,y) Map(rbind, x, y),Gamma.aux))
}

#####M-step####
###M
##MoM
M.pMOM.cda.3<-function(x,Ez.x,Ezz.list,D,M.old,p.gamma,sigma,b){
  n <- nrow(x)
  p <- ncol(x)
  q <- ncol(M.old)
  p_b<-ncol(b)
  W.M <- 2*p.gamma
  if(!is.null(b)){
    XZaux = llply(.data = 1:p_b, .fun = function(y){
      sigma.j<-diag(1/sigma[,y])
      sigma.j%*%t(x[c(1:n)[b[,y]==1],])%*%Ez.x[c(1:n)[b[,y]==1],]
    }, .parallel = FALSE)
    xz = Reduce('+', XZaux)

    Maux = llply(.data = 1:p, .fun = function(y){
      Ezz.aux<-lapply(1:p_b, function(k) (Ezz.list[[k]])/sigma[y,k])
      Ezz.x = Reduce('+', Ezz.aux)
      b.j<-diag(D[y,])+Ezz.x
      aux_vec<-quadratic.solve(-b.j[1,1],(xz[y,1]-sum(M.old[y,-1]*Ezz.x[1,-1])),W.M[y,1])
      for(k in 2:q){
      #for(k in (q-1):1){
        aux<-quadratic.solve(-b.j[k,k],(xz[y,k]-sum(M.old[y,-k]*Ezz.x[k,-k])),W.M[y,k])
        aux_vec<-c(aux_vec,aux)
      }
      return(aux_vec)
    }, .parallel = FALSE)
  }
  else{
    xz <- t(x)%*%Ez.x
    Maux = llply(.data = 1:p, .fun = function(y){
      Ezz.x = Ezz.list
      b.j<-diag(sigma[y]*D[y,])+Ezz.x
      aux_vec<-quadratic.solve(-b.j[1,1],(xz[y,1]-sum(M.old[y,-1]*Ezz.x[1,-1])),sigma[y]*W.M[y,1])
      M.old[y,1]<-aux_vec
      for(k in 2:q){
        aux<-quadratic.solve(-b.j[k,k],(xz[y,k]-sum(M.old[y,-k]*Ezz.x[k,-k])),sigma[y]*W.M[y,k])
        M.old[y,k]<-aux
        aux_vec<-c(aux_vec,aux)}
      return(aux_vec)
    }, .parallel = FALSE)
  }
  Mnew = do.call("rbind", Maux)
  return(Mnew)
}

##Normal
M.new<-function(x,Ez.x,Ezz.list,D,b,sigma){
    n <- nrow(x)
    p <- ncol(x)
    p_b<-ncol(b)
    #q <- ncol(Ez.x)
    #xz <- t(x)%*%Ez.x
    if(!is.null(b)){
      XZaux = llply(.data = 1:p_b, .fun = function(y){
        sigma.j<-diag(1/sigma[,y])
        sigma.j%*%t(x[c(1:n)[b[,y]==1],])%*%Ez.x[c(1:n)[b[,y]==1],]
      }, .parallel = FALSE)
      xz = Reduce('+', XZaux)

      Maux = llply(.data = 1:p, .fun = function(y){
        Ezz.aux<-lapply(1:p_b, function(k) (Ezz.list[[k]])/sigma[y,k])
        Ezz.x = Reduce('+', Ezz.aux)
        xz[y,]%*%solve(diag(D[y,])+Ezz.x)
        #b.j<-diag(D[y,])+Ezz.x
        #sapply(1:q,function(k) (xz[y,k]-sum(M.old[y,-k]*b.j[k,-k]))/b.j[k,k])
      }, .parallel = FALSE)
    }
    else{
      xz <- t(x)%*%Ez.x
      Maux = llply(.data = 1:p, .fun = function(y){
        Ezz.x = Ezz.list
        xz[y,]%*%solve(diag(sigma[y]*D[y,])+Ezz.x)
        #b.j<-diag(D[y,])+Ezz.x
        #sapply(1:q,function(k) (xz[y,k]-sum(M.old[y,-k]*b.j[k,-k]))/b.j[k,k])
      }, .parallel = FALSE)
    }
    Mnew = do.call("rbind", Maux)
    #Mnew <- t(x)%*%Ez.x%*%solve(Ezz.x)
    return(Mnew)
}

##Flat
M.new.NS<-function(x,Ez.x,Ezz.x){
  M <- t(x)%*%Ez.x%*%solve(Ezz.x)
  return(M)
}

M.new.NS.2<-function(x,Ez.x,Ezz.list,b,sigma){
  n <- nrow(x)
  p <- ncol(x)
  p_b<-ncol(b)
  if(!is.null(b)){
    XZaux = llply(.data = 1:p_b, .fun = function(y){
      sigma.j<-diag(1/sigma[,y])
      sigma.j%*%t(x[c(1:n)[b[,y]==1],])%*%Ez.x[c(1:n)[b[,y]==1],]
    }, .parallel = FALSE)
    xz = Reduce('+', XZaux)

    Maux = llply(.data = 1:p, .fun = function(y){
      Ezz.aux<-lapply(1:p_b, function(k) (Ezz.list[[k]])/sigma[y,k])
      Ezz.x = Reduce('+', Ezz.aux)
      xz[y,]%*%solve(Ezz.x)
    }, .parallel = FALSE)
    Mnew = do.call("rbind", Maux)
  }
  else{
    #Ezz.x = Reduce('+', Ezz.list)
    Ezz.x = Ezz.list
    Mnew <- t(x)%*%Ez.x%*%solve(Ezz.x)
  }
  return(Mnew)
}

###Sigma
sigma.new<-function(x,Ez.x,Ezz.x,M,a,b){
  nl = nrow(x)
  diag((1/(nl+a-2))*diag(t(x)%*%x-2*t(x)%*%Ez.x%*%t(M)+M%*%Ezz.x%*%t(M)+a*b))
}

###Theta
Theta.new.prior<-function(x,v,b,Ez.x,M,sigma,sigmaT){
  p=nrow(M)
  q=ncol(M)
  n=nrow(x)
  p_v=ncol(v)
  if(!is.null(b)){
    p_b=ncol(b)
    theta_list = llply(.data = 1:p, .fun = function(y){
      vv_inv_aux<-lapply(1:p_b, function(k) (as.numeric(1/sigma[y,k])*(t(v[c(1:n)[b[,k]==1],])%*%v[c(1:n)[b[,k]==1],])))
      vv_inv<-solve(Reduce('+', vv_inv_aux)+diag(1/sigmaT,p_v))
      theta_aux<-lapply(1:p_b, function(k) (as.numeric(1/sigma[y,k])*(t(v[c(1:n)[b[,k]==1],])%*%(x[c(1:n)[b[,k]==1],y]-Ez.x[c(1:n)[b[,k]==1],]%*%M[y,]))))
      numerator<-Reduce('+', theta_aux)
      t(numerator)%*%vv_inv
    }, .parallel = FALSE)
    t(do.call("rbind", theta_list))
  }else{
    solve(t(v)%*%v+diag(1/sigmaT,p_v))%*%t(v)%*%(x-Ez.x%*%t(M))
  }
}

###Zeta
theta.gamma.new.2 <- function(P,a,b){
  p = nrow(P)
  q = ncol(P)
  j = c(1:q)
  theta = (apply(P,2,sum)+a/j-1)/(a/j+b+p-1)
  #Restrict values between 0 and 1
  theta[theta>1]<-1
  theta[theta<0]<-0
  return(theta)
}

#####Auxiliars####
###delta
delta<-function(a,b,c){
  b^2-4*a*c
}
##Solver
quadratic.solve<-function(a,b,c){
  if(delta(a,b,c) > 0){ # first case D>0
    pos = (-b+sqrt(delta(a,b,c)))/(2*a)
    neg = (-b-sqrt(delta(a,b,c)))/(2*a)
    #x = ifelse(abs(pos)>abs(neg),pos,neg)
    x = ifelse(b<0,pos,neg)
    result = x
  }
  else if(delta(a,b,c) == 0){ # second case D=0
    x = -b/(2*a)
    result =x
  }
  else {"There are no real roots."} # third case D<0
}
##lof
lof<-function(matrix,reference){
  matrix[,order(apply(reference,2,sum),decreasing=TRUE)]
}

######Logposterior####
likelihoodFA <- function(x,M,sigma){
  p = ncol(x)
  n = nrow(x)
  q = ncol(M)
  C = M %*% t(M) + sigma
  #l = -n/2 *(p*log(2 * pi )+log(det(C)) + tr(solve(C)%*%cov(x)))
  l = sum(dmvnorm(x, mean = rep(0, p), C, log = TRUE))
  return(l)
}

likelihoodFA2 <- function(x,mu,sigma){
  l = sum(dmnorm(x, mean = mu, sigma, log = TRUE))
  return(l)
}

logpriorSigma<- function(sigma,a,b){
  -sum((a+1)*log(sigma)+b/(sigma))
}

logpriorGamma<-function(gamma,at,bt,p,q){
  sum.gamma=apply(gamma,2,sum)
  #sum.gamma=rep(p,q)
  BetaBin.aux = llply(.data = 1:q, .fun = function(y){
    lbeta(sum.gamma[y]+at/y,p-sum.gamma[y]+bt)-lbeta(at/y,bt)
  }, .parallel = FALSE)
  BetaBin.aux = do.call("rbind", BetaBin.aux)
  #Stirling approximation
  BetaBin.approx = llply(.data = 1:q, .fun = function(y){
    0.5*log(2*pi)+(sum.gamma[y]+at/y-0.5)*log(sum.gamma[y]+at/y)+
      +(p-sum.gamma[y]+bt-0.5)*log(p-sum.gamma[y]+bt)-
      (p+at/y+bt-0.5)*log(p+at/y+bt)
  }, .parallel = FALSE)
  BetaBin.approx = do.call("rbind", BetaBin.approx)
  BetaBin.aux[is.infinite(BetaBin.aux)]=BetaBin.approx[is.infinite(BetaBin.aux)]
  l.gamma.gtheta=sum(BetaBin.aux)
  return(l.gamma.gtheta)
}

logpriorTheta<-function(theta,MuT,sigmaT){
  sum(dnorm(theta,MuT,sigmaT,log=TRUE))
}

logpriorFAall <- function(M,sigma,theta,gtheta,gamma,D,hyper,prior){
  q=ncol(M)
  p=nrow(M)
  l.sigma=logpriorSigma(sigma,hyper$aS/2,hyper$aS*hyper$bS/2)
  l.theta = logpriorTheta(theta,hyper$MuT,hyper$sigmaT)
  if(prior=="Flat"){
    l.M =0
    l.gamma.gtheta = 0
  }
  if(prior=="Normal"){
    # l.M = -.5*sum(M^2*D)
    l.M = sum(do.call("rbind",llply(.data = 1:p, .fun = function(y){
      -0.5*logdet(2*pi*diag((1-gamma[y,])*hyper$l0+gamma[y,]*hyper$l1))-0.5*sum(M[y,]^2*D[y,])
    }, .parallel = FALSE)))
    l.gamma.gtheta =logpriorGamma(gamma,hyper$at,hyper$bt,p,q)}
  if(prior=="N.MOM"){
    # l.M = -.5*sum(M^2*D) + sum(log(M^2)*gamma)
    l.M = sum(do.call("rbind",llply(.data = 1:p, .fun = function(y){
      -0.5*logdet(2*pi*diag((1-gamma[y,])*hyper$l0+gamma[y,]*hyper$l1))-0.5*sum(M[y,]^2*D[y,])
    }, .parallel = FALSE))) + sum(log(M^2)*gamma)
    l.gamma.gtheta =logpriorGamma(gamma,hyper$at,hyper$bt,p,q)}
  l.total=l.sigma+l.theta+l.gamma.gtheta+l.M
  return(list(lsigma=l.sigma,
              lTheta = l.theta,
              lGammaGtheta=l.gamma.gtheta,
              lM=l.M,
              lTotal=l.total))
}


#####Posterior Evaluation####
PostM<-function(FA){
  index<-FA$gamma>0.5
  q<-ncol(FA$M)
  p<-nrow(FA$M)
  M.new<-matrix(ncol=q,nrow=p,rep(0,p*q))
  if(sum(apply(index,2,sum)!=0)>2){
    for(k in 1:q){
      M.new[index[,k],k]<-FA$M[index[,k],k]
    }
  }else{
    M.new[,1:2]<-FA$M[,1:2]
  }
  return(M.new)
}

PostM2<-function(FA){
  index<-(apply(FA$gamma>0.5,2,sum)!=0)
  q<-ncol(FA$M)
  p<-nrow(FA$M)
  M.new<-matrix(ncol=q,nrow=p,rep(0,p*q))
  if(sum(index)>2){
    for(k in 1:q){
      if(index[k]==TRUE){M.new[,k]<-FA$M[,k]}
    }
  }else{
    M.new[,1:2]<-FA$M[,1:2]
  }
  return(M.new)
}

PostGamma<-function(FA){
  index<-(apply(FA$gamma>0.5,2,sum)!=0)
  q<-ncol(FA$M)
  p<-nrow(FA$M)
  Gamma.new<-matrix(ncol=q,nrow=p,rep(0,p*q))
  if(sum(index)>2){
    for(k in 1:q){
      if(index[k]==TRUE){Gamma.new[,k]<-FA$gamma[,k]}
    }
  }else{
    Gamma.new[,1:2]<-FA$M[,1:2]
  }
  return(Gamma.new)
}

PostD<-function(FA){
  index<-(apply(FA$gamma>0.5,2,sum)!=0)
  q<-ncol(FA$M)
  p<-nrow(FA$M)
  D.new<-matrix(ncol=q,nrow=p,rep(0,p*q))
  if(sum(index)>2){
    for(k in 1:q){
      if(index[k]==TRUE){D.new[,k]<-FA$D[,k]}
    }
  }else{
    D.new[,1:2]<-FA$M[,1:2]
  }
  return(D.new)
}


PostML<-function(FA){
  index<-(apply(FA$P_star>0.5,2,sum)!=0)
  q<-ncol(FA$B)
  p<-nrow(FA$B)
  M.new<-matrix(ncol=q,nrow=p,rep(0,p*q))
  if(sum(index)>2){
    for(k in 1:q){
      if(index[k]==TRUE){M.new[,k]<-FA$B[,k]}
    }
  }else{
    M.new[,1:2]<-FA$B[,1:2]
  }
  return(M.new)
}

PostMfun<-function(FA){
  gamma<-FA$gamma
  M<-FA$M
  index<-gamma>0.5
  q<-ncol(M)
  p<-nrow(M)
  M.new<-matrix(ncol=q,nrow=p,rep(0,p*q))
  if(sum(apply(index,2,sum)!=0)>2){
    for(k in 1:q){
      M.new[index[,k],k]<-M[index[,k],k]
    }
  }else{
    M.new[,1:2]<-M[,1:2]
  }
  return(M.new)
}

FN.CV.L.weighted<-function(FA,Y,Omega,B,group,q){
  k.cv<-length(FA)
  I.q <- diag(q)
  Z=llply(.data = 1:k.cv, .fun = function(y){
    E.Z(Y[group[[y]],],FA[[y]]$B,FA[[y]]$sigma,I.q)$z
    #E.Z(Y[group[[y]],],PostM2(FA[[y]]),FA[[y]]$sigma,I.q)$z
  }, .parallel = FALSE)
  FN.X.test=unlist(llply(.data = 1:k.cv, .fun = function(y){
    norm((Y[group[[y]],]-Z[[y]]%*%t(FA[[y]]$B))%*%diag(1/FA[[y]]$sigma),type ="F")
  }, .parallel = FALSE))
  FN.X.train=unlist(llply(.data = 1:k.cv, .fun = function(y){
    norm((Y[-group[[y]],]-FA[[y]]$Omega%*%t(FA[[y]]$B))%*%diag(1/FA[[y]]$sigma),type ="F")
  }, .parallel = FALSE))
  return(list(X.test=FN.X.test,
              X.train=FN.X.train))
}

FN.CV.L<-function(FA,Y,Omega,B,group,q){
  k.cv<-length(FA)
  I.q <- diag(q)
  Z=llply(.data = 1:k.cv, .fun = function(y){
    E.Z(Y[group[[y]],],FA[[y]]$B,FA[[y]]$sigma,I.q)$z
    #E.Z(Y[group[[y]],],PostM2(FA[[y]]),FA[[y]]$sigma,I.q)$z
  }, .parallel = FALSE)
  FN.X.test=unlist(llply(.data = 1:k.cv, .fun = function(y){
    norm(Y[group[[y]],]-Z[[y]]%*%t(FA[[y]]$B),type ="F")
  }, .parallel = FALSE))
  FN.X.train=unlist(llply(.data = 1:k.cv, .fun = function(y){
    norm(Y[-group[[y]],]-FA[[y]]$Omega%*%t(FA[[y]]$B),type ="F")
  }, .parallel = FALSE))
  return(list(X.test=FN.X.test,
              X.train=FN.X.train))
}

q.hat.fun<-function(FA){
  k.cv<-length(FA)
  if(is.null(FA[[1]]$gamma)){
    if(is.null(FA[[1]]$v)){
      unlist(lapply(1:k.cv, function(k) (sum(apply(FA[[k]]$P_star>.5,2,sum)!=0))))
    }else{
      unlist(lapply(1:k.cv, function(k) (sum(apply(FA[[k]]$v,2,sum)!=0))))
    }
  }else{
    unlist(lapply(1:k.cv, function(k) (sum(apply(FA[[k]]$gamma>.5,2,sum)!=0))))
  }
}

m.hat.fun<-function(FA,q=10){
  k.cv<-length(FA)
  if(is.null(FA[[1]]$gamma)){
    if(is.null(FA[[1]]$v)){
      unlist(lapply(1:k.cv, function(k) (sum(apply(FA[[k]]$P_star>.5,2,sum)))/q))
    }else{
      unlist(lapply(1:k.cv, function(k) (sum(apply(FA[[k]]$v,2,sum)))/q))
    }
  }else{
    unlist(lapply(1:k.cv, function(k) (sum(apply(FA[[k]]$gamma>.5,2,sum)))/q))
  }
}

it.fun<-function(FA,q=10){
  k.cv<-length(FA)
  if(is.null(FA[[1]]$gamma)){
    if(is.null(FA[[1]]$v)){
      unlist(lapply(1:k.cv, function(k) (FA[[k]]$niter)))
    }else{
      rep(5,k.cv)
    }
  }else{
    unlist(lapply(1:k.cv, function(k) (FA[[k]]$iterations)))
  }
}

FN.CV.BE.weighted<-function(FA,X,V,B,group){
  q <- ncol(FA[[1]]$M)
  k.cv <- length(FA)
  I.q <- diag(q)
  n <- nrow(X)
  if(!is.null(B)){
    p_b <- ncol(B)
    EZ=llply(.data = 1:k.cv, .fun = function(y){
      #print(y)
      if(is.matrix(V)){
        X.hat<-X[group[[y]],]-cbind(V[group[[y]],],B[group[[y]],])%*%t(FA[[y]]$Theta)
      }else{
        X.hat<-X[group[[y]],]-cbind(V[group[[y]]],B[group[[y]],])%*%t(FA[[y]]$Theta)
      }
      B.group<-B[group[[y]],]
      n.test<-nrow(B.group)
      index.aux<-c(1:n.test)
      index.aux.sort<-unlist(lapply(1:p_b, function(k) index.aux[B.group[,k]==1]))
      z.list<-lapply(1:p_b, function(k) E.Z(X.hat[c(1:n.test)[B.group[,k]==1],],FA[[y]]$M,FA[[y]]$sigma[,k],I.q)$z)
      z<-Reduce('rbind', z.list)
      z<-z[order(index.aux.sort),]
      return(z)
    }, .parallel = FALSE)
    FN.X.test=unlist(llply(.data = 1:k.cv, .fun = function(y){
      if(is.matrix(V)){
        X.dif<-X[group[[y]],]-
          (EZ[[y]]%*%t(PostM2(FA[[y]]))+cbind(V[group[[y]],],B[group[[y]],])%*%t(FA[[y]]$Theta))
      }else{
        X.dif<-X[group[[y]],]-
          (EZ[[y]]%*%t(PostM2(FA[[y]]))+cbind(V[group[[y]]],B[group[[y]],])%*%t(FA[[y]]$Theta))
      }
      B.group<-B[group[[y]],]
      n.test<-nrow(B.group)
      norm.test<-lapply(1:p_b, function(k)
        norm(X.dif[c(1:n.test)[B.group[,k]==1],]%*%diag(1/FA[[y]]$sigma[,k]),type ="F"))
      norm.test<-Reduce('+', norm.test)
      return(norm.test)
    }, .parallel = FALSE))
    FN.X.train=unlist(llply(.data = 1:k.cv, .fun = function(y){
      if(is.matrix(V)){
        X.dif<-X[-group[[y]],]-
          (FA[[y]]$Ez%*%t(PostM2(FA[[y]]))+cbind(V[-group[[y]],],B[-group[[y]],])%*%t(FA[[y]]$Theta))
      }else{
        X.dif<-X[-group[[y]],]-
          (FA[[y]]$Ez%*%t(PostM2(FA[[y]]))+cbind(V[-group[[y]]],B[-group[[y]],])%*%t(FA[[y]]$Theta))
      }
      B.group<-B[-group[[y]],]
      n.train<-nrow(B.group)
      norm.train<-lapply(1:p_b, function(k)
        norm(X.dif[c(1:n.train)[B.group[,k]==1],]%*%diag(1/FA[[y]]$sigma[,k]),type ="F"))
      norm.train<-Reduce('+', norm.train)
      return(norm.train)
    }, .parallel = FALSE))
  }else{
    Z=llply(.data = 1:k.cv, .fun = function(y){
      E.Z(X[group[[y]],],FA[[y]]$M,FA[[y]]$sigma,I.q)$z
    }, .parallel = FALSE)
    FN.X.test=unlist(llply(.data = 1:k.cv, .fun = function(y){
      norm((X[group[[y]],]-Z[[y]]%*%t(PostM2(FA[[y]])))%*%diag(1/FA[[y]]$sigma),type ="F")
    }, .parallel = FALSE))
    FN.X.train=unlist(llply(.data = 1:k.cv, .fun = function(y){
      norm((X[-group[[y]],]-FA[[y]]$Ez%*%t(PostM2(FA[[y]])))%*%diag(1/FA[[y]]$sigma),type ="F")
    }, .parallel = FALSE))
  }
  return(list(X.test=FN.X.test,
              X.train=FN.X.train))
}

FN.CV.BE.jk.weighted<-function(FA,X,V,B,group){
  q <- ncol(FA[[1]]$M)
  k.cv <- length(FA)
  I.q <- diag(q)
  n <- nrow(X)
  if(!is.null(B)){
    p_b <- ncol(B)
    EZ=llply(.data = 1:k.cv, .fun = function(y){
      #print(y)
      if(is.matrix(V)){
        X.hat<-X[group[[y]],]-cbind(V[group[[y]],],B[group[[y]],])%*%t(FA[[y]]$Theta)
      }else{
        X.hat<-X[group[[y]],]-cbind(V[group[[y]]],B[group[[y]],])%*%t(FA[[y]]$Theta)
      }
      B.group<-B[group[[y]],]
      n.test<-nrow(B.group)
      index.aux<-c(1:n.test)
      index.aux.sort<-unlist(lapply(1:p_b, function(k) index.aux[B.group[,k]==1]))
      z.list<-lapply(1:p_b, function(k) E.Z(X.hat[c(1:n.test)[B.group[,k]==1],],FA[[y]]$M,FA[[y]]$sigma[,k],I.q)$z)
      z<-Reduce('rbind', z.list)
      z<-z[order(index.aux.sort),]
      return(z)
    }, .parallel = FALSE)
    FN.X.test=unlist(llply(.data = 1:k.cv, .fun = function(y){
      if(is.matrix(V)){
        X.dif<-X[group[[y]],]-
          (EZ[[y]]%*%t(PostM(FA[[y]]))+cbind(V[group[[y]],],B[group[[y]],])%*%t(FA[[y]]$Theta))
      }else{
        X.dif<-X[group[[y]],]-
          (EZ[[y]]%*%t(PostM(FA[[y]]))+cbind(V[group[[y]]],B[group[[y]],])%*%t(FA[[y]]$Theta))
      }
      B.group<-B[group[[y]],]
      n.test<-nrow(B.group)
      norm.test<-lapply(1:p_b, function(k)
        norm(X.dif[c(1:n.test)[B.group[,k]==1],]%*%diag(1/FA[[y]]$sigma[,k]),type ="F"))
      norm.test<-Reduce('+', norm.test)
      return(norm.test)
    }, .parallel = FALSE))
    FN.X.train=unlist(llply(.data = 1:k.cv, .fun = function(y){
      if(is.matrix(V)){
        X.dif<-X[-group[[y]],]-
          (FA[[y]]$Ez%*%t(PostM(FA[[y]]))+cbind(V[-group[[y]],],B[-group[[y]],])%*%t(FA[[y]]$Theta))
      }else{
        X.dif<-X[-group[[y]],]-
          (FA[[y]]$Ez%*%t(PostM(FA[[y]]))+cbind(V[-group[[y]]],B[-group[[y]],])%*%t(FA[[y]]$Theta))
      }
      B.group<-B[-group[[y]],]
      n.train<-nrow(B.group)
      norm.train<-lapply(1:p_b, function(k)
        norm(X.dif[c(1:n.train)[B.group[,k]==1],]%*%diag(1/FA[[y]]$sigma[,k]),type ="F"))
      norm.train<-Reduce('+', norm.train)
      return(norm.train)
    }, .parallel = FALSE))
  }else{
    Z=llply(.data = 1:k.cv, .fun = function(y){
      E.Z(X[group[[y]],],FA[[y]]$M,FA[[y]]$sigma,I.q)$z
    }, .parallel = FALSE)
    FN.X.test=unlist(llply(.data = 1:k.cv, .fun = function(y){
      norm((X[group[[y]],]-Z[[y]]%*%t(PostM(FA[[y]])))%*%diag(1/FA[[y]]$sigma),type ="F")
    }, .parallel = FALSE))
    FN.X.train=unlist(llply(.data = 1:k.cv, .fun = function(y){
      norm((X[-group[[y]],]-FA[[y]]$Ez%*%t(PostM(FA[[y]])))%*%diag(1/FA[[y]]$sigma),type ="F")
    }, .parallel = FALSE))
  }
  return(list(X.test=FN.X.test,
              X.train=FN.X.train))
}

Ez.test.ov<-function(X,V,B,FA,group,q){
  k.cv <- length(FA)
  I.q <- diag(q)
  p_b <- ncol(B)
  n <- nrow(X)
  EZ=llply(.data = 1:k.cv, .fun = function(y){
    #print(y)
    X.hat<-X[group[[y]],]-cbind(V[group[[y]]],B[group[[y]],])%*%t(FA[[y]]$Theta)
    B.group<-B[group[[y]],]
    n.test<-nrow(B.group)
    index.aux<-c(1:n.test)
    index.aux.sort<-unlist(lapply(1:p_b, function(k) index.aux[B.group[,k]==1]))
    #z.list<-lapply(1:p_b, function(k) E.Z(X.hat[c(1:n.test)[B.group[,k]==1],],PostM(FA[[y]]),FA[[y]]$sigma[,k],I.q)$z)
    z.list<-lapply(1:p_b, function(k) E.Z(X.hat[c(1:n.test)[B.group[,k]==1],],FA[[y]]$M,FA[[y]]$sigma[,k],I.q)$z)
    z<-Reduce('rbind', z.list)
    z<-z[order(index.aux.sort),]
    return(z)
  }, .parallel = FALSE)
}

Ez.test.COMBAT.ov<-function(X,FA,group){
  k.cv <- length(FA)
  n <- nrow(X)
  EZ=llply(.data = 1:k.cv, .fun = function(y){
    q<-ncol(FA[[y]]$Gamma)
    tilde.Gamma <- sqrt(1/FA[[y]]$Sigma) * FA[[y]]$Gamma
    M <- diag(q) + t(tilde.Gamma) %*% tilde.Gamma
    eigenM <- eigen(M, T)
    YSG <- X[group[[y]],] %*% (1/FA[[y]]$Sigma * FA[[y]]$Gamma)
    varZ <- eigenM$vectors %*% (1/eigenM$values * t(eigenM$vectors))
    z <- YSG %*% varZ
    return(z)
  }, .parallel = FALSE)
}

Ez.test.COMBAT.ov.2<-function(X,FA,group){
  k.cv <- length(FA)
  n <- nrow(X)
  EZ=llply(.data = 1:k.cv, .fun = function(y){
    q<-ncol(FA[[y]]$Gamma)
    tilde.Gamma <- sqrt(1/FA[[y]]$Sigma) * FA[[y]]$Gamma
    M <- diag(q) + t(tilde.Gamma) %*% tilde.Gamma
    eigenM <- eigen(M, T)
    YSG <- tryCatch(t(X[[y]]) %*% (1/FA[[y]]$Sigma * FA[[y]]$Gamma),
                    error = function(e) cat("Error for ", y))
    varZ <- eigenM$vectors %*% (1/eigenM$values * t(eigenM$vectors))
    z <- tryCatch(YSG %*% varZ, error = function(e) cat(""))
    return(z)
  }, .parallel = FALSE)
}

FNorm<-function(E.z,M_hat,Z,M){
  norm(E.z%*%t(M_hat)-Z%*%t(M),type ="F")
}

FNorm2<-function(Sigma_hat,M_hat,Sigma){
  C_hat = M_hat%*%t(M_hat)+diag(Sigma_hat)
  C = Sigma
  FNorm = norm(C_hat-C,type ="F")
  return(FNorm)
}

FNorm2BE<-function(Sigma_hat,M_hat,Sigma){
  C = Sigma
  if(is.vector(Sigma_hat)){
    C_hat = M_hat%*%t(M_hat)+diag(Sigma_hat)
  }else{
    C_hat = M_hat%*%t(M_hat)+diag(apply(Sigma_hat,1,sum))
  }
  FNorm = norm(C_hat-C,type ="F")
  return(FNorm)
}

FNorm3<-function(E.z,M_hat,x){
  #round(M_hat,1)
  C_hat = cov(E.z%*%t(M_hat))
  C = cov(x)
  FNorm = norm(C_hat-C,type ="F")
  return(FNorm)
}

FNorm4<-function(E.z,M_hat,Z,M){
  #round(M_hat,1)
  C_hat = cov(E.z%*%t(M_hat))
  C = cov(Z%*%t(M))
  FNorm = norm(C_hat-C,type ="F")
  return(FNorm)
}

FNormXhat<-function(FA,z,v,b,Theta,Beta,M,Postj=NULL){
  x.hat = z%*%t(M)+v%*%t(Theta)+b%*%t(Beta)
  if(is.null(Postj)){
    X.rec = FA$Ez%*%t(FA$M)+cbind(v,b)%*%t(FA$Theta)
  }
  else{
    if(Postj==TRUE){
      X.rec = FA$Ez%*%t(PostM2(FA))+cbind(v,b)%*%t(FA$Theta)
    }else{X.rec = FA$Ez%*%t(PostM(FA))+cbind(v,b)%*%t(FA$Theta)}
  }
  FNorm = norm(X.rec-x.hat,type ="F")
  return(FNorm)
}

FNormXhatL<-function(Ez,M_hat,z,v,b,Theta,Beta,M){
  x.hat = z%*%t(M)+v%*%t(Theta)+b%*%t(Beta)
  X.rec = Ez%*%t(M_hat)
  FNorm = norm(X.rec-x.hat,type ="F")
  return(FNorm)
}

FNormXhatC<-function(Ez,M_hat,X){
  x.hat = X
  X.rec = Ez%*%t(M_hat)
  FNorm = norm(X.rec-x.hat,type ="F")
  return(FNorm)
}

#Ez.std
FA.ov<-function(FA,B,q=100){
  p_b<-ncol(B)
  n<-nrow(FA$Ez)
  FA_z<-Reduce('rbind', llply(.data = 1:p_b, .fun = function(y){
      FA$Ez[c(1:n)[B[,y]==1],]%*%solve(t(FA$M/FA$sigma[,y])%*%FA$M+diag(q))
    }, .parallel = FALSE))
}

#Scale X
x.scale<-function(x){
  p<-ncol(x)
  x_var_col<-apply(x,2,sd)
  x_mean_col<-apply(x,2,mean)
  x_scale=Reduce("cbind",llply(.data = 1:p, .fun = function(k){
    (x[,k]-x_mean_col[k])/x_var_col[k]
  }, .parallel = FALSE))
  return(x_scale)
}
