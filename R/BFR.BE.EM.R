#' Bayesian factor regression with batch effect correction.
#'
#' @param x n x p matrix of p observations for n individuals
#' @param v n x p_v matrix of known covariates.
#' @param b n x p_b matrix of batch indicator.
#' @param q latent cardinality.
#' @param eps margin for convergence.
#' @param it maximum number of iterations.
#' @param epsM margin for convergence for M.
#' @param prior Spike prior to use for the laodings: "N.MOM", "Normal" or "Flat".
#' @param varianceBE TRUE if variance batch correction is modeled; FALSE if not.
#' @param init list of initial parameters.
#' @param hyper list of hyper parameters.
#' @param seed random seed to use.
#' @param scaling mean centering the data?
#' @param intercept adds intercept to the Bayesian factor regression.
#' @param varimax TRUE for initial varimax rotated loadings.
#' @return A list with the optimised parameters.
#' @example /examples/example.R
#'
#####General algorithm###
BFR.BE.EM <- function(x,v=NULL,b=NULL,q=2,eps=0.001,it=100,epsM=0.05,
                      prior="N.MOM",varianceBE=TRUE,init=NULL,hyper=NULL,
                      seed=4,scaling=FALSE,intercept=FALSE,varimax=FALSE){

  set.seed(seed)

  if (scaling==TRUE) x<-scale(x)

  if(is.null(hyper)){
    if(prior=="N.MOM"){
      l0<-.1/qnorm(.975)^2
      l1<-0.284215
    }else{
      l0<-.1/qnorm(.975)^2
      l1<-0.8526451
    }
    hyper=list(theta=0.5,
               aS = 1,
               bS = 1,
               l0 = l0,
               l1 = l1,
               at = 1,
               bt = 1,
               MuT = 0,
               sigmaT = 1)
  }

  n<-nrow(x)
  p<-ncol(x)
  p_b<-1
  w_b<-1
  sorted.index<-c(1:n)
  if (!is.null(b)) {p_b=ncol(b)
  n_b=apply(b,2,sum)
  w_b=n_b/n
  sorted.index.aux <- 1
  for(i in 1:p_b){
    sorted.index.aux <- c(sorted.index.aux,sorted.index[b[,i]==1])
  }
  sorted.index<-sorted.index.aux[-1]
  if (intercept==FALSE){
    if(!is.null(v)){
      v=cbind(v,b)
    }else{
      v=b
    }
  }
  }
  if (intercept==TRUE){
    if(!is.null(v)){
      v_0=rep(1,n)
      if(!is.null(b)){
        v=cbind(v_0,v,b[,-p_b])
      }else{
        v=cbind(v_0,v)
      }
    }
  }
  if(!is.null(v)){
    p_v<-ncol(v)
  }else{p_v<-1}

  #Hyperparameters
  g.theta = hyper$theta
  aS = hyper$aS
  bS = hyper$bS
  l0 = hyper$l0
  l1 = hyper$l1
  at = hyper$at
  bt = hyper$bt
  MuT = hyper$MuT
  sigmaT = hyper$sigmaT

  ## initialization of M sigma and theta
  if (is.null(init)) {
    if(!is.null(v)){
      fm1 <- lm(x ~ v)
      theta <- t(fm1$coefficients[-1,])
      #Replacing NA for zeros
      theta[is.na(theta)] <- 0
      svd.X <- svd(fm1$residuals)
      svd.U <- svd.X$u
      svd.V <- svd.X$v
      Lambda <- diag(svd.X$d)
      if(n>=q){
        sigma2 <- 1/(p-q) * sum(svd.U[(q+1):p])                                                          # variance; average variance associated with discarded dimensions
        aux_p <- min(n,p)
        if(is.positive.definite(Lambda-sigma2*diag(aux_p))==FALSE){
          aux.chol <- matrix(ncol=n,nearPD(Lambda-sigma2*diag(aux_p))$mat)
        }else{aux.chol <- Lambda-sigma2*diag(aux_p)}
        M<-(svd.V %*% chol(aux.chol) %*% diag(aux_p))[,1:q]
      }else{
        #Reconstruction x=svd.U%*%Lambda%*%t(svd.V)%*%svd.V%*%t(svd.V) PC=svd.U%*%Lambda%*%t(svd.V) M'=svd.V%*%t(svd.V)
        M<-svd.V%*%t(svd.V)[,1:q]
      }
    }else{
      theta <- rep(0,p)
      svd.X <- svd(x)
      svd.U <- svd.X$u
      svd.V <- svd.X$v
      Lambda <- diag(svd.X$d)
      if(n>=q){
        sigma2 <- 1/(p-q) * sum(svd.U[(q+1):p])                                                          # variance; average variance associated with discarded dimensions
        aux_p <- min(n,p)
        if(is.positive.definite(Lambda-sigma2*diag(aux_p))==FALSE){
          aux.chol <- matrix(ncol=n,nearPD(Lambda-sigma2*diag(aux_p))$mat)
        }else{aux.chol <- Lambda-sigma2*diag(aux_p)}
        M<-(svd.V %*% chol(aux.chol) %*% diag(aux_p))[,1:q]
      }else{
        #Reconstruction x=svd.U%*%Lambda%*%t(svd.V)%*%svd.V%*%t(svd.V) PC=svd.U%*%Lambda%*%t(svd.V) M'=svd.V%*%t(svd.V)
        M<-svd.V%*%t(svd.V)[,1:q]
      }
    }
    if(varimax==TRUE){
      M<-varimax(M+0.0000001)$loadings
    }
    z<-mvrnorm(n = n, rep(0,q), diag(q))
    sigma<-diag(p)
    #g.theta <- rep(g.theta,p)
    g.theta <- rep(g.theta,q)
    if(prior=="N.MOM"){
      gamma <- E.gamma.pMOM(M,l0,l1,g.theta)
    }else{
      gamma <- E.gamma.2(M,l0,l1,g.theta)
    }
    D <- gamma$d
    p.gamma <- gamma$p

  } else {
    theta<-init$theta
    sigma<-init$sigma
    M<-init$M
    z<-init$z
    D <- init$D
    p.gamma <- init$p.gamma
    g.theta <- init$g.theta}

  ##traces
  trace_sigma <- matrix(nrow =1, ncol = p*p_b,diag(sigma))
  trace_theta <- matrix(nrow =1, ncol = p*p_v,c(theta))
  trace_M <- matrix(nrow = 1, ncol = p*q,c(M))
  trace_D <- matrix(nrow = 1, ncol = p*q,c(D))
  trace_p <- matrix(nrow = 1, ncol = p*q,c(p.gamma))
  trace_g.theta <- matrix(nrow = 1, ncol = q,c(g.theta))

  #tracesaux
  trace_priorGammaGtheta <- 0
  trace_priorMpMOM <- 0

  #####EM#####
  count<-1
  likelihood<-NULL
  change<-1000
  changeM<-1000
  wait<-0
  #lik<--100000000000
  if(!is.null(v)){
    x_o <- x
    x<-x_o-v%*%t(theta)
    tv <- t(v)
    I.pv <- diag(p_v)
    }
  I.q <- diag(q)

  likelihood<-likelihoodFA(x,M,sigma)
  logprior_aux<-logpriorFAall(M,diag(sigma),theta,g.theta,p.gamma,D,hyper,prior)
  logprior<-logprior_aux$lTotal
  logpriorM <- logprior_aux$lM
  logpriorsigma <- logprior_aux$lsigma
  logpriorTheta <- logprior_aux$lTheta
  logpriorGammaGtheta <- logprior_aux$lGammaGtheta
  logpost<-logprior + likelihood
  lik<-logpost

  trace_like <- likelihood
  trace_like2 <- likelihood
  trace_prior <- logprior
  trace_post <- logpost
  trace_priorM <- logpriorM
  trace_priorsigma <- logpriorsigma
  trace_priorTheta <- logpriorTheta
  if(!prior=="Flat"){
    trace_priorGammaGtheta <- logpriorGammaGtheta
  }

  if(varianceBE==FALSE){
    sigma_aux=rep(1,p)
  }else{sigma_aux=matrix(1,p,p_b)}

  while ((count < it) & ((changeM > epsM ) & (change > eps ))) {
    #print(count)
    ##E step
    if (varianceBE==FALSE){
      Ez<-E.Z(x,as.matrix(M),sigma_aux,I.q)
      Ez.x<-Ez$z
      Ezz.x<-Ez$zz
      Ezz.list<-Ezz.x
    }
    if (varianceBE==TRUE){
      Ez.x<-c(1:q)
      Ezz.list<-list()
      for (i in 1:p_b){
        Ez<-E.Z(x[c(1:n)[b[,i]==1],],as.matrix(M),sigma_aux[,i],I.q)
        Ez.x<-rbind(Ez.x,Ez$z)
        Ezz.list[[i]]<-Ez$zz
      }
      #sort them
      Ez.x<-Ez.x[-1,]
      Ez.x<-Ez.x[order(sorted.index),]
      #Add Ezz
      Ezz.x = Reduce('+', Ezz.list)
    }

    if(prior=="N.MOM"){
      gamma <- E.gamma.pMOM(M,l0,l1,g.theta)
      D <- gamma$d
      p.gamma <- gamma$p
    }
    if(prior=="Normal"){
      gamma <- E.gamma.2(M,l0,l1,g.theta)
      D <- gamma$d
      p.gamma <- gamma$p
    }

    ##M-step
    M.old<-M
    if(prior=="N.MOM"){
      M<-M.pMOM.cda.3(x,Ez.x,Ezz.list,D,M.old,p.gamma,sigma_aux,b)
    }
    if(prior=="Normal"){
      M<-M.new(x,Ez.x,Ezz.list,D,b,sigma_aux)
    }
    if(prior=="Flat"){
      M<-M.new.NS.2(x,Ez.x,Ezz.list,b,sigma_aux)
    }

    tM <- t(M)

    if (varianceBE==FALSE){
      sigma<-sigma.new(x,Ez.x,Ezz.x,M,aS,bS)
      sigma_aux<-diag(sigma)
    }else{
      for (i in 1:p_b){
        x_b<-x[c(1:n)[b[,i]==1],]
        z_b<-Ez.x[c(1:n)[b[,i]==1],]
        sigmaInv <- diag(1/sigma_aux[,i])
        W_b<-solve(I.q + tM%*%sigmaInv%*%M)
        zz_b<-t(z_b)%*%z_b+n_b[i]*W_b
        sigma_aux[,i]<-diag(sigma.new(x_b,z_b,zz_b,M,aS,bS))
      }
      sigma<-diag(c(sigma_aux %*% w_b))
    }

    if(!is.null(v)){
      theta<-Theta.new.prior(x_o,v,b,Ez.x,M,sigma_aux,sigmaT)
      x<-x_o-v%*%theta
      trace_theta <- rbind(trace_theta,c(t(theta)))
    }

    if(!prior=="Flat"){
      g.theta <- theta.gamma.new.2(p.gamma,at,bt)
    }

    #Saving values
    trace_sigma <- rbind(trace_sigma,c(sigma_aux))
    trace_M <- rbind(trace_M,c(M))
    if(!prior=="Flat"){
      trace_D <- rbind(trace_D,c(D))
      trace_p <- rbind(trace_p,c(p.gamma))
      trace_g.theta <- rbind(trace_g.theta,c(g.theta))
    }

    ##Likelihood
    #likelihood<-likelihoodFA(x,M,sigma)
    if (varianceBE==FALSE){
      likelihood<-likelihoodFA(x,M,sigma)
      likelihood2<-likelihood
    }else{
      like.aux<-llply(.data = 1:p_b, .fun = function(y){
        likelihoodFA(x[c(1:n)[b[,y]==1],],M,diag(sigma_aux[,y]))
        #likelihoodFA2(x_o[c(1:n)[b[,y]==1],],v[c(1:n)[b[,y]==1],]%*%theta+Ez.x[c(1:n)[b[,y]==1],]%*%t(M),diag(sigma_aux[,y]))
      }, .parallel = FALSE)
      likelihood<-sum(do.call("rbind", like.aux))
      likelihood2<-likelihood
    }

    logprior_aux<-logpriorFAall(M,sigma_aux,theta,g.theta,p.gamma,D,hyper,prior)

    logprior<-logprior_aux$lTotal
    logpriorM <- logprior_aux$lM
    logpriorsigma <- logprior_aux$lsigma
    logpriorTheta <- logprior_aux$lTheta
    logpriorGammaGtheta <- logprior_aux$lGammaGtheta
    logpost <- likelihood+logprior
    #change<-likelihood-lik
    #lik<-likelihood
    #change<-abs(logpost-lik)
    #change<-logpost-lik
    change<-abs(logpost-lik)
    changeM<-max(abs(M.old-M))
    #change<-logpost-lik
    #if(change<0) change<-(eps+1)
    lik<-logpost
    trace_like<-c(trace_like,likelihood)
    trace_like2<-c(trace_like2,likelihood2)
    trace_prior<-c(trace_prior,logprior)
    trace_post<-c(trace_post,logpost)
    trace_priorM <- c(trace_priorM,logpriorM)
    trace_priorsigma <- c(trace_priorsigma,logpriorsigma)
    trace_priorTheta <- c(trace_priorTheta,logpriorTheta)
    trace_priorGammaGtheta <- c(trace_priorGammaGtheta,logpriorGammaGtheta)
    count<-count+1
    #print(count)
    #print(logpriorM)
    #print(round(change,3))
    #print(round(changeM,3))

  }

  return(list(M=M,
              sigma=sigma_aux,
              Ez=Ez.x,
              Ezz=Ezz.x,
              Theta=t(theta),
              gTheta = g.theta,
              gamma = p.gamma,
              D = D,
              tracesigma=trace_sigma,
              traceTheta=trace_theta,
              traceM=trace_M,
              traceD=trace_D,
              traceGtheta=trace_g.theta,
              traceP=trace_p,
              iterations=count,
              like=trace_like,
              like2=trace_like2,
              prior=trace_prior,
              post=trace_post,
              priorM = trace_priorM,
              priorsigma = trace_priorsigma,
              priorGammaGtheta = trace_priorGammaGtheta,
              priorTheta = trace_priorTheta))
}
