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
#' @example /examples/exampleCV.R
#'
#####General algorithm###
BFR.BE.EM.CV <- function(x,v=NULL,b=NULL,q=2,eps=0.001,it=100,epsM=0.05,
                         prior="N.MOM",varianceBE=TRUE,init=NULL,hyper=NULL,
                         seed=10,scaling=TRUE,intercept=FALSE,varimax=FALSE,
                         CV=TRUE,folds=10){

  set.seed(seed)

  if(CV==FALSE){
    model=BFR.BE.EM(x,v,b,q,eps,it,epsM,prior,varianceBE,init,hyper,seed,scaling,intercept,varimax)
    result_var=NULL
    result=NULL
  }else{
    #parameters
    n<-nrow(x)
    index.cv<-c(1:n)
    group<-list()
    #creating groups
    for (i in 1:folds){
      group[[i]]<-sample(index.cv,n/folds,replace=FALSE)
      index.cv<-index.cv[!index.cv%in%group[[i]]]
    }
    #CV
    print("CV varimax rotation")
    result_var<-llply(.data = 1:folds, .fun = function(y){
      print(y)
      if(is.matrix(v)){
        BFR.BE.EM(x=x[-group[[y]],], v=v[-group[[y]],], b=b[-group[[y]],],q,eps,it,
                  epsM,prior,varianceBE,init,hyper,seed,scaling,intercept,varimax=TRUE)
      }else{
        BFR.BE.EM(x=x[-group[[y]],], v=v[-group[[y]]], b=b[-group[[y]],],q,eps,it,
                  epsM,prior,varianceBE,init,hyper,seed,scaling,intercept,varimax=TRUE)
      }
    }, .parallel = FALSE)
    print("CV no rotation")
    result<-llply(.data = 1:folds, .fun = function(y){
      print(y)
      if(is.matrix(v)){
        BFR.BE.EM(x=x[-group[[y]],], v=v[-group[[y]],], b=b[-group[[y]],],q,eps,it,
                  epsM,prior,varianceBE,init,hyper,seed,scaling,intercept,varimax=FALSE)
      }else{
        BFR.BE.EM(x=x[-group[[y]],], v=v[-group[[y]]], b=b[-group[[y]],],q,eps,it,
                  epsM,prior,varianceBE,init,hyper,seed,scaling,intercept,varimax=FALSE)
      }
    }, .parallel = FALSE)
    #Choosing the model
    fnNR=mean(FN.CV.BE.weighted(result,X=x,V=v,B=b,group)$X.test)
    fnR=mean(FN.CV.BE.weighted(result_var,X=x,V=v,B=b,group)$X.test)
    fnNRjk=mean(FN.CV.BE.jk.weighted(result,X=x,V=v,B=b,group)$X.test)
    fnRjk=mean(FN.CV.BE.jk.weighted(result_var,X=x,V=v,B=b,group)$X.test)
    #Running chosen model
    model=BFR.BE.EM(x,v,b,q,eps,it,epsM,prior,varianceBE,init,hyper,seed,scaling,intercept,
                    varimax=(min(fnR,fnRjk)<=min(fnNR,fnNRjk)))
    if(min(fnNR,fnR)<=min(fnNRjk,fnRjk)){
      Mpost=PostM2(model)
    }else{
      Mpost=PostM(model)
    }
  }
  return(list(varimax=result_var,
              novarimax=result,
              group=group,
              M=model$M,
              sigma=model$sigma,
              Ez=model$Ez,
              Ezz=model$Ezz,
              Theta=model$Theta,
              gTheta = model$gTheta,
              gamma = model$gamma,
              D = model$D,
              iterations=model$iterations,
              Mpost=Mpost,
              rotation=(min(fnR,fnRjk)<=min(fnNR,fnNRjk)),
              columns=(min(fnR,fnNR)<=min(fnRjk,fnNRjk))))
}
