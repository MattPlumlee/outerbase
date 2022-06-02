#' outerbase
#' @aliases outerbase-package
#'  
#' @docType package
#' @import methods
#' @import Rcpp
#' @importFrom utils relist
#' @importFrom stats quantile
#' @useDynLib outerbase, .registration = TRUE
#' @exportPattern "^[[:alpha:]]+"
"_PACKAGE"
Rcpp::loadModule("obmod", TRUE)

#' BFGS standard
#' 
#' Do generic minimization of a function \code{funcw} that takes
#' a list rho using the Broyden-Fletcher-Goldfarb-Shanno algorithm.
#' Useful for hyperparameter optimization because it handles infs 
#' fairly easily.  
#' 
#' @param funcw An object to optimize
#' @param rho An initial point
#' @param ... additional parameters passed to \code{funcw}
#' @param verbose Integer from 0-3 where larger prints more information
#' @return A list of information from optimization
#' @export
BFGS_std <- function(funcw, rho, ...,verbose = 0){
  #linesearchparameters
  Bs = NULL #initial inverse hessian
  lr0 = 0.1 #initial learning rate
  c1 = 0.0001
  c2 = 0.9
  numatte0 = 5
  
  lr = lr0
  rhov = unlist(rho)
  
  optid = funcw(relist(rhov,rho), ...)
  valo = optid$val
  go = unlist(optid$gval)
  resetB = TRUE
  
  if(is.null(Bs)) Bs = diag(1/sqrt(go^2+0.001))
  else resetB = FALSE
  
  twice = F
  
  B = Bs
  lr00 = lr0
  numtimes = 0
  if (verbose > 0) print('doing opt...')
  if(!sum(is.na(go))){
    for(k in 1:100){
      dirc = as.vector(-B %*% go)
      
      st = lr*dirc
      rhovp = rhov + st
      optid = funcw(relist(rhovp,rho), ...)
      wolfcond1 = (optid$val-valo) - c1*lr*(sum(dirc*go))
      wolfcond2 = -(sum(dirc*unlist(optid$gval)))+(c2*(sum(dirc*go)))
      
      numatte = numatte0
      lrlb = 0
      lrub = Inf
      lrh = lr
      wolfcond2o = wolfcond2
      numatte = numatte0
      if (verbose > 1) {
        print('Wolfe conditions')
        print(c(wolfcond1,wolfcond2))
      }
      while(numatte > 0 && 
            ((is.na(wolfcond1) || is.na(wolfcond2)) ||
             ((wolfcond1>0) || (wolfcond2>0)))){
        if (is.na(wolfcond1) || is.na(wolfcond2) || wolfcond1>0){
          lrub = lrh
          lrh = 1/2*(lrlb+lrub)
        } else {
          lrlb = lrh
          if(is.finite(lrub)) lrh = 1/2*(lrlb+lrub)
          else lrh = 2*lrlb
        }
        rhovp = rhov + lrh*dirc
        optidh = funcw(relist(rhovp,rho), ...)
        wolfcond1 = (optidh$val-valo) - c1*lrh*(sum(dirc*go))
        wolfcond2 = -(sum(dirc*unlist(optidh$gval)))+(c2*(sum(dirc*go)))
        numatte = numatte-1
      }
      if (is.na(wolfcond1) || is.na(wolfcond2)){
        stop('something is very wrong... stuck on NAs')
      }
      if(wolfcond1>0){
        if(resetB){
          c2 = c2^0.5
          lr0 = lr0/10
          lr = lr0
        }
        if(lr0 < lr00/(10^2+1)){
          break
        }
        optid = funcw(relist(rhov,rho), ...) #do not feed it extra info
        valo = optid$val
        go = unlist(optid$gval)
        B = diag(1/sqrt(0.01+go^2))
        resetB = TRUE
      } else {
        if(lr != lrh){
          lr = lrh
          st = rhovp-rhov
          rhov = rhovp
          optid = optidh
        } else{
          rhov = rhovp
        }
        
        if (k > 2 && sum(st*go) > -length(go)/4 && twice){
          break
        } else if (k > 2 && sum(st*go) > -length(go)/4){
          twice=T
          resetB = TRUE
        }
        
        goo = go
        valo = optid$val
        go = unlist(optid$gval)
        yv = go-goo
        
        if(resetB){
          B = sum(st*yv)/sum(yv*yv)* diag(length(rhov))
          resetB = FALSE
        }
        
        cvh = 1/sum(st*yv)
        M1 = diag(length(go)) - cvh*outer(st,yv)
        B = M1 %*% B %*% t(M1) + cvh*outer(st,st)
        lr = lr * 1.05
      }
      if (verbose > 0) print(valo)
    }} else{
      stop('initial gradient was undefined, stopping.')
    }
  
  funcw(relist(rhov,rho), ...) #finish by evaluating
  if (verbose > 0) print('finished opt...')
  list(vec=relist(rhov,rho), B=B, lr=lr, optid=optid)
}

#' BFGS lpdf
#' 
#' A wrapper for code{\link{BFGS_std}} that is useful for easily calling 
#' parameter optimization for this package with as few lines as possible.
#' Note that \code{om} and \code{logpdf} will be set to optimal 
#' parameters, the return is simply for information.
#' 
#' @param om an \code{\link{outermod}} object
#' @param logpdf a \code{\link{lpdf}} object
#' @param rho an initial point, initialized from objects if needed
#' @param newt bool for if Newtons method should be used
#' @param ... additional parameters passed to \code{\link{BFGS_std}}
#' @return A list of information from optimization
#' @export
BFGS_lpdf <- function(om, logpdf, rho=list(), newt=F, ...){
  if(is.null(rho$hyp)) rho$hyp = gethyp(om)
  if(is.null(rho$para)) rho$para = getpara(logpdf)
  
  .lpdfwrapper(rho, om, logpdf, newt=newt)
  
  optsum = BFGS_std(.lpdfwrapper, rho, om=om, newt=newt,
                    logpdf=logpdf, ...)
  
  .lpdfwrapper(optsum$vec, om, logpdf, newt=newt)
  
  optsum
}

# Optimization wrapper
# 
# returns list of \code{val} and \code{grad}
.lpdfwrapper = function(hypl, om, logpdf, newt = F) {
  regpara = logpdf$paralpdf(hypl$para)
  reghyp = om$hyplpdf(hypl$hyp)
  if(is.finite(regpara) && is.finite(reghyp)) { 
    om$updatehyp(hypl$hyp)
    logpdf$updateom()
    logpdf$updatepara(hypl$para)
    if(newt) logpdf$optnewton()
    else logpdf$optcg(0.001, 100)
    
    gval = hypl
    gval$hyp = -logpdf$gradhyp-om$hyplpdf_grad(hypl$hyp)
    gval$para = -logpdf$gradpara-logpdf$paralpdf_grad(hypl$para)#
    list(val = -logpdf$val-reghyp-regpara, gval = gval)#
  } else list(val = Inf, gval = hypl) 
}