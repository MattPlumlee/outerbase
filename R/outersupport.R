#' outerbase
#' @aliases outerbase-package
#'  
#' @docType package
#' @import methods
#' @import Rcpp
#' @importFrom utils relist
#' @importFrom stats quantile sd
#' @useDynLib outerbase, .registration = TRUE
#' @exportPattern "^[[:alpha:]]+"
"_PACKAGE"
Rcpp::loadModule("obmod", TRUE)

#' BFGS standard
#' 
#' Do generic minimization of a function \code{funcw} that takes
#' a list par using the "Broyden-Fletcher-Goldfarb-Shanno" (BFGS) algorithm.
#' Useful for hyperparameter optimization because it handles infinite returns
#' fairly easily.  
#' 
#' @param funcw An object to optimize
#' @param par An initial point as a list
#' @param B An initial Hessian to start from
#' @param lr An initial learning rate to start from
#' @param ... additional parameters passed to \code{funcw}
#' @param verbose an integer from 0-3 where larger prints more information
#' @return a list of information from optimization, with the value stored in
#' `par`
#' @export
BFGS_std <- function(funcw, par, 
                     B = NULL, lr = 0.1, 
                     ...,verbose = 0){
  #linesearchparameters
  c1 = 0.0001
  c2 = 0.9
  numatte0 = 5
  
  parv = unlist(par)
  
  optid = funcw(relist(parv,par), ...)
  valo = optid$val
  go = unlist(optid$gval)
  resetB = TRUE
  
  if(is.null(B)) B = diag(1/sqrt(go^2+0.001))
  else resetB = FALSE
  
  twice = F
  lr0 = lr
  lr00 = lr
  numtimes = 0
  
  optdf = data.frame('iter.no'=0,'obj.value'=valo,
                     'wolfe.cond.1'=NA,'wolfe.cond.2'=NA,
                     'learning.rate'=lr)
  
  rownames(optdf) <- NULL
  if (verbose > 0) cat('\n', '########started BFGS#######', '\n', sep = "")
  if (verbose > 1) print(optdf[1,], row.names = FALSE, 
                         digits = 6)
  if(!sum(is.na(go))){
    for(k in 1:100){
      dirc = as.vector(-B %*% go)
      
      st = lr*dirc
      parvp = parv + st
      optid = funcw(relist(parvp,par), ...)
      wolfcond1 = (optid$val-valo) - c1*lr*(sum(dirc*go))
      wolfcond2 = -(sum(dirc*unlist(optid$gval)))+(c2*(sum(dirc*go)))
      
      numatte = numatte0
      lrlb = 0
      lrub = Inf
      lrh = lr
      wolfcond2o = wolfcond2
      numatte = numatte0
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
        parvp = parv + lrh*dirc
        optidh = funcw(relist(parvp,par), ...)
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
          lr0 = lr0 / 10
          lr = lr0
        }
        if(lr0 < lr00/(10^2+1)){
          break
        }
        optid = funcw(relist(parv,par), ...) #do not feed it extra info
        valo = optid$val
        go = unlist(optid$gval)
        B = diag(1/sqrt(0.001+go^2))
        resetB = TRUE
        if (verbose > 0) cat('restarted hessian','\n', sep = "")
        
        optdf[k+1,] <- data.frame('iter.no'=k,
                                  'obj.value'=NA,
                                  'wolfe.cond.1'=NA,
                                  'wolfe.cond.2'=NA,
                                  'learning.rate'=lr)
      } else {
        
        
        if(lr != lrh){
          lr = lrh
          st = parvp-parv
          parv = parvp
          optid = optidh
        } else{
          parv = parvp
        }
        
        if (k > 2 && (sum(st*go) > -length(go)/4) && twice){
          break
        } else if (k > 2 && (sum(st*go) > -length(go)/4)){
          twice = T
        }
        
        
        
        goo = go
        valo = optid$val
        go = unlist(optid$gval)
        yv = go-goo
        
        optdf[k+1,] <- data.frame('iter.no'=k,
                                  'obj.value'=valo,
                                  'wolfe.cond.1'=wolfcond1,
                                  'wolfe.cond.2'=wolfcond2,
                                  'learning.rate'=lr)
        if (verbose > 1) print(optdf[k+1,], row.names = FALSE, 
                               digits = 6)
        
        if (resetB) {
          B = sum(st*yv)/sum(yv*yv)* diag(length(parv))
          resetB = FALSE
        }
        
        cvh = 1/sum(st*yv)
        M1 = diag(length(go)) - cvh*outer(st,yv)
        B = M1 %*% B %*% t(M1) + cvh*outer(st,st)
        lr = lr^0.9 #drift toward 1
      }
    }} else{
      stop('initial gradient was undefined, stopping.')
    }
  
  optid = funcw(relist(parv,par), ...) #finish by evaluating
  if (verbose > 0){
    cat('num iter: ',k,'  ',
        'obj start: ', optdf[1,2],'  ',
        'obj end: ',   optid$val, '\n', sep = "")
    cat('final learning rate: ', lr, '\n', sep = "")
    cat('approx lower bound (not achieved): ',
        optid$val+sum(dirc*go), '\n', sep = "")
    if (verbose > 0) cat('#########finished BFGS########', '\n\n')
  } 
  list(par=relist(parv,par), B=B, lr=lr, optid=optid)
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
#' @param par an initial point, initialized from objects if needed
#' @param newt boolean for if Newtons method should be used
#' @param cgsteps max number of cg iterations, if \code{newt=FALSE}
#' @param cgtol cg tolerance, if \code{newt=FALSE}
#' @param ... additional parameters passed to \code{\link{BFGS_std}}
#' @return A list of information from optimization
#' @export
BFGS_lpdf <- function(om, logpdf, par=list(), newt=F, 
                      cgsteps=100, cgtol=0.001, ...){
  if(is.null(par$hyp)) par$hyp = gethyp(om)
  if(is.null(par$para)) par$para = getpara(logpdf)
  
  .lpdfwrapper(par, om, logpdf, newt=newt) #start by aligning para with actual
  optsum = BFGS_std(.lpdfwrapper, par, om=om, newt=newt,
                    logpdf=logpdf, ...)
  
  optsum
}

# Optimization wrapper
# 
# returns list of \code{val} and \code{grad}
.lpdfwrapper = function(parlist, om, logpdf, newt = F,
                        cgsteps=100, cgtol=0.001) {
  regpara = logpdf$paralpdf(parlist$para)
  reghyp = om$hyplpdf(parlist$hyp)
  if(is.finite(regpara) && is.finite(reghyp)) { 
    om$updatehyp(parlist$hyp)
    logpdf$updateom()
    logpdf$updatepara(parlist$para)
    if(newt) logpdf$optnewton()
    else logpdf$optcg(cgtol, cgsteps)
    
    gval = parlist
    gval$hyp = -logpdf$gradhyp-om$hyplpdf_grad(parlist$hyp)
    gval$para = -logpdf$gradpara-logpdf$paralpdf_grad(parlist$para)#
    list(val = -logpdf$val-reghyp-regpara, gval = gval)#
  } else list(val = Inf, gval = NULL) 
}