#' list all covariance functions
#' 
#' @return list all names of covariance functions recommend
#' as of this edition.  The first is the default.
#' @export
listcov = function() {
  c("mat25pow", "mat25", "mat25ang")
}

#' fit an outerbase
#' 
#' @param x a n by d sized matrix of inputs
#' @param y a n length vector of outputs
#' @param numb size of basis to use
#' @param covnames a d length vector of covariance names
#' @return Saving important model information
#' @export
obfit = function(x, y, numb=100, covnames=NULL) {
  if(dim(x)[1] != length(y)) stop("\n x and y dims do not align")
  if(dim(x)[1] < dim(x)[2]) 
    stop('\n dimension larger than sample size has not been tested')
  if(dim(x)[1] > 10^6) 
    stop('\n sample size should be less than 1000000')
  if(dim(x)[2] > 200) 
    stop('\n dimension should be less than 200')
  if(dim(x)[1] > 10^5) 
    warning('\n sample size is larger than has been tested')
  if(dim(x)[2] > 20) 
    warning('\n more than 20 dimensions has not been tested')
  if(dim(x)[2] == 1) 
    stop('\n dimension must be larger than 1')
  if(dim(x)[2] == 2) 
    stop('\n dimension 2 has not been tested')
  if(numb < 2*dim(x)[2]) 
    stop('\n number of basis functions should be over twice the dimension')
  if(numb > 2000) 
    warning('\nnumber of basis functions is large')
  if(numb > 10000) 
    stop('\n number of basis functions is beyond testing')
  if(numb > dim(x)[1]) 
    warning('\n number of basis functions larger than sample size, not tested')
  
  y_cent = mean(y)
  y_sca = sd(y)
  y = (y - y_cent) / y_sca
  
  d = dim(x)[2]
  
  if(!is.null(covnames)){
    if(length(covnames) != d){
      stop("\n cov names must be same size as columns in x")
    }
  }
  
  if(is.null(covnames)) covnames = rep(listcov()[1],d)
  
  for(k in 1:d) .checkcov(covnames[k], x[,k]) 
  
  om = new(outermod)
  setcovfs(om, covnames)
  setknot(om, .genknotlist(rep(40,d), x)) #40 knot point for each dim
  
  terms = om$selectterms(min(numb, 100*d)) #small number of terms
  
  logpr = new(logpr_gauss, om, terms) #initial parameter
  loglik = new(loglik_gda, om, terms, y, x) #initial parameter
  loglik$dodiag = T # for the first round, go ahead and go the diagonal
  logpdf = new(lpdfvec, logpr, loglik)
  
  BFGS_lpdf(om, logpdf,verbose=0) 
  
  terms = om$selectterms(numb) #get new terms
  bassize = ceiling(pmax(16,
                         pmin(70,
                              2*apply(terms,2,max))))
  setknot(om, .genknotlist(bassize, x))
  loglik_faster = new(loglik_gauss, om, terms, y, x) #initial parameter
  logpdf_faster = new(lpdfvec, logpr, loglik_faster)
  logpdf_faster$domarg = T
  
  for(k in 1:2){
    terms = om$selectterms(numb)
    logpdf_faster$updateterms(terms)
    BFGS_lpdf(om, logpdf_faster,verbose=0) 
  }
  obmodel = list()
  obmodel$y_cent = y_cent
  obmodel$y_sca = y_sca
  obmodel$om = om
  obmodel$predobj = new(predictor,loglik_faster)
  obmodel
}

#' pred from an outerbase
#' 
#' @param obmodel output from \code{\link{obfit}}
#' @param x a new m by d sized matrix of inputs
#' @return A list with \code{mean} and \code{var} at new x
#' @export
obpred = function(obmodel, x){
  obmodel$predobj$update(x)
  out = list()
  out$mean = (obmodel$y_cent+obmodel$y_sca*obmodel$predobj$mean())
  out$var = (obmodel$y_sca^2) * obmodel$predobj$var()
  out
}


.checkcov = function(covname, x) {
  if (!(covname %in% listcov())) {
    stop('\n covariances must be from listcov()')
  }
  covobj = new(get(paste("covf_",covname,sep="")))
  
  if (min(x) <covobj$lowbnd || max(x) > covobj$uppbnd){
    stop(paste('\n x ranges exceed limits of covariance functions \n',
               'the limits are between', covobj$lowbnd, 'and', covobj$uppbnd,
               ' \n try rescaling'))
  }
  
  if (diff(range(x)) < 1/20*(covobj$uppbnd-covobj$lowbnd)){
    stop(paste('\n x are too small for ranges\n',
               'the limits are between', covobj$lowbnd, 'and', covobj$uppbnd,
               ' \n try rescaling'))
  }
}

.genknotlist <- function(bassize, x) {
  knotlist = list()
  d = dim(x)[2]
  for(k in 1:d) 
    knotlist[[k]] = quantile(x[,k],
                             seq(0,1,length=bassize[k])*bassize[k]/
                               (bassize[k]+1)+0.5/(bassize[k]+1))
  knotlist
}