#' list all covariance functions
#' 
#' @return list all names of covariance functions recommend
#' as of this edition.  The first is the default.
#' @export
listcov = function() {
  c("mat25pow", "mat25", "mat25ang")
}

#' Outerbase model fit
#' 
#' This function fits an outerbase model for prediction and hides most of the
#' actual object-oriented aspects of the package.
#' 
#' @param x a n by d sized matrix of inputs
#' @param y a n length vector of outputs
#' @param numb size of basis to use
<<<<<<< HEAD
#' @param hyp initial covariance hyperparameters
=======
#' @param hyp initial covariance hyperparmaeters
>>>>>>> main
#' @param verbose 0-3, how much information on optimization to print to console
#' @param covnames a d length vector of covariance names, ignored if \code{omst}
#' @param numberopts number of optimizations done for hyperparameters, must be
#' larger than 1
#' is provided
#' @return Saving important model information to be used with 
#' \code{\link{obpred}}
#' @export
obfit = function(x, y, numb=100, verbose = 0,
                 covnames=NULL, hyp=NULL,
                 numberopts=2) {
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
    stop('\n number of basis functions should be less than twice the dimension')
  if(numb > 5000) 
    warning('\n number of basis functions is large, might take time to fit.')
  if(numb > 100000) 
    stop('\n number of basis functions is beyond testing')
  if(numb > dim(x)[1]) 
    warning(paste0('\n number of basis functions larger than sample size, \n',
    'this has not been thoroughly tested'))
  
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
  
  if (length(hyp) == length(gethyp(om))) om$updatehyp(hyp)
  
  setknot(om, .genknotlist(rep(40,d), x)) #40 knot point for each dim
  
  numbr = min(c(floor(length(y)/2), numb, 80*d))
  terms = om$selectterms(min(c(numbr))) #small number of terms
  ssr = min(c(length(y),3*numbr))
  logpr = new(logpr_gauss, om, terms) #initial parameter
  subsetinds = sample(length(y),ssr)
  yr = y[subsetinds]
  xr = x[subsetinds,]
  loglik = new(loglik_gda, om, terms, yr, xr) #
  loglik$dodiag = T # for the first round, go ahead and go the diagonal
  logpdf = new(lpdfvec, logpr, loglik)
  
  if (verbose >0) {
    cat('doing partial optimization ', '\n')
    if (verbose >1) cat('max number of cg steps set to', 100, '\n')
  }
  optinfo = BFGS_lpdf(om, logpdf, verbose = verbose,
                      cgsteps=100) 
  
  terms = om$selectterms(numb) #get new terms
  bassize = ceiling(pmax(16,
                         pmin(70,
                              2*apply(terms,2,max))))
  setknot(om, .genknotlist(bassize, x))
  
  loglik_faster = new(loglik_gauss, om, terms, y, x) #initial parameter
  logpdf_faster = new(lpdfvec, logpr, loglik_faster)
  logpdf_faster$domarg = T 
  #one fewer para, so we will strip that one off
  optinfo$B = optinfo$B[,-ncol(optinfo$B)]
  optinfo$B = optinfo$B[-nrow(optinfo$B),]
  #decrease scale
  optinfo$B = length(yr)/length(y)*optinfo$B
  logpdf_faster$updatepara(getpara(logpdf)[1:2])
  
  for(k in 1:numberopts){
    nsteps = .getsteps(numb, length(y),
                       var(y) / exp(2*getpara(logpdf_faster)[2]))
    if (verbose >0) {
      cat('doing optimization ', k, '\n')
      if (verbose >1) cat('max number of cg steps set to', nsteps, '\n')
    }
    terms = om$selectterms(numb)
    logpdf_faster$updateterms(terms)
    optinfo = BFGS_lpdf(om, logpdf_faster, verbose=verbose,
                        B = optinfo$B, lr = optinfo$lr/2,
                        cgsteps=nsteps) 
  }
  obmodel = list()
  obmodel$y_cent = y_cent
  obmodel$y_sca = y_sca
  obmodel$om = om
  obmodel$predobj = new(predictor,loglik_faster)
  obmodel
}

#' Prediction from outerbase
#' 
#' This function allows for turning an \code{obmodel} into predictions with 
#' mean and variance.
#' 
#' @param obmodel output from \code{\link{obfit}}
#' @param x a new m by d sized matrix of inputs
#' @return A list with \code{mean} and \code{var} at new x
#' @seealso obfit
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
  
  if (min(x) < covobj$lowbnd || max(x) > covobj$uppbnd){
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


.getsteps = function(numb, sampsize,
                     sigtonoiseratio=10^(-3), tol = 0.001) {
  kapp = (1+sqrt(numb/sampsize))^2 / (1-sqrt(numb/sampsize))^2 #semi circle law
  kapp = min(1000,kapp) # in case ratio is off
  # cg complexity computes below
  iterest = 1/2 * sqrt(kapp)* log(2*sampsize*sigtonoiseratio/tol) 
  ceiling(1.5*iterest) # 1.5 is safety feature
}
