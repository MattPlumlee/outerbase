#include "customconfig.h"
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp; 

#include "covfuncs.h"
#include "modandbase.h"
#include "fit.h"

//' @name setcovfs
//' @title Set covariance functions
//' @description
//' \preformatted{
//' setcovfs(om, covnames)
//' }
//' This function sets the covariance functions for an outermod object,
//' which is first thing one does when creating an outermod object.//'
//' @param om \code{\link{outermod}} object
//' @param covnames vector of strings of the covariance functions
//' @examples
//' om = new(outermod)
//' setcovfs(om, c("mat25", "mat25", "mat25"))
//' setcovfs(om, c("mat25", "mat25pow", "mat25", "mat25ang"))
//' @seealso \code{\link{outermod}} 
void setcovfs(outermod& om, StringVector covstr){
  
  om.d = covstr.size();
    
  om.covflist.clear();
  for (uword k = 0; k < om.d; ++k) {  
    if (covstr(k)=="mat25"){
      om.covflist.push_back(new covf_mat25());
    } else if (covstr(k)=="mat25pow") {
      om.covflist.push_back(new covf_mat25pow());
    } else if (covstr(k)=="mat25ang") {
      om.covflist.push_back(new covf_mat25ang());
    } else {
    //  throw std::range_error("need to choose one of the existing cov functions");
    }
  }
  
  om.hyp_init();
  om.setcovfs = true;
  om.setknots = false;
}

//' @name setknot
//' @title Set knot points
//' @description
//' \preformatted{
//' setknot(om, knotslist)
//' }
//' This function sets the knot points to estimate the eigenfunctions
//' and eigenvalues. It will naturally check if the knot points have the
//' same dimension as the covariance functions.  It will also check if the 
//' knot points are within reasonable bounds for the covariance functions.//' 
//' @param om \code{\link{outermod}} object
//' @param knotslist list of one dimensional vectors
//' @examples
//' om = new(outermod)
//' setcovfs(om, c("mat25", "mat25", "mat25"))
//' setknot(om,
//'          list(seq(0,1,by=0.01),seq(0,1,by=0.01),seq(0,1,by=0.01)))
//' @seealso \code{\link{outermod}}, \code{\link{setcovfs}}
void setknot(outermod& om, List L){
  /* Need to figure out how to do error checking later
  if (!om.setcovfs) {
    throw std::range_error("Need to set cov. funcs before setting knots.");
    return;
  }
  
  if (L.length() != om.d) {
    std::string errormsgdims = "dim needs to match" ;
    errormsgdims += std::to_string(om.d);
    errormsgdims += ".";
    throw std::range_error(errormsgdims);
    return;
  }
  
  for (unsigned int l = 0; l < om.d; ++l) {
    bool inbnds = (*om.covflist[l]).inputcheck(as<vec>(L[l]));
    if(!inbnds){
      std::string errormsgbnd = std::to_string(l+1);
      errormsgbnd += "knot point needs to be between " ;
      errormsgbnd += std::to_string((*om.covflist[l]).lowbnd);
      errormsgbnd += " and ";
      errormsgbnd += std::to_string((*om.covflist[l]).uppbnd);
      throw std::range_error(errormsgbnd);
      return;
    }
  }
   */
  
  om.knotptst.resize(om.d+1);
  int currst = 0;
  for (unsigned int l = 0; l < om.d; ++l) {
    om.knotptst[l] = currst;
    currst += (as<vec>(L[l])).n_elem;
  }
  om.knotptst[om.d] = currst;
  om.knotpt.resize(currst);
  for (unsigned int l = 0; l < om.d; ++l)
    om.knotpt.subvec(om.knotptst[l],om.knotptst[l+1]-1) = as<vec>(L[l]);
  om.setknots = true;
  
  om.knotptstge.resize(om.d+1);
  om.gest.resize(om.hypst[om.d]+1);
  currst = 0;
  int currstalt = 0;
  for (unsigned int l = 0; l < om.d; ++l){
    om.knotptstge[l] = currst;
    for (unsigned int k = 0; k < (om.hypst[l+1]-om.hypst[l]); ++k){
      om.hypmatch[currstalt] = l;
      om.gest[currstalt] = currst;
      currst += om.knotptst[l+1]-om.knotptst[l];
      currstalt += 1;
    }}
  om.knotptstge[om.d] = currst;
  om.gest[om.hypst[om.d]] = currst;
  
  om.build();
}

//' @name gethyp
//' @title Get the hyperparameters
//' @param om \code{\link{outermod}} object
//' @description
//' \preformatted{
//' hyp = gethyp(om)
//' }
//' This function gets the current hyperparameters from an outermod object. It 
//' formats them in a way that makes reading in `R` helpful. 
//' @examples
//' om = new(outermod)
//' setcovfs(om, c("mat25", "mat25", "mat25"))
//' hyp = gethyp(om)
//' print(hyp)
//' @seealso \code{\link{outermod}}
NumericVector gethyp(outermod& om) {
  NumericVector out(om.hyp.n_elem);
  CharacterVector hypnames(out.length());
  for (uword k = 0; k < out.length(); ++k) {  
    out[k] = om.hyp[k];
    hypnames[k] = "inpt";
    hypnames[k] += std::to_string(1+om.hypmatch[k]);
    hypnames[k] += ".";
    hypnames[k] += (*om.covflist[om.hypmatch[k]])\
      .hypnames[k-om.hypst[om.hypmatch[k]]];
  }
  out.names() = hypnames;
  return out;
} 


//' @name getpara
//' @title Get the model parameters
//' @description
//' \preformatted{
//' para = getpara(loglik)
//' }
//' This function gets the current parameters from an \code{\link{lpdf}} object.  
//' It formats them in a way that makes reading in `R` helpful.
//' @param lpdf logpdf object
NumericVector getpara(lpdf& logpdf) {
  NumericVector out(logpdf.para.n_elem);
  out.names() = logpdf.paranames;
  for (uword k = 0; k < out.length(); ++k)  
    out[k] = logpdf.para[k];
  return out;
} 

RCPP_EXPOSED_CLASS(outerbase)
RCPP_EXPOSED_CLASS(outermod)
RCPP_EXPOSED_CLASS(lpdf)
  
  
//' @name outermod
//' @aliases
//' Rcpp_outermod-class Rcpp_outermod
//' @title Outer product-type model
//' @description Type the name of the class to see its methods.
//' @field \link{setcovfs} to set covariance functions
//' @field \link{setknot} to set knot points
//' @field \link{gethyp} to get hyperparameters
//' @field \link{outermod$updatehyp} to update hyperparameters
//' @field \link{outermod$selectterms} to find best terms
//' @field \link{outermod$getvar} to find variances of coefficients
//' @examples
//' om = new(outermod)
//' setcovfs(om, c("mat25", "mat25", "mat25"))
//' setknot(om,
//'          list(seq(0,1,by=0.01),seq(0,1,by=0.01),seq(0,1,by=0.01)))
//' terms = om$selectterms(40)
//' coeffvar =om$getvar(terms)
//' hyp = gethyp(om)
//' hyp[1:2] = 0.5
//' om$updatehyp(hyp)
//' coeffvar = om$getvar(terms)
//' @seealso \code{\link{outerbase}} the main product from outermods

//' @name outermod$updatehyp
//' @title Update hyperparameters
//' @description
//' \preformatted{
//' outermod$updatehyp(hyp)
//' }
//' Updates the hyperparameters
//' @param hyp A vector of hyperparameters
//' @seealso \code{\link{outermod}}

//' @name outermod$selectterms
//' @title Select optimal terms
//' @description
//' \preformatted{
//' terms = om$selectterms(numterms)
//' }
//' Selects the best \code{terms} given the current \code{outermod}
//' @param numterms Number of basis \code{terms} desired
//' @returns terms A matrix of \code{terms}
//' @seealso \code{\link{outermod}}


//' @name outermod$getvar
//' @title Get variance of coefficients
//' @description
//' \preformatted{
//' coeffvar = outermod$getvar(terms)
//' }
//' Returns the variance of the coefficients associated with \code{terms}.
//' @param terms A matrix of \code{terms}
//' @returns coeffvar A vector of variances of each coefficient
//' @seealso \code{\link{outermod}}


//' @name outerbase
//' @aliases 
//' Rcpp_outerbase-class Rcpp_outerbase
//' @title Outer product-type basis
//' @description 
//' \preformatted{
//' ob = new(outerbase, om, x)
//' }
//' Object that handles the basis for a given set of points 
//' \code{x}.
//' @param x a matrix of predictors, must have as many columns as dims in 
//' \code{om}
//' @field nthreads number of threads for \code{omp} to use
//' @field \link{outerbase$getbase}(k) to get each dimensions basis 
//' functions
//' @field \link{outerbase$getmat}(terms) to get the basis matrix at 
//' \code{terms}
//' @field \link{outerbase$build}() to (re)build the basis object
//' @field \link{outerbase$matmul}(terms,a) matrix multiply without 
//'building the basis matrix
//' @field \link{outerbase$tmatmul}(terms,a) transpose matrix multiply 
//' without building the basis matrix
//' @examples
//' om = new(outermod)
//' setcovfs(om, c("mat25", "mat25", "mat25"))
//' setknot(om,
//'          list(seq(0,1,by=0.025),seq(0,1,by=0.025),seq(0,1,by=0.025)))
//' x = matrix(runif(10*3),ncol=3)
//' ob = new(outerbase, om, x)
//' terms = om$selectterms(40)
//' basismat = ob$getmat(terms)
//' @seealso \code{\link{outermod}} the core element that controls outerbase

//' @name outerbase$getbase
//' @title Get base functions
//' @description
//' \preformatted{
//' basis_func = outerbase$getbase(k)
//' }
//' Returns the basis for a dimension
//' @param k An integer from that corresponds to the dimension.
//' @returns basis_func A matrix of evaluated basis functions for that 
//' dimension.  Designed mostly for visualization.
//' @seealso \code{\link{outerbase}}

//' @name outerbase$getmat
//' @title Get basis matrix
//' @description
//' \preformatted{
//' basismat = outerbase$getmat(terms)
//'  }
//' Returns the basis matrix
//' @param terms A matrix of terms
//' @returns basismat A matrix of evaluated basis functions based on 
//' \code{terms}.
//' @seealso \code{\link{outerbase}}

//' @name outerbase$build
//' @title Builds the outerbase
//' @description
//' \preformatted{
//' outerbase$build()
//'  }
//' Build (or re-build) a basis based on the recent evaluation 
//' of \code{\link{outermod}}.
//' @seealso \code{\link{outerbase}}

//' @name outerbase$matmul
//' @title Matrix multiply
//' @description
//' \preformatted{
//'  b = outerbase$matmul(terms, a)
//' }
//' Multiplies the basis times a vector without building the basis 
//' matrix.
//' @param terms A matrix of \code{terms}
//' @param a A vector of length the same as the rows in \code{terms}
//' @returns b A vector resulting from the matrix multiplication
//' @seealso \code{\link{outerbase}}

//' @name outerbase$tmatmul
//' @title Transpose Matrix multiply
//' @description
//' \preformatted{
//'   b = outerbase$tmatmul(terms, a)
//' }
//' Multiplies the transpose of the basis times a vector without 
//'  building the basis matrix.
//' @param terms A matrix of \code{terms}
//' @param a A vector of length the same as the rows in \code{outerbase}
//' @returns b A vector resulting from the matrix multiplication
//' @seealso \code{\link{outerbase}}

//' @name lpdf
//' @aliases 
//' Rcpp_lpdf-class Rcpp_lpdf
//' @title Log probability density function class
//' @description This is a base class designed to handle the learning of 
//' the underlying coefficients, hyperparameters, and parameters associated with
//' a specific learning instance.  Polymorphism allows for the implied methods to 
//' be used across several similar classes.
//'
//' @field lpdf$val current value
//' @field lpdf$para current model parameters
//' @field lpdf$coeff current coefficients
//' @field lpdf$compute_val on calling \code{update}, compute value and store in 
//' \code{val}
//' @field lpdf$grad current gradient with respect to coefficients
//' @field lpdf$gradhyp current gradient with respect to covariance hyperparameters
//' @field lpdf$gradpara current gradient with respect to model parameters
//' @field lpdf$compute_grad on calling \code{update}, compute gradient with 
//' respect to coefficients and store in \code{grad}
//' @field lpdf$compute_gradhyp on calling \code{update}, compute gradient
//' with respect to covariance hyparameters and store in \code{gradhyp}
//' @field lpdf$compute_gradpara on calling \code{update}, compute gradient
//' with respect to model parameters and store in \code{gradpara}
//' @field lpdf$update(coeff) update using new coefficients
//' @field \link{lpdf$optcg}(tol,epoch) do optimization with respect to coefficients 
//' via conjugate gradient
//' @field \link{lpdf$optnewton}() do optimization via matrix inversion, one Newton 
//' step
//' @field lpdf$updateom() update based on recent version of \code{\link{outermod}}
//' @field lpdf$updatepara(para) update using new model parameters
//' @field lpdf$updateterms(terms) update using new \code{terms}
//' @field lpdf$hess() returns the hessian with respect 
//' to coefficients
//' @field lpdf$hessgradhyp() returns gradient of \code{hess()} with respect to 
//' covariance hyperparameters
//' @field lpdf$hessgradpara() returns the gradient of \code{hess()} with respect to 
//' model parameters
//' @field lpdf$diaghess() returns the diagonal of the hessian with 
//' respect to coefficients
//' @field lpdf$diaghessgradhyp() returns the gradient of \code{diaghess()} with 
///' respect to  covariance hyperparameters
//' @field lpdf$diaghessgradpara() returns the gradient of \code{diaghess()} with 
//' respect to model parameters
//' @field lpdf$paralpdf(para) compute the log-prior on the parameters, useful for 
//' fitting
//' @field lpdf$paralpdf_grad(para) gradient of \code{paralpdf(para)}
//' @seealso container class: \code{\link{lpdfvec}}
//' @seealso derived classes: \code{\link{loglik_std}}, 
//' \code{\link{loglik_gauss}}, \code{\link{loglik_lda}}, 
//' \code{\link{logpr_gauss}}

//' @name lpdfvec
//' @aliases
//' Rcpp_lpdfvec-class Rcpp_lpdfvec 
//' @title Vector of lpdf objects
//' @description 
//' \preformatted{
//' logpdf = new(lpdfvec, loglik, logpr)
//' }
//' This is a class that contains two \code{\link{lpdf}} object and can be 
//' manipulated as a single object.  It presumes both are based on the same
//' \code{\link{outermod}} object, thus they share hyperparameters.  However
//' the model parameters are concatenated.  Currently also includes variations
//' on marginal adjustments.  
//'
//' Currently it is designed only for a pair, but the ordering is arbitrary.
//'
//' @param loglik one reference to a \code{lpdf} object
//' @param logpr another reference to a \code{lpdf} object that shares 
//' \code{\link{outermod}}  with \code{loglik}
//' @field lpdfvec$domarg  A boolean that controls if marginal adjustment is 
//' done
//' @seealso base class: \code{\link{lpdf}}

//' @name lpdf$optcg
//' @title Optimization via Conjugate Gradient
//' @description 
//' \preformatted{
//' lpdf$optcg(tol,epoch)
//' }
//' This optimizes the coefficient vector \code{coeff} using conjugate gradient. 
//' It currently is designed only for quadratic \code{\link{lpdf}} objects.
//' @param tol A positive double representing tolerance, default is 
//' \code{0.001}.
//' @param epoch A positive integer representing the maximum number of steps 
//' conjugate gradient will take.
//' @seealso \code{\link{lpdf}}

//' @name lpdf$optnewton
//' @title Optimization via Newton's Method
//' @description 
//' \preformatted{
//' lpdf$optnewton()
//' }
//' This optimizes the coefficient vector \code{coeff} using Newton's Method.  
//' It currently is designed only for quadratic \code{\link{lpdf}} objects.  
//' It should take a single step.
//' @seealso \code{\link{lpdf}}  

//' @name loglik_std
//' @aliases 
//' Rcpp_loglik_std-class Rcpp_loglik_std 
//' @title Gaussian errors
//' @description 
//' \preformatted{
//' loglik = new(loglik_std, om, terms, y, x)
//' }
//' This is a standard model which has the form
//' \deqn{y = \langle \phi(x), \theta \rangle + \varepsilon, \varepsilon \sim 
//' N(0,\sigma^2)}
//' where \eqn{\phi(x)} is the basis, \eqn{\theta} is the coefficient vector,
//' \eqn{\varepsilon} is an unseen noise vector.
//' The parameter vector is of length 1 where 
//' \code{para} \eqn{= \log(\sigma)}.  It is a slightly slower (sometimes) 
//' version of \code{\link{loglik_gauss}} but allows for complete marginal 
//' inference.
//' @param om an \code{\link{outermod}} object to be referred to
//' @param terms a matrix of \code{terms}, must have as many columns as dims in 
//' \code{om}
//' @param y a vector of observations
//' @param x a matrix of predictors, must have as many columns as dims in 
//' \code{om} and the same number of rows as \code{y}
//' @inherit lpdf description
//' @seealso base class: \code{\link{lpdf}}

//' @name loglik_gauss
//' @aliases 
//' Rcpp_loglik_gauss-class Rcpp_loglik_gauss 
//' @title Gaussian errors, large scale
//' @description 
//' \preformatted{
//' loglik = new(loglik_gauss, om, terms, y, x)
//' }
//' This is a standard model which has the form
//' \deqn{y = \langle \phi(x), \theta \rangle + \varepsilon, \varepsilon \sim 
//' N(0,\sigma^2)}
//' where \eqn{\phi(x)} is the basis, \eqn{\theta} is the coefficient vector,
//' \eqn{\varepsilon} is an unseen noise vector. 
//' The parameter vector is of length 1 where 
//' \code{para} \eqn{= \log(\sigma)}.  It is a slightly (sometimes) version of
//' \code{\link{loglik_std}}  but allows can only handle diagonal variational 
//' inference.
//' @param om an \code{\link{outermod}} object to be referred to
//' @param terms a matrix of \code{terms}, must have as many columns as dims in 
//' \code{om}
//' @param y a vector of observations
//' @param x a matrix of predictors, must have as many columns as dims in 
//' \code{om} and the same number of rows as \code{y}
//' @inherit lpdf description
//' @seealso base class: \code{\link{lpdf}}

//' @name loglik_gda
//' @aliases
//' Rcpp_loglik_gda-class Rcpp_loglik_gda 
//' @title Gaussian errors with diagonal adjustment
//' @description 
//' \preformatted{
//' loglik = new(loglik_gda, om, terms, y, x)
//' }
//' This is a standard model which has the form
//' \deqn{y = \langle \phi(x), \theta \rangle + \delta(x) + \varepsilon,
//' \delta(x) \sim N(0, \lambda g(x)), \varepsilon \sim N(0,\sigma^2)}
//' where \eqn{\phi(x)} is the basis, \eqn{\theta} is the coefficient vector,
//' \eqn{\delta(x)} is unseen vector corresponding to unmodeled 
//' variance \eqn{\lambda g(x)}, \eqn{\varepsilon} is an unseen noise vector.
//' The parameter vector is of length 2 where 
//' \eqn{\sigma=} \code{exp(para[0])} and \eqn{\lambda=}\code{exp(2 para[1])}.
//' @param om an \code{\link{outermod}} object to be referred to
//' @param terms a matrix of \code{terms}, must have as many columns as dims in 
//' \code{om}
//' @param y a vector of observations
//' @param x a matrix of predictors, must have as many columns as dims in 
//' \code{om} and the same number of rows as \code{y}
//' @inherit lpdf description
//' @seealso base class: \code{\link{lpdf}}

//' @name logpr_gauss
//' @aliases 
//' Rcpp_logpr_gauss-class Rcpp_logpr_gauss
//' @title Gaussian prior
//' @description 
//' \preformatted{
//' logpr = new(logpr_gauss, om, terms)
//' }
//' This is a standard model of coefficents which has them as drawn independently
//' from
//' \deqn{ \theta_i \sim N(0, \rho c_i)}
//' where \eqn{c_i} is the variance supplied by \code{om} for the $i$th term. 
//' The parameter vector is of length 1 where 
//' \eqn{\rho=} \code{exp(para[0])}.
//' @param om an \code{\link{outermod}} object to be referred to
//' @param terms a matrix of \code{terms}, must have as many columns as dims in 
//' \code{om}
//' @inherit lpdf description
//' @seealso base class: \code{\link{lpdf}}

//' @name covf
//' @aliases 
//' Rcpp_covf-class Rcpp_covf 
//' @title covariance function class
//' @description This is a base class designed to handle the specific features of 
//' covariances needed for outerbase.  Polymorphism allows for the implied 
//' methods to be used across several similar classes.
//' @field covf$hyp hyperparameters for this specific correlation function
//' @field covf$lowbnd,covf$uppbnd upper and lower bounds for the inputs to the 
//' covariance function.
//' @field covf$cov(x1,x2) returns the covariance matrix between two vectors of 
//' inputs \code{x1} and \code{x2}
//' @field covf$covdiag(x1) returns the diagonal of the covariance matrix between 
//' \code{x1} and itself
//' @field covf$cov_gradhyp(x1,x2) returns a cube of the gradient the \code{cov} 
//' with respect to the covariance hyperparameters
//' @seealso derived class: \code{\link{covf_mat25}}, 
//' \code{\link{covf_mat25pow}}, 
//' \code{\link{covf_mat25ang}}

//' @name covf_mat25
//' @aliases 
//' Rcpp_covf_mat25-class Rcpp_covf_mat25 
//' @title Matern covariance function
//' @description 
//' \preformatted{
//' covf = new(covf_mat25)
//' }
//' This is the standard Matern covariance function which has form
//' \deqn{c(x_1,x_2) = (1+ |h| + h^2/3) \exp(-|h|) }
//' where \eqn{h = (x_1-x_2)/\rho} and \eqn{\rho}=\code{exp(2*hyp[0])}.
//' @seealso base class: \code{\link{covf}}

//' @name covf_mat25pow
//' @aliases 
//' Rcpp_covf_mat25pow-class Rcpp_covf_mat25pow 
//' @title Matern covariance function with power transform
//' @description 
//' \preformatted{
//' covf = new(covf_mat25pow)
//' }
//' This is the standard Matern covariance function with a power transformation
//' which has form
//' \deqn{c(x_1,x_2) = (1+ |h| + h^2/3) \exp(-|h|) }
//' where \eqn{h = (x_1^\alpha-x_2^\alpha)/\rho} and \code{hyp} is a two 
//' dimensional vector with \eqn{\rho}=\code{exp(2*hyp[0]+0.25*hyp[1])}
//' and \eqn{\alpha}=\code{exp(2*hyp[0])}.
//' @seealso base class: \code{\link{covf}}

//' @name covf_mat25ang
//' @aliases 
//' Rcpp_covf_mat25ang-class Rcpp_covf_mat25ang 
//' @title Matern covariance function with angular transform
//' @description 
//' \preformatted{
//' covf = new(covf_mat25ang)
//' }
//' This is the standard Matern covariance function with a power transformation
//' which has form
//' \deqn{c(x_1,x_2) = (1+ |h| + h^2/3) \exp(-|h|) }
//' where \deqn{h = (\sin(x_1)-\sin(x_2))^2/\rho_s + 
//' (\cos(x_1)-\cos(x_2))^2/\rho_c,.}
//' \code{hyp} is a two dimensional vector with 
//' \eqn{\rho_s}=\code{exp(2*hyp[0])} and \eqn{\rho_c}=\code{exp(2*hyp[1])}.
//' @seealso base class: \code{\link{covf}}


//' @name predictor
//' @aliases 
//' Rcpp_predictor-class Rcpp_predictor 
//' @title prediction class
//' @description 
//' \preformatted{
//' pred = new(predictor, logpdf)
//' }
//' This is a base class design to allow for coherent building of
//' predictions across multiple models.  Unlike many base classes in this 
//' package, it is meant to be directly used.
//'
//' @param logpdf An \code{\link{lpdf}} instance to build the predictor
//' @field predictor$update(x) update the current input for prediction
//' @field predictor$mean() return the current mean of the prediction
//' @field predictor$var() return the current var of the prediction

RCPP_MODULE(obmod){
  using namespace Rcpp;
  using namespace arma;
  
  function("setcovfs", &setcovfs, "type ?setcovfs");
  function("setknot", &setknot, "type ?setknot");
  function("gethyp", &gethyp, "type ?gethyp");
  function("getpara", &getpara,"type ?getpara");
  
  class_<outermod>("outermod")
    .constructor()
    .method("updatehyp",&outermod::hyp_set)
    .method("selectterms",&outermod::selectterms)
    .method("getvar",&outermod::getvar)
    .method("getlvar_gradhyp",&outermod::getlvar_gradhyp)
    .method("hyplpdf",&outermod::hyplpdf)
    .method("hyplpdf_grad",&outermod::hyplpdf_grad)
  ;
  
  class_<outerbase>("outerbase")
    .constructor<const outermod &, mat>()
    .field("nthreads", &outerbase::nthreads)
    .method("getbase",&outerbase::getbase)
    .method("getmat",&outerbase::getmat)
    .method("build",&outerbase::build)
    .method("matmul",&outerbase::mm_out)
    .method("tmatmul",&outerbase::tmm_out)
    .method("getmat_gradhyp",&outerbase::getmat_gradhyp)
    .method("matmul_gradhyp",&outerbase::mm_gradhyp_out)
    .method("tmatmul_gradhyp",&outerbase::tmm_gradhyp_out)
  ;
  
  class_<lpdf>("lpdf")  
    .constructor()
    .field("compute_val", &lpdf::compute_val)
    .field("compute_grad", &lpdf::compute_grad)
    .field("compute_gradpara", &lpdf::compute_gradhyp)
    .field("compute_gradhyp", &lpdf::compute_gradpara)
    .field_readonly("fullhess", &lpdf::fullhess)
    .field_readonly("val", &lpdf::val)
    .field_readonly("coeff", &lpdf::coeff)
    .field_readonly("grad", &lpdf::grad)
    .field_readonly("gradhyp", &lpdf::gradhyp)
    .field_readonly("gradpara", &lpdf::gradpara)
    .field_readonly("para", &lpdf::para)
    .field_readonly("nterms", &lpdf::nterms)
    .method("setnthreads", &lpdf::setnthreads)
    .method("optcg", &lpdf::optcg)
    .method("optnewton", &lpdf::optnewton)
    .method("update", &lpdf::update)
    .method("updateom", &lpdf::updateom)
    .method("updatepara", &lpdf::updatepara)
    .method("updateterms", &lpdf::updateterms)
    .method("hessmult", &lpdf::hessmult)
    .method("diaghess", &lpdf::diaghess)
    .method("diaghessgradhyp", &lpdf::diaghessgradhyp)
    .method("diaghessgradpara", &lpdf::diaghessgradpara)
    .method("paralpdf", &lpdf::paralpdf)
    .method("paralpdf_grad", &lpdf::paralpdf_grad)
  ;
  //.method("hess", &lpdf::hess)
  //.method("hessgradhyp", &lpdf::hessgradhyp)
  //.method("hessgradpara", &lpdf::hessgradpara)
  //.field_readonly("totdiaghess", &lpdf::totdiaghess)
  //.field_readonly("tothess", &lpdf::tothess)
  
  class_<predictor>("predictor")
    .constructor<const lpdf&>()
    .method("update",&predictor::update)
    .method("mean",&predictor::mean)
    .method("var",&predictor::var)
  ;
  
  class_<loglik_std>("loglik_std")
    .derives<lpdf>("lpdf")
    .constructor<const outermod&, umat, vec, mat>()
    .field_readonly("yhat",&loglik_std::yhat)
  ;
  
  class_<loglik_gauss>("loglik_gauss")
    .derives<lpdf>("lpdf")
    .constructor<const outermod&, umat, vec, mat>()
  ;
  
  class_<loglik_gda>("loglik_gda")
    .derives<lpdf>("lpdf")
    .constructor<const outermod&, umat, vec, mat>()
    .field("dodiag", &loglik_gda::doda)
  ;
  
  class_<logpr_gauss>("logpr_gauss")
    .derives<lpdf>("lpdf")
    .constructor<const outermod&, umat>()
    .field_readonly("thetasd", &logpr_gauss::coeffsd)
  ;
  
  class_<lpdfvec>("lpdfvec")
    .derives<lpdf>("lpdf")
    .constructor<lpdf&, lpdf&>()
    .field("domarg", &lpdfvec::domargadj)
  ;
  /*
   * 
   .field("val_margadj", &lpdfvec::val_margadj)
   .field("gradhyp_margadj", &lpdfvec::gradhyp_margadj)
   .field("gradpara_margadj", &lpdfvec::gradpara_margadj)
   */
  
  class_<covf>("covf")
    .constructor()
    .field("hyp", &covf::hyp)
    .field_readonly("hyplb", &covf::hyplb)
    .field_readonly("hypub", &covf::hypub)
    .field_readonly("hyp0", &covf::hyp0)
    .field_readonly("hypvar", &covf::hypvar)
    .field_readonly("lowbnd", &covf::lowbnd)
    .field_readonly("uppbnd", &covf::uppbnd)
    .method("cov", &covf::cov)
    .method("covdiag", &covf::covmdiag)
    .method("cov_gradhyp", &covf::cov_gradhyp)
  ;
  
  class_<covf_mat25>("covf_mat25")
    .derives<covf>("covf")
    .constructor()
  ;
  
  class_<covf_mat25pow>("covf_mat25pow")
    .derives<covf>("covf")
    .constructor()
  ;
  
  class_<covf_mat25ang>("covf_mat25ang")
    .derives<covf>("covf")
    .constructor()
  ;
  
}

