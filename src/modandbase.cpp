/**************************************************************************** 
 * 
 * outerbase 
 * Copyright 2022, Matthew Plumlee 
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy 
 * of this software and associated documentation files (the "Software"), to deal 
 * in the Software without restriction, including without limitation the rights 
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
 * copies of the Software, and to permit persons to whom the Software is 
 * furnished to do so, subject to the following conditions: 
 * 
 * The above copyright notice and this permission notice shall be included in 
 * all copies or substantial portions of the Software. 
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN 
 * THE SOFTWARE. 
 * 
 * Questions? Contact Matthew Plumlee 
 *  
 ****************************************************************************/ 

#include "customconfig.h"  
#include <RcppArmadillo.h>  

// JRR on SO suggested this
#ifdef _OPENMP
  #include <omp.h>
#else
  // for machines with compilers void of openmp support
  #define omp_get_num_threads()  1
  #define omp_get_thread_num()   0
  #define omp_get_max_threads()  1
  #define omp_get_thread_limit() 1
  #define omp_get_num_procs()    1
  #define omp_set_nested(a)   // empty statement to remove the call
  #define omp_get_wtime()        0
  #define omp_in_parallel()      true
#endif

using namespace arma;  
using namespace Rcpp;  

#include "covfuncs.h"
#include "modandbase.h"  
#include "linalg.h"  

/* 
 **************************************************************************** 
 **************************************************************************** 
 *************************** outermod *************************************** 
 **************************************************************************** 
 **************************************************************************** 
 */ 

/* 
 * outermod::setsize 
 * 
 * set the sizes of the outermod artifacts 
 */ 

void outermod::setsizes_() 
{ 
  uword indm = diff(knotptst).index_max();  
  uword mmax = knotptst[indm+1]-knotptst[indm];  
  if (rotmat.n_rows != mmax || rotmat.n_cols !=knotpt.n_elem) 
    rotmat.set_size(mmax,knotpt.n_elem);  
  if (rotmat_gradhyp.n_rows != mmax || rotmat_gradhyp.n_cols != knotptstge[d]) 
    rotmat_gradhyp.set_size(mmax,knotptstge[d]);  
  if (basisvar.n_elem != knotpt.n_elem)  
    basisvar.set_size(knotpt.n_elem);  
  if (logbasisvar_gradhyp.n_elem != knotptstge[d]) 
    logbasisvar_gradhyp.set_size(knotptstge[d]);  
  if (maxlevel.n_elem != knotptstge[d]) 
    maxlevel.set_size(d);  
} 

/* 
 * outermod::hypcheck
 * 
 * check if the hyperparameter vector is reasonable
 */ 

double outermod::hyplpdf(vec hypp) const { 
  double out = 0;
  
  if(hyp.n_elem != hypp.n_elem) out = -datum::inf;
  else {
  for (unsigned int l = 0; l <  d; ++l) 
    out += (*covflist[l]).lpdf(hypp.subvec(hypst[l], hypst[l+1]-1)); 
  }
  
  return out;
} 

/* 
 * outermod::hypcheck_gradhyp
 * 
 * check if the hyperparameter vector is reasonable
 */ 

vec outermod::hyplpdf_grad(vec hypp) const { 
  vec out;
  out.set_size(hyp.n_elem);
  out.zeros();
  
  if(hyp.n_elem == hypp.n_elem){
    for (unsigned int l = 0; l <  d; ++l) 
      out.subvec(hypst[l], hypst[l+1]-1) = 
        (*covflist[l]).lpdf_gradhyp(hypp.subvec(hypst[l], hypst[l+1]-1)); 
  }
  
  return out;
} 


/* 
 * outermod::hypcheck
 * 
 * check if the hyperparameter vector is reasonable
 */ 

void outermod::hyp_init() { 
  
  hypst.resize(d+1);
  unsigned int currst = 0;
  for (unsigned int l = 0; l < d; ++l){
    hypst[l] = currst;
    currst += (*covflist[l]).numhyp;
  }
  hypst[d] = currst;
  
  hyp.set_size(currst);
  for (unsigned int l = 0; l < d; ++l)
    hyp.subvec(hypst[l],hypst[l+1]-1) = (*covflist[l]).hyp0;
  
  hypmatch.resize(hypst[d]);
  int currstalt = 0;
  for (unsigned int l = 0; l < d; ++l){
    for (unsigned int k = 0; k < (hypst[l+1]-hypst[l]); ++k){
      hypmatch[currstalt] = l;
      currstalt++;
    }}
  
  outermod::hyp_set(hyp);
  
  return;
} 

/* 
 * outermod::hypcheck
 * 
 * check if the hyperparameter vector is reasonable
 */ 

void outermod::hyp_set(vec hyp_) { 
  
  hypst.resize(d+1);
  unsigned int currst = 0;
  for (unsigned int l = 0; l < d; ++l){
    hypst[l] = currst;
    currst += (*covflist[l]).numhyp;
  }
  hypst[d] = currst;
  
  /*
  if (currst != hyp_.n_elem) {
    throw std::range_error("wrongsized vector");
    return;
  }
  */
  hyp = hyp_;
  
  for (unsigned int l = 0; l <  d; ++l)
    (*covflist[l]).hyp = hyp.subvec(hypst[l],hypst[l+1]-1);
  
  if (setknots) {
  knotptstge.resize(d+1);
  hypmatch.resize(hypst[d]);
  gest.resize(hypst[d]+1);
  currst = 0;
  int currstalt = 0;
  for (unsigned int l = 0; l < d; ++l){
    knotptstge[l] = currst;
    for (unsigned int k = 0; k < (hypst[l+1]-hypst[l]); ++k){
      hypmatch[currstalt] = l;
      gest[currstalt] = currst;
      currst += knotptst[l+1]-knotptst[l];
      currstalt += 1;
    }}
  knotptstge[d] = currst;
  gest[hypst[d]] = currst;
  
  build();
  }
  return;
} 

/* 
 * outermod::build 
 * 
 * build a model for outerbase based on a Gaussian process model 
 */ 

void outermod::build() 
{ 
  /*
  if (!setknots || !setcovfs) {
    throw std::range_error("Need to set covfs and knots before building.");
    return;
  }
   */
  
  outermod::setsizes_(); 
  
  // make sure dims agree 
  // assign a global correlation function for now 
  
  // building all matrices to induce approximate orthgonality 
  for (uword k = 0; k < d; ++k) {  
    
    // get x's for this specific dim  
    vec xsh = knotpt.subvec(knotptst[k],knotptst[k+1]-1);  
    uword lenh = xsh.n_elem; 
    vec hyp_ = hyp.subvec(hypst[k],hypst[k+1]-1); 
    // compute covariance 
    mat R = (*covflist[k]).cov(xsh, xsh); 
    // Do eigenvalue decomposition 
    vec sr(lenh);  
    mat U(lenh,lenh);  
    eig_sym(sr, U, R);  
    sr = reverse(sr); //bigger should always be first... 
    U = fliplr(U); 
    
    /* Removing sign ambiguity */ 
    uword halfw = lenh/2; 
    U.each_row() %= sign(U.row(halfw) + 2.71828*U.row(halfw+1)); 
    
    /* Hardcoded fixes to insure stable linear algebra */ 
    double minsv = 0.00000000001*mean(sr); 
    uvec issmall = find(-diff(sr) < minsv, 1);  
    if (issmall.n_elem > 0) maxlevel(k) = issmall.min();  
    else maxlevel(k) = sr.n_elem-1; 
    sr = sr+linspace(minsv/1000, lenh*minsv/1000, lenh); 
    
    /* Save for later usage */ 
    rotmat.submat(0,knotptst[k],lenh-1,knotptst[k]+lenh-1) = U;  
    rotmat.submat(0,knotptst[k],lenh-1,knotptst[k]+lenh-1).each_row() /= 
      (sr / sqrt(lenh)).as_row(); // making sure the vector is stable 
    basisvar.subvec(knotptst[k],knotptst[k]+lenh-1) = log(sr/lenh);  
    
    /* Compute gradient matrices */ 
    cube Rge = (*covflist[k]).cov_gradhyp(xsh,xsh); 
    mat Fm(lenh,lenh);  
    Fm.each_row() = sr.as_row();  
    Fm.diag().zeros(); 
    Fm.each_col() -= sr;  
    Fm = 1/Fm; 
    mat Ah(lenh,lenh);  
    mat UtdRV(R.n_rows,lenh); 
    for (uword l = 0; l < (hypst[k+1]-hypst[k]); ++l) {  
      UtdRV = U.t() * Rge.slice(l) * U;  
      logbasisvar_gradhyp.subvec(knotptstge[k]+l*lenh, 
                                 knotptstge[k]+(l+1)*lenh-1) = diagvec(UtdRV)/sr;  
      Ah = U * (UtdRV % Fm);  
      Ah.each_row() /= (sr / sqrt(lenh)).as_row(); 
      rotmat_gradhyp.submat(0,knotptstge[k]+l*lenh, 
                            lenh-1,knotptstge[k]+(l+1)*lenh-1) = Ah;  
    }  
  }  
} 


/* 
 * outermod::buildob 
 * 
 * build a outerbase class from an outermod object 
 */ 

void outermod::buildob(mat& R, const mat& xp, const uword& k) const{ 
  uword lenh = knotptst[k+1]-knotptst[k]; 
  
  //mat (*covftemp)(const vec&, const vec&, const vec&);
  //cube (*covftemp_gradhyp)(const vec&, const vec&, const vec&);
  
  //covftemp = covfs[k];
  R = (*covflist[k]).cov(xp.col(k),  
               knotpt.subvec(knotptst[k],knotptst[k+1]-1)); 
  R *= rotmat.submat(0,knotptst[k],lenh-1,knotptst[k+1]-1);  
  
  /* Standardizing for our linear algebra */ 
  R.cols(1,lenh-1).each_col() /= R.col(0); //all divided by first column 
} 

/* 
 * outermod::buildob 
 * 
 * build a outerbase class from an outermod object 
 */ 

void outermod::buildob(mat& R, cube& Rt, 
                       const mat& xp, const uword& k) const{ 
  uword lenh = knotptst[k+1]-knotptst[k]; 
  
  R = (*covflist[k]).cov(xp.col(k),  
               knotpt.subvec(knotptst[k],knotptst[k+1]-1)); 
  Rt =(*covflist[k]).cov_gradhyp(xp.col(k),  
                        knotpt.subvec(knotptst[k],knotptst[k+1]-1)); 
  /* Orthogonalizing with rotmat */ 
  for (uword l = hypst[k]; l <hypst[k+1]; ++l) {   
    Rt.slice(l-hypst[k]) *=\
      rotmat.submat(0, knotptst[k],lenh-1,knotptst[k+1]-1);   
    Rt.slice(l-hypst[k]) += R *\
      rotmat_gradhyp.submat(0, gest[l],lenh-1,gest[l+1]-1);   
  }   
  R *= rotmat.submat(0,knotptst[k],lenh-1,knotptst[k+1]-1);  
  
  /* Standardizing for our linear algebra */ 
  for (uword l = hypst[k]; l < hypst[k+1]; ++l)  
    Rt.slice(l-hypst[k]).each_col() /=  R.col(0); //all divided by first column 
  R.cols(1,lenh-1).each_col() /= R.col(0); //all divided by first column 
} 


/* 
 * outermod::buildvar 
 * 
 * description 
 */ 

vec outermod::totvar(const mat& xp) const{ 
  vec out(xp.n_rows); 
  out.ones(); 
  for (uword k = 0; k < d; ++k) 
    out %= (*covflist[k]).covmdiag(xp.col(k)); 
  return out; 
} 

/* 
 * outermod::getvar 
 * 
 * given a set of terms, return the variance using basisvar 
 */ 

vec outermod::getvar(const umat& terms) const {   
  vec out(terms.n_rows);   
  for (uword k = 0; k < terms.n_rows; ++k) {   
    out(k) = sum(basisvar.elem(knotptst.head(d).as_row()+terms.row(k)));   
  }   
  return exp(out);   
}   

/* 
 * outermod::getlvar_gradhyp 
 * 
 * given a set of terms, return the variance using basisvar 
 */ 

mat outermod::getlvar_gradhyp(const umat& terms) const {   
  vec outa(terms.n_rows);   
  for (uword k = 0; k < terms.n_rows; ++k) {   
    outa(k) = sum(basisvar.elem(knotptst.head(d).as_row()+terms.row(k)));   
  }   
  
  vec outh(terms.n_rows);   
  mat out(terms.n_rows, hypmatch.n_elem);   
  out.zeros();   
  
  for (uword k = 0; k < d; ++k) {   
    for (uword l = hypst[k]; l < hypst[k+1]; ++l) {   
      out.col(l) += logbasisvar_gradhyp.elem(gest(l)+terms.col(k));   
    }}   
  return out;   
}   

/* 
 * outermod::selectterms 
 * 
 * given a current outermod, choose the most value terms 
 */ 

umat outermod::selectterms(const unsigned int numele) const { 
  imat terms(numele,d); //preallocate terms 
  terms.zeros(); 
  vec tv(numele); //value of these terms 
  tv.zeros(); 
  unsigned int nd = 0; //number of current chosen terms 
  
  imat pterms(10*numele,d); //preallocate potential terms 
  pterms.zeros(); 
  vec ptv(10*numele);//value of these terms 
  ptv.zeros(); 
  
  uvec indsh(d); //getting starting with a value of 0s 
  indsh = conv_to<uvec>::from(knotptst.head(d).as_row()+pterms.row(0));  
  ptv(0) = sum(basisvar.elem(indsh));  
  unsigned int np = 1; //number of current potential terms 
  
  unsigned int kstar = 0; //keeping which one we point to 
  for (unsigned int k = 0; k < numele; ++k) { 
    double mval = -0.1+max(ptv.head(np)); 
    uvec islarge = find(ptv.head(np) > mval);  
    islarge = shuffle(islarge); 
    kstar = islarge(0);// choose the maximum 
    terms.row(nd) = pterms.row(kstar); 
    tv(nd) = ptv(kstar); 
    nd += 1; 
    np -= 1; 
    if(np > kstar){ //move things back  
      pterms.row(kstar) = pterms.row(np); 
      ptv.row(kstar) = ptv.row(np); 
    } 
    
    /* Below: add new potential terms */ 
    imat pt = terms.row(nd-1) -(terms.head_rows(nd)).each_row(); 
    uvec gt = (sum(pt.head_rows(nd),1) == 0); 
    gt %= (sum(abs(pt.head_rows(nd)),1) == 2); 
    for (unsigned int l = 0; l < d; ++l) { 
      if(terms(nd-1,l) < maxlevel(l)){//read max from om later 
        uvec gt2 = gt; 
        gt2 %= ((pt.col(l)).head(nd) == -1); 
        unsigned int h3 = sum((terms.row(nd-1)>0))+(terms(nd-1,l)<1); 
        unsigned int h4 = 1+sum(gt2); 
        if (h3==h4) { 
          pterms.row(np) = terms.row(nd-1); 
          pterms(np,l) += 1; 
          indsh = conv_to<uvec>::from(knotptst.head(d).as_row()+pterms.row(np)); 
          ptv(np) = sum(basisvar.elem(indsh)); 
          np +=1; 
        } 
      } 
    } 
  } 
  return conv_to<umat>::from(terms); 
} 


/* 
 **************************************************************************** 
 **************************************************************************** 
 *************************** outerbase ************************************** 
 **************************************************************************** 
 **************************************************************************** 
 */ 


/* 
 * outerbase::outerbase 
 * 
 * initialize the basis, which is tied to outermod and the prediction point 
 */ 

outerbase::outerbase(const outermod& om_, mat xp_) :  
  om(om_), xp(xp_)
{ 
  dograd = true;
  n_row = xp.n_rows; 
  nthreads = omp_get_num_procs(); 
  
  outerbase::build();  
}  


/* 
 * outerbase::outerbase 
 * 
 * initialize the basis, which is tied to outermod and the prediction point 
 */ 

outerbase::outerbase(const outermod& om_, mat xp_, bool dograd_) :  
  om(om_), xp(xp_)
{ 
  dograd = dograd_;
  n_row = xp.n_rows; 
  nthreads = omp_get_num_procs(); 
  
  outerbase::build();  
}  


/* 
 * outerbase::setvals_ 
 * 
 * (re)size the basis matrix and other artifacts using current outermod 
 */ 

void outerbase::setvals_(){ 
  // borrowing current info from outermod 
  d = om.d;   
  hypmatch = om.hypmatch;   
  hypst = om.hypst;   
  gest = om.gest; 
  knotptst = om.knotptst;   
  n_hyp = om.hypmatch.n_elem; 
  
  //reset loopsize in case chunk size changed 
  loopsize = (n_row+chunksize-1)/chunksize; // used for omp 
  vertpl = loopsize > 20; //hardcoded after some testing 
  //only split vertically for very large datasets 
} 

/* 
 * outerbase::setsizes_ 
 * 
 * (re)size the basis matrix and other artifacts using current outermod 
 */ 

void outerbase::setsizes_(){ 
  if (basemat.n_rows != xp.n_rows || basemat.n_cols != om.knotpt.n_elem)  
    basemat.set_size(xp.n_rows, om.knotpt.n_elem);   
  if (basemat_gradhyp.n_rows != xp.n_rows ||  
      basemat_gradhyp.n_cols != om.knotptstge[d])   
    basemat_gradhyp.set_size(xp.n_rows, om.knotptstge[d]);  
  if (basematsq.n_rows != xp.n_rows || basematsq.n_cols != om.knotpt.n_elem)  
    basematsq.set_size(xp.n_rows, om.knotpt.n_elem);
  if (basematsq_gradhyp.n_rows != xp.n_rows ||  
      basematsq_gradhyp.n_cols != om.knotptstge[d])   
    basematsq_gradhyp.set_size(xp.n_rows, om.knotptstge[d]);  
  if (basescalemat.n_rows != xp.n_rows ||  
      basescalemat.n_cols != d)
    basescalemat.set_size(xp.n_rows, d);  
  
  basescale.set_size(xp.n_rows);  
  basescale.ones();  
  basescalesq.set_size(xp.n_rows);  
} 

/* 
 * outerbase::rebuild 
 * 
 * (re)build the basis matrix and other artifacts using current outermod 
 */ 

void outerbase::build() {   
  // set values from source 
  setvals_(); // setting values 
  setsizes_(); // resizing 
  
#pragma omp parallel num_threads(nthreads)//start up an omp region 
{ 
  mat R;  
  cube Rt; 
  if (loopsize > d) { //if we have a short outerbase, just loop over d 
    mat xp_; //local xp 
    uword startind; 
    uword endind; 
    #pragma omp for nowait
    for (uword j = 0; j < loopsize; ++j) {  // if we have a tall n_row,  
      // just loop over d 
      startind = j*chunksize;  
      endind = std::min((j+1)*chunksize-1,n_row-1); //same as above,  
      //but chunked out by chunksize 
      xp_ = xp.rows(startind, endind); 
      for (uword k = 0; k < d; ++k) { 
        if(dograd) om.buildob(R, Rt, xp_, k); 
        else om.buildob(R, xp_, k);
        
        basescalemat.col(k).rows(startind, endind) = R.col(0);
        basescale.rows(startind, endind) %= R.col(0); 
        R.col(0).ones();
        basemat.submat(startind, knotptst[k],
                       endind, knotptst[k+1]-1) = R;
        //save the squared vers.
        basematsq.submat(startind, knotptst[k],
                         endind, knotptst[k+1]-1) = square(R);
        if(dograd){
        for (uword l = om.hypst[k]; l < hypst[k+1]; ++l) {
          basemat_gradhyp.submat(startind, gest[l],
                                 endind, gest[l+1]-1) =\
                                   Rt.slice(l-hypst[k]);
          basematsq_gradhyp.submat(startind, gest[l],
                                   endind, gest[l+1]-1) = 2*
                                     (Rt.slice(l-hypst[k]) % R); 
        }
        }
      }
      basescalesq.rows(startind, endind) =\
        square(basescale.rows(startind, endind));//save the squared vers.
    }
  } else { 
    #pragma omp for nowait
    for (uword k = 0; k < d; ++k) {   
      if(dograd) om.buildob(R, Rt, xp, k); 
      else om.buildob(R, xp, k);
      
      basescalemat.col(k) = R.col(0);
      basescale %= R.col(0); //keep the first column multipled for scaling 
      R.col(0).ones();
      basemat.cols(om.knotptst[k],om.knotptst[k+1]-1) = R; //put it in the 
        //right  place 
      basematsq.cols(om.knotptst[k],om.knotptst[k+1]-1) = square(R); 
      if(dograd){
      for (uword l = om.hypst[k]; l < hypst[k+1]; ++l) {
        basemat_gradhyp.cols(gest[l],gest[l+1]-1) =\
                                 Rt.slice(l-hypst[k]);
        basematsq_gradhyp.cols(gest[l], gest[l+1]-1) =\
          2*(Rt.slice(l-hypst[k]) % R); 
      }
      }
    }
    basescalesq = square(basescale);//save the squared vers. 
  } 
} 
} 

/* 
 * outerbase::getbase
 * 
 *
 */ 

mat outerbase::getbase(const uword dim) const {
  uword k = dim-1;
  mat out = basemat.cols(knotptst[k],knotptst[k+1]-1);
  out.each_col() %= basescalemat.col(k);
  return out;
} 



/* 
 * outerbase::getmat
 * 
 *
 */ 

mat outerbase::getmat(const umat& terms) const {
  mat out;
  getm_(out, terms, basemat, basescale, knotptst, 
          vertpl, chunksize, loopsize, nthreads); 
  return out;
} 


/* 
 * outerbase::getmat
 * 
 *
 */ 

cube outerbase::getmat_gradhyp(const umat& terms) const {
  cube outge;
  getmge_(outge, terms, basemat, basescale, knotptst, 
          basemat_gradhyp, gest, hypmatch, 
        vertpl, chunksize, loopsize, nthreads); 
  return outge;
} 

/* 
 * outerbase::mm 
 * 
 * in-place matrix multiply to out for basis with terms by a 
 */ 

void outerbase::mm(vec& out, const umat& terms, const vec& a) const {   
  prodmm_(out, terms, a, basemat, basescale, knotptst, 
          vertpl, chunksize, loopsize, nthreads); 
} 

/* 
 * outerbase::mm_out 
 * 
 * mm with return 
 */ 

vec outerbase::mm_out(const umat& terms, const vec& a) const {   
  vec out; 
  outerbase::mm(out, terms, a); 
  return out; 
} 

/* 
 * outerbase::tmm 
 * 
 * in-place transpose matrix multiply to out for basis with terms by a 
 */ 

void outerbase::tmm(vec& out, const umat& terms, const vec& a) const {  
  tprodmm_(out, terms, a, basemat, basescale, knotptst, 
           vertpl, chunksize, loopsize, nthreads); 
} 


/* 
 * outerbase::tmm_out 
 * 
 * tmm with return 
 */ 

vec outerbase::tmm_out(const umat& terms, const vec& a) const {   
  vec out; 
  outerbase::tmm(out, terms, a); 
  return out; 
} 

/* 
 * outerbase::mm_gradhyp 
 * 
 * in-place gradiant of matrix multiply to outge for basis with terms by a 
 * out is the matrix multiply 
 */ 

void outerbase::mm_gradhyp(vec& out, mat& out_gradhyp, 
                           const umat& terms, const vec& a) const {  
  prodmmge_(out, out_gradhyp, terms, a, 
            basemat, basescale, knotptst, 
            basemat_gradhyp, gest, hypmatch, 
            vertpl, chunksize, loopsize, nthreads); 
} 

/* 
 * outerbase::mm_gradhyp_out 
 * 
 * outerbase::mm_gradhyp with return 
 */ 

mat outerbase::mm_gradhyp_out(const umat& terms, const vec& a) const {   
  vec out; 
  mat out_gradhyp; 
  outerbase::mm_gradhyp(out, out_gradhyp, terms, a); 
  return out_gradhyp; 
} 

/* 
 * outerbase::tmm_gradhyp 
 * 
 * in-place gradiant of transpose matrix multiply to outge for basis with terms  
 * by a.  out is the matrix multiply. 
 * 
 */ 


void outerbase::tmm_gradhyp(vec& out, mat& out_gradhyp, 
                            const umat& terms, const vec& a) const {  
  tprodmmge_(out, out_gradhyp, terms, a, 
             basemat, basescale, knotptst, 
             basemat_gradhyp, gest, hypmatch, 
             vertpl, chunksize, loopsize, nthreads); 
} 


/* 
 * outerbase::mm_gradhyp_out 
 * 
 * outerbase::mm_tgradhyp with return 
 */ 


mat outerbase::tmm_gradhyp_out(const umat& terms, const vec& a) const {  
  vec out; 
  mat out_gradhyp; 
  outerbase::tmm_gradhyp(out, out_gradhyp, terms, a); 
  return out_gradhyp; 
} 

/* 
 * outerbase::sqmm 
 * 
 * square the basis matrix at terms and multiply by a 
 */ 

vec outerbase::sqmm(const umat& terms, const vec& a) const {  
  vec out; 
  prodmm_(out, terms, a, 
          basematsq, basescalesq, knotptst, 
          vertpl, chunksize, loopsize, nthreads); 
  return out; 
} 

/* 
 * outerbase::sqmm_gradhyp 
 * 
 * gradient with respect to hyperparameters of sqmm 
 */ 

mat outerbase::sqmm_gradhyp(const umat& terms, const vec& a) const {  
  vec out; 
  mat outge;  
  
  prodmmge_(out, outge, terms, a, 
            basematsq, basescalesq, knotptst, 
            basematsq_gradhyp, gest, hypmatch, 
            vertpl, chunksize, loopsize, nthreads); 
  
  return outge; 
} 

/* 
 * outerbase::sqtmm 
 * 
 * square the basis matrix at terms, transpose it and multiply by a 
 */ 

vec outerbase::sqtmm(const umat& terms, const vec& a) const {  
  vec out; 
  tprodmm_(out, terms, a, 
           basematsq, basescalesq, knotptst, 
           vertpl, chunksize, loopsize, nthreads); 
  return out; 
} 


/* 
 * outerbase::sqtmmm 
 * 
 * square the basis matrix at terms, transpose it and multiply by matrix a 
 */ 

mat outerbase::sqtmmm(const umat& terms, const mat& a) const {  
  mat out; 
  tprodmm_(out, terms, a, 
           basematsq, basescalesq, knotptst, 
           vertpl, chunksize, loopsize, nthreads); 
  return out; 
} 

/* 
 * outerbase::sqtmm_gradhyp
 * 
 * gradient with respect to hyperparameters of sqtmm 
 */ 

mat outerbase::sqtmm_gradhyp(const umat& terms, const vec& a) const {  
  vec out; 
  mat outge;  
  
  tprodmmge_(out, outge, terms, a, 
             basematsq, basescalesq, knotptst, 
             basematsq_gradhyp, gest, hypmatch, 
             vertpl, chunksize, loopsize, nthreads); 
  
  return outge; 
} 

/* 
 * outerbase::sqcolsums 
 * 
 * square the basis matrix at terms and take column sums 
 */ 

vec outerbase::sqcolsums(const umat& terms) const {  
  vec ro(n_row); 
  ro.ones(); 
  return outerbase::sqtmm(terms, ro); 
} 

/* 
 * outerbase::sqcolsums_gradhyp 
 * 
 * gradient with respect to hyperparameters of sqcolsums 
 */ 

mat outerbase::sqcolsums_gradhyp(const umat& terms) const { 
  vec ro(n_row); 
  ro.ones(); 
  return outerbase::sqtmm_gradhyp(terms, ro);
} 



/* 
 * outerbase::residvar 
 * 
 * gradient with respect to hyperparameters of sqtmm 
 */ 

vec outerbase::residvar(const umat& terms) const {  
  vec out; 
  out = 1 - outerbase::sqmm(terms, om.getvar(terms)); 
  // assumes that it is a correlation function
  // om.totvar(xp) should be used otherwise
  return out; 
} 


/* 
 * outerbase::residvar_gradhyp 
 * 
 * gradient with respect to hyperparameters of sqtmm 
 */ 

mat outerbase::residvar_gradhyp(const umat& terms) const {   
  mat outge;  
  
  vec varc = om.getvar(terms);  
  
  outge = -outerbase::sqmm_gradhyp(terms, varc);  
  // assumes that it is a correlation function
  
  mat lvarc_gradhyp2 = om.getlvar_gradhyp(terms); 
  lvarc_gradhyp2.each_col() %= varc;
  mat lvarc_gradhyp3; 
  prodmm_(lvarc_gradhyp3, terms, lvarc_gradhyp2, 
          basematsq, basescalesq, knotptst, 
          vertpl, chunksize, loopsize, nthreads); 
  outge -= lvarc_gradhyp3; 
  // assumes that it is a correlation function
  
  return outge;  
}  
