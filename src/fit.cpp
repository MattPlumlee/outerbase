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
using namespace Rcpp;
using namespace arma;  

#include "covfuncs.h"
#include "modandbase.h"
#include "fit.h"

void lpdf::optcg(double tol, unsigned int maxepch) {
  fullhess = false;
  
  compute_val = true;
  compute_grad = true;
  compute_gradhyp = false;
  compute_gradpara= false;
  
  if(coeff.n_elem != nterms){
    coeff.set_size(nterms);
    coeff.zeros();
  }
  update(coeff);
  
  vec m = diaghess();
  
  if(!(m.is_finite()) && !(grad.is_finite())){
    val = (-datum::inf);
    return;
  }
  vec rm = grad/m;
  
  double num = 0;
  
  vec p = rm;
  vec q = hessmult(p);
  
  double beta = 0;
  double denom = 1;
  double alpha = 0;
  double num2 = 0;
  double valo = val;
  double valdiff = 10;
  unsigned int k;
  for (k = 0; k < maxepch; k++) {
    num = accu(grad % rm);
    if (num < tol && valdiff < tol) break;
    denom = accu(q%p);
    alpha = num/denom;
    coeff +=  alpha*p;
    valo = val;
    update(coeff);
    valdiff = val-valo;
    rm = grad/m;
    num2 = -accu((alpha*q) % rm);
    beta = num2/num;
    p = rm+beta*p;
    q = hessmult(p);
  }  
  
  compute_gradhyp = true;
  compute_gradpara= true;
  
  update(coeff);
  
  compute_gradhyp = false;
  compute_gradpara= false;
  
  return;
}

void lpdf::optnewton() { 
  fullhess = true;
  
  compute_val = true;
  compute_grad = true;
  compute_gradhyp = false;
  compute_gradpara= false;
  
  if(coeff.n_elem != nterms){
    coeff.set_size(nterms);
    coeff.zeros();
  }
  update(coeff);
  
  mat h = hess();
  vec r = grad;
  
  if(!(h.is_finite()) && !(r.is_finite())){
    val = (-datum::inf);
    return;
  }
  
  coeff += solve(h,r);
  
  compute_gradhyp = true;
  compute_gradpara= true;
  
  update(coeff);
  
  compute_gradhyp = false;
  compute_gradpara= false;
  
  return;
}

double lpdf::paralpdf(vec parap) const {
  double out = 0;
  
  // return infs if needed
  if(npara != parap.n_elem) return (-datum::inf);
  
  out -= 0.5*accu(square(parap-para0)/paravar);
  
  return out;
}


vec lpdf::paralpdf_grad(vec parap) const {
  vec out;
  out.set_size(para.n_elem);
  out.zeros();
  
  // return infs if needed
  if(npara != parap.n_elem) return out;
  
  out -= (parap-para0)/paravar;
  
  return out;
}

/* 
 **************************************************************************** 
 **************************************************************************** 
 ******************************* lpdfvec ************************************ 
 **************************************************************************** 
 **************************************************************************** 
 */ 


/*
 * lpdfvec::lpdfvec
 *
 * constructor for a vector of logpdfs, useful for combining prior and 
 * likelihoods
 */

lpdfvec::lpdfvec(lpdf& loglik_, lpdf& logpr_) //names are arbitrary...
  :  lpdf()
{
  parasrt.set_size(2); //only two entries in the vector for now
  paraend.set_size(2);
  
  terms = loglik_.terms;
  nterms = loglik_.nterms; //match the number of terms
  lpdflist.push_back(loglik_);
  parasrt(0) = 0;
  paraend(0) = loglik_.npara-1;
  
  lpdflist.push_back(logpr_);
  parasrt(1) = paraend(0)+1;
  paraend(1) = parasrt(1)+logpr_.npara-1;
  
  para.set_size(1+paraend(1));
  
  para.subvec(parasrt[0],paraend[0]) = loglik_.para;
  paranames = loglik_.paranames;
  para.subvec(parasrt[1],paraend[1]) = logpr_.para;
  paranames.insert(paranames.end(),
                   logpr_.paranames.begin(),
                   logpr_.paranames.end());
}


/*
 * lpdfvec::updateom
 *
 * update the lpdfs with new outermod
 */

void lpdfvec::updateom() {
  for (lpdf& lpdfptr : lpdflist) lpdfptr.updateom();
  redohess = true; //redo the hessian if we update the model
}


/*
 * lpdfvec::updateom
 *
 * update the lpdfs with new para
 */

void lpdfvec::updatepara(vec para_) {
  uword cnt = 0;
  for (lpdf& lpdf_ : lpdflist){ 
    para.subvec(parasrt[cnt],paraend[cnt]) = 
      para_.subvec(parasrt[cnt],paraend[cnt]);
    lpdf_.updatepara(para_.subvec(parasrt[cnt],paraend[cnt])); 
    cnt++;
  }
  redohess = true; //redo the hessian if we update the model
}


/*
 * lpdfvec::updateterms
 *
 * update the lpdfs with new terms
 */

void lpdfvec::updateterms(umat terms_) {
  terms = terms_;
  for (lpdf& lpdfptr : lpdflist){ //loop over all things in vector
    lpdfptr.updateterms(terms_);
    nterms = lpdfptr.nterms;
  }
  redohess = true; //redo the hessian if we update the model
}

/*
 * lpdfvec::buildhess
 *
 * build the hessian when needed
 */

void lpdfvec::buildhess(){
  if (redohess) { //if needed, rebuild the hessian
    diaghessv = lpdfvec::diaghess_(); //first hull hessian
    settotdiaghess(diaghessv);
    if (domargadj) { // save this stuff if we need marginal adjustment
                     // otherwise it is of no use.
      diaghessgradhypv = lpdfvec::diaghessgradhyp_();
      diaghessgradparav = lpdfvec::diaghessgradpara_();
      
      val_margadj = -0.5*sum( log(diaghessv) );
      mat tempm = diaghessgradhypv;
      tempm.each_col() /= diaghessv;
      gradhyp_margadj = -0.5*sum(tempm,0).as_col();
      tempm = diaghessgradparav;
      tempm.each_col() /= diaghessv;
      gradpara_margadj = -0.5*sum(tempm,0).as_col();
    }
    
    if (fullhess) {     
      hessv = lpdfvec::hess_(); //first hull hessian
      settothess(hessv);
      if (domargadj) { // save this stuff if we need marginal adjustment
                       // otherwise it is of no use.
        hessgradhypv = lpdfvec::hessgradhyp_();
        hessgradparav = lpdfvec::hessgradpara_();
        vec heigval; 
        mat heigvec; 
        eig_sym( heigval, heigvec, hessv);
        mat hessi = heigvec.t();
        hessi.each_col() /= heigval;
        hessi = heigvec * hessi;
        
        val_margadj = -0.5*sum( log(heigval) );
        mat tempm;
        gradhyp_margadj.set_size(hessgradhypv.n_slices);
        for (uword l = 0; l < hessgradhypv.n_slices; l++){
          tempm = hessgradhypv.slice(l);
          tempm %= hessi;
          gradhyp_margadj(l) = -0.5*accu(tempm);
        }
        gradpara_margadj.set_size(hessgradparav.n_slices);
        for (uword l = 0; l < hessgradparav.n_slices; l++){
          tempm = hessgradparav.slice(l);
          tempm %= hessi;
          gradpara_margadj(l) = -0.5*accu(tempm);
        }
      }
    }  
  }
  redohess = false; //reset until needed
}

/*
 * lpdfvec::setnthreads
 *
 * 
 */

void lpdfvec::setnthreads(int k) {
  for (lpdf& lpdfptr : lpdflist){ //loop over all things in vector
    lpdfptr.setnthreads(k);
  }
}

  
/*
 * lpdfvec::margadj
 *
 * do the marginal adjustment when combining models
 */

void lpdfvec::update(const vec& coeff_){
  coeff = coeff_;
  for (lpdf& lpdfptr : lpdflist) { //sync up computing
    lpdfptr.compute_val = compute_val;
    lpdfptr.compute_grad = compute_grad;
    lpdfptr.compute_gradhyp = compute_gradhyp;
    lpdfptr.compute_gradpara = compute_gradpara;
  }
  
  // loop through and update all the lpdfs
  for (lpdf& lpdfptr : lpdflist)  lpdfptr.update(coeff); 
  
  // reset and resize all computations that are needed
  if(compute_val) val = 0;
  if(compute_grad){
    grad.copy_size(lpdflist[0].get().grad);
    grad.zeros();
  }
  if(compute_gradhyp){
    gradhyp.copy_size(lpdflist[0].get().gradhyp);
    gradhyp.zeros();
  } 
  if(compute_gradpara){
    gradpara.copy_size(para);
    gradpara.zeros();
  }
  buildhess();
  
  uword cnt = 0;
  for (lpdf& lpdfptr : lpdflist) {  //loop through all lpdfs
    if(compute_val) val += lpdfptr.val;
    if(compute_grad) grad += lpdfptr.grad;
    if(compute_gradhyp) gradhyp += lpdfptr.gradhyp;
    if(compute_gradpara)
      gradpara.subvec(parasrt[cnt],paraend[cnt]) += lpdfptr.gradpara;
    cnt++;
  }
  
  if (domargadj) margadj(); // do the marginal adjustment of requested
                            // this is only for the diagonal
}


/*
 * lpdfvec::margadj
 *
 * do the marginal adjustment when combining models
 */

void lpdfvec::margadj () {
  
  //
  if(compute_val) val += val_margadj;
  if(compute_gradhyp) gradhyp += gradhyp_margadj;
  if(compute_gradpara){
    gradpara += gradpara_margadj;
  }
}


/*
 * pdfvec::hessmult
 *
 * multiply by hessian
 */

vec lpdfvec::hessmult(const vec& g) {
  vec out;
  uword cnt = 0;
  for (lpdf& lpdfptr : lpdflist) {  //loop through all lpdfs
    if(cnt==0) out = lpdfptr.hessmult(g);
    else out += lpdfptr.hessmult(g);
    cnt++;
  }
  return out;
}

/*
 * lpdfvec::diaghess
 *
 * return the diagonal of the hessian
 */

vec lpdfvec::diaghess(){
  return diaghessv;  
}

/*
 * lpdfvec::diaghessgradhyp
 *
 * return the diagonal of the hessian gradient with respect to hyp
 */

mat lpdfvec::diaghessgradhyp(){
  return  diaghessgradhypv;
}

/*
 * lpdfvec::diaghessgradpara
 *
 * return the diagonal of the hessian gradient with respect to para
 */

mat lpdfvec::diaghessgradpara(){
  return diaghessgradparav;
}


/*
 * lpdfvec::diaghess
 *
 * return the diagonal of the hessian
 */

mat lpdfvec::hess(){
  return hessv;  
}

/*
 * lpdfvec::diaghessgradhyp
 *
 * return the diagonal of the hessian gradient with respect to hyp
 */

cube lpdfvec::hessgradhyp(){
  return hessgradhypv;
}

/*
 * lpdfvec::diaghessgradpara
 *
 * return the diagonal of the hessian gradient with respect to para
 */

cube lpdfvec::hessgradpara(){
  return hessgradparav;
}


/*
 * 
 *
 * 
 */

double lpdfvec::paralpdf(vec parap) const {
  double out = 0;
  uword cnt = 0;
  for (lpdf& lpdfptr : lpdflist) {
    out += lpdfptr.paralpdf(parap.subvec(parasrt[cnt],paraend[cnt]));
    cnt++;
  }
  return out;
}


/*
 * 
 *
 * 
 */

vec lpdfvec::paralpdf_grad(vec parap) const {
  vec out(parap.n_elem);
  uword cnt = 0;
  for (lpdf& lpdfptr : lpdflist) {
    out.subvec(parasrt[cnt],paraend[cnt]) = 
      lpdfptr.paralpdf_grad(parap.subvec(parasrt[cnt],paraend[cnt]));
    cnt++;
  }
  return out;
}


/*
 * lpdfvec::hess_
 *
 * compute the hessian for use later
 */

mat lpdfvec::hess_(){
  mat out;
  uword cnt = 0;
  for (lpdf& lpdfptr : lpdflist) {  //loop through all lpdfs
    if(cnt==0) out = lpdfptr.hess();
    else out += lpdfptr.hess();
    cnt++;
  }
  return out;
}

/*
 * lpdfvec::diaghessgradhyp_
 *
 * compute the hessian_gradhyp for use later
 */

cube lpdfvec::hessgradhyp_() {
  cube out;
  uword cnt = 0;
  out.zeros();
  for (lpdf& lpdfptr : lpdflist) { //loop through all lpdfs
    if(cnt==0) out = lpdfptr.hessgradhyp();
    else out += lpdfptr.hessgradhyp();
    cnt++;
  }
  return out;
}

/*
 * lpdfvec::diaghessgradpara_
 *
 * compute the hessian_gradpara for use later
 */

cube lpdfvec::hessgradpara_() {
  cube out;
  out.set_size(lpdflist[0].get().nterms, lpdflist[0].get().nterms, para.n_elem);
  out.zeros();
  uword cnt = 0;
  for (lpdf& lpdfptr : lpdflist) {
    out.slices(parasrt[cnt],paraend[cnt]) = lpdfptr.hessgradpara();
    cnt++;
  }
  return out;
}


/*
 * lpdfvec::diaghess_
 *
 * compute the hessian for use later
 */

vec lpdfvec::diaghess_(){
  vec out;
  uword cnt = 0;
  for (lpdf& lpdfptr : lpdflist) {  //loop through all lpdfs
    if(cnt==0) out = lpdfptr.diaghess();
    else out += lpdfptr.diaghess();
    cnt++;
  }
  return out;
}

/*
 * lpdfvec::diaghessgradhyp_
 *
 * compute the hessian_gradhyp for use later
 */

mat lpdfvec::diaghessgradhyp_(){
  mat out;
  uword cnt = 0;
  out.zeros();
  for (lpdf& lpdfptr : lpdflist) { //loop through all lpdfs
    if(cnt==0) out = lpdfptr.diaghessgradhyp();
    else out += lpdfptr.diaghessgradhyp();
    cnt++;
  }
  return out;
}

/*
 * lpdfvec::diaghessgradpara_
 *
 * compute the hessian_gradpara for use later
 */

mat lpdfvec::diaghessgradpara_(){
  mat out;
  out.set_size(lpdflist[0].get().nterms, para.n_elem);
  out.zeros();
  uword cnt = 0;
  for (lpdf& lpdfptr : lpdflist) {
    out.cols(parasrt[cnt],paraend[cnt]) = lpdfptr.diaghessgradpara();
    cnt++;
  }
  return out;  
}

void lpdfvec::settotdiaghess(vec diaghess){
  totdiaghess = diaghess;
  for (lpdf& lpdfptr : lpdflist) lpdfptr.settotdiaghess(diaghess);
}

void lpdfvec::settothess(mat hess){
  tothess = hess;
  for (lpdf& lpdfptr : lpdflist) lpdfptr.settothess(hess);
}


#include "lpdfs/loglik_std.cpp"   
#include "lpdfs/logpr_gauss.cpp"
#include "lpdfs/loglik_gauss.cpp"
#include "lpdfs/loglik_gda.cpp"



/*

*/

