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

/* 
 * This file is a stub!  It cannot be used alone and is only used to offer a 
 * shorter format to make files easier to edit.  See fit.h and fit.cpp in the
 * main directory.
 * 
 */

/* 
 **************************************************************************** 
 **************************************************************************** 
 *************************** loglik_gda ************************************* 
 **************************************************************************** 
 **************************************************************************** 
 */ 

/* 
 * loglik_gda::loglik_gda 
 * 
 */ 

loglik_gda::loglik_gda(const outermod& om_, umat terms_, vec y_,
                       mat x_)
                       :  lpdf(), om(om_), ob(om,x_,true) ,  y(y_), x(x_), 
                          doda(true), redostd(true)
{
  terms = terms_;
  npara = 2;
  
  para0.set_size(2);
  para0[0] = 0.5*log(0.01*var(y));
  para0[1] = 0;
  paravar.set_size(2);
  paravar[0] = 4;
  paravar[1] = 4;
  
  paranames = {"noisescale","lik.coeffscale"};
  para = para0;
  vec obsvar;
  buildstd();
  
  nterms = terms.n_rows;
}

/* 
 * loglik_gda::setnthreads
 * 
 */ 

void loglik_gda::setnthreads(int k) {
  ob.nthreads = k;
}

/* 
 * loglik_gda::updateom
 * 
 */ 

void loglik_gda::updateom() {
  ob.build();
  if(doda) redostd = true;
}

/* 
 * loglik_gda::updatepara
 * 
 */ 

void loglik_gda::updatepara(vec para_) {
  para = para_;
  redostd = true;
}

/* 
 * loglik_gda::updateterms
 * 
 */ 

void loglik_gda::updateterms(umat terms_) {
  terms = terms_;
  nterms = terms.n_rows;
  if(doda) redostd = true;
}


/*
 * loglik_gda::update
 *
 */

void loglik_gda::update(const vec& coeff_) {
  yhat.copy_size(y);
  residtemp.copy_size(y);
  coeff = coeff_;
  
  if(compute_gradhyp){
    yhatge.set_size(y.n_elem, ob.n_hyp);
    ob.mm_gradhyp(yhat, yhatge, terms, coeff);
  } else ob.mm(yhat, terms, coeff);
  
  loglik_gda::buildstd();
  
  residtemp = (yhat-y)/obssd;
  residtemp2 = square(residtemp);
  if(compute_val) val = -0.5*sum(square(residtemp)) - sum(log(obssd));
  if(compute_grad){
    grad.copy_size(coeff);
    residtemp = -1.*(residtemp/obssd);
    residtemp2 /= obssd;
    ob.tmm(grad, terms, residtemp);
    if(compute_gradhyp){
      gradhyp = (residtemp.t() * yhatge).as_col();
      if (doda) {
        gradhyp += (residtemp2.t() * obssd_gradhyp).as_col();
        gradhyp -= ((1/obssd).t() * obssd_gradhyp).as_col();
      }
    }
    if(compute_gradpara){ 
      gradpara = (residtemp2.t() * obssd_gradpara).as_col();
      gradpara -= ((1/obssd).t() * obssd_gradpara).as_col();
    }
  }
}

/*
 * loglik_gda::hessmult
 *
 */

vec loglik_gda::hessmult(const vec& g) {
  yhattemp.copy_size(y);
  ob.mm(yhattemp, terms, g);
  yhattemp /= obssd;
  yhattemp /= obssd;
  gradtemp.copy_size(grad);
  ob.tmm(gradtemp, terms, yhattemp);
  return gradtemp;
}


/* 
 * loglik_gda::loglik_gda::diaghess()
 * 
 */ 

vec loglik_gda::diaghess() {
  vec lh = ob.sqtmm(terms, 1 / square(obssd));
  return lh;
}

/* 
 * loglik_gda::diaghessgradhyp()
 * 
 */ 

mat loglik_gda::diaghessgradhyp() {
  vec temp = 1 / square(obssd);
  mat lh = ob.sqtmm_gradhyp(terms, temp);
  temp %= -2/obssd;
  
  
  if (doda) {
    mat obssd_gradhyp_temp = obssd_gradhyp;
    obssd_gradhyp_temp.each_col() %= temp;
    lh += ob.sqtmmm(terms, obssd_gradhyp_temp);
  }
  
  return lh;
}

/* 
 * loglik_gda::diaghessgradpara()
 * 
 */ 

mat loglik_gda::diaghessgradpara() {
  vec temp = 1 / square(obssd);
  temp %= -2/obssd;
  mat obssd_gradpara_temp = obssd_gradpara;
  obssd_gradpara_temp.each_col() %= temp;
  mat lh = ob.sqtmmm(terms, obssd_gradpara_temp);
  
  return lh;
}

/* 
 * loglik_gda::buildstd()
 * 
 */ 

void loglik_gda::buildstd(){
  if(redostd){
    vec obsvar;
    obsvar.set_size(y.n_elem);
    obsvar.fill(exp(2*para(0)));
    vec rterms = ob.residvar(terms);
    if(doda) obsvar += exp(2*para(1))*rterms;
    obssd = sqrt(obsvar);
    if(doda) {
      obssd_gradhyp = ob.residvar_gradhyp(terms);
      obssd_gradhyp.each_col() %= (exp(2*para(1))*0.5)/obssd;
    }
    obssd_gradpara.set_size(y.n_elem, 2);
    obssd_gradpara.col(0) = exp(2*para(0))/obssd;
    if(doda) obssd_gradpara.col(1) = exp(2*para(1))*rterms/obssd;
    else obssd_gradpara.col(1).zeros();
  }
  redostd = false;
}

predf* loglik_gda::pred() const {
  return (new pred_gda(*this));
}

/* 
 **************************************************************************** 
 **************************************************************************** 
 ***************************** pred_gda ************************************* 
 **************************************************************************** 
 **************************************************************************** 
 */ 

pred_gda::pred_gda(const loglik_gda& loglik) : 
om(loglik.om), para(loglik.para), terms(loglik.terms),
x(loglik.x), ob(om,x,false) //private
{
  nthreads = loglik.ob.nthreads;
  ob.nthreads = nthreads;
  doda = loglik.doda;
  coeff = loglik.coeff;
  
  if(!loglik.didnotothess) {
    if(loglik.didfulltothess) coeffvar = 1/loglik.totdiaghess.diag();
    else {
      coeffvar = 1/loglik.totdiaghess;
    }
  } else {
    coeffvar = 0*coeff;
  }
}

void pred_gda::update(const mat& x_) {
  x = x_;
  new (&ob) outerbase(om,x_,false);
  ob.nthreads = nthreads;
}
vec pred_gda::mean() const {
  return ob.mm_out(terms, coeff);
}
vec pred_gda::var() const {
  vec out = ob.sqmm(terms, coeffvar);
  out += exp(2*para(0));
  if (doda) out += exp(2*para(1))*ob.residvar(terms);
  
  return out;
}

