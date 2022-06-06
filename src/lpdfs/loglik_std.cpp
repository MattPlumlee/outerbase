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
 **************************************************************************** 
 **************************************************************************** 
 *************************** loglik_std *********************************** 
 **************************************************************************** 
 **************************************************************************** 
 */ 

/*
 * loglik_std::loglik_std
 *
 */

loglik_std::loglik_std(const outermod& om_, umat terms_, vec y_, mat x_)
  :  lpdf(), om(om_), ob(om,x_,true),  y(y_), x(x_)
{
  terms = terms_;
  npara = 1;
  
  basismat = ob.getmat(terms);
  basismat_gradhyp = ob.getmat_gradhyp(terms);
  
  para0.set_size(1);
  para0[0] = log(0.01*var(y));
  paravar.set_size(1);
  paravar[0] = 1;
  
  paranames = {"noisescale"};
  para = para0;
  nterms = terms.n_rows;
}


/*
 * loglik_std::updateom
 *
 */

void loglik_std::updateom() {
  ob.build();
  basismat = ob.getmat(terms);
  basismat_gradhyp = ob.getmat_gradhyp(terms);
}


/*
 * loglik_gauss::updatepara
 *
 */

void loglik_std::updatepara(vec para_) {
  para = para_;
}

/*
 * loglik_gauss::updateterms
 *
 */

void loglik_std::updateterms(umat terms_) {
  terms = terms_;
  nterms = terms.n_rows;
  basismat = ob.getmat(terms);
  basismat_gradhyp = ob.getmat_gradhyp(terms);
}


/*
 * loglik_gauss::update
 *
 */

void loglik_std::update(const vec& coeff_) {
  coeff = coeff_;
  yhat = basismat * coeff;
  mat yhatge;
  if(compute_gradhyp){
    yhatge.set_size(y.n_elem, basismat_gradhyp.n_slices);
    for(uword l = 0; l < basismat_gradhyp.n_slices; l++) 
      yhatge.col(l) = basismat_gradhyp.slice(l) * coeff;
  }
  
  vec residtemp = exp(-para[0])*(yhat-y);
  vec residtemp2 = square(residtemp);
  if(compute_val) val = -0.5*sum(residtemp2) - y.n_elem*(para[0]);
  if(compute_grad){
    grad.copy_size(coeff);
    residtemp = -exp(-para[0])*(residtemp);
    ob.tmm(grad, terms, residtemp);
    if(compute_gradhyp) gradhyp = (residtemp.t() * yhatge).as_col();
    if(compute_gradpara) gradpara = sum(residtemp2) - y.n_elem;
  }
}

/*
 * loglik_std::hessmult
 *
 */

vec loglik_std::hessmult(const vec& g) {
  vec lh = (basismat.t() * (basismat * g));
  return exp(-2*para[0])*lh;
}

/*
 * loglik_std::diaghess
 *
 */

vec loglik_std::diaghess() {
  vec lh = sum(square(basismat),0).t();
  return ((exp(-2*para(0)))*lh);
}

/*
 * loglik_std::diaghessgradhyp
 *
 */

mat loglik_std::diaghessgradhyp() {
  cube basismatsq_gradhyp = basismat_gradhyp;
  basismatsq_gradhyp.each_slice() %= 2*basismat;
  mat lh = sum(basismatsq_gradhyp,0);
  return ((exp(-2*para(0)))*lh);
}

/*
 * loglik_std::diaghessgradpara
 *
 */

mat loglik_std::diaghessgradpara() {
  vec lh = sum(square(basismat),0).t();
  return ((-2*exp(-2*para(0)))*lh);
}


/*
 * loglik_std::hess
 *
 */

mat loglik_std::hess() {
  mat lh = basismat.t() * basismat;
  return (exp(-2*para(0)))*lh;
}

/*
 * loglik_std::hessgradhyp
 *
 */

cube loglik_std::hessgradhyp() {
  cube hess_gradhyp(basismat.n_cols,
                    basismat.n_cols,
                    basismat_gradhyp.n_slices);
  
  for(uword l = 0; l < basismat_gradhyp.n_slices; l++) {
    hess_gradhyp.slice(l) = (exp(-2*para(0))) *\
      (basismat.t() * basismat_gradhyp.slice(l));
    hess_gradhyp.slice(l) += hess_gradhyp.slice(l).t();
  }
  
  return hess_gradhyp;
}

/*
 * loglik_std::hessgradpara
 *
 */

cube loglik_std::hessgradpara() {
  cube lh(basismat.n_cols, basismat.n_cols, 1);
  lh.slice(0) = basismat.t() * basismat;
  return (-2*exp(-2*para(0)))*lh;
}

predf* loglik_std::pred() const {
  return (new predr_std(*this));
}


/* 
 **************************************************************************** 
 **************************************************************************** 
 ***************************** predr_std ************************************* 
 **************************************************************************** 
 **************************************************************************** 
 */ 

predr_std::predr_std(const loglik_std& loglik) : 
om(loglik.om), para(loglik.para), terms(loglik.terms),
x(loglik.x), ob(om,x,false) //private
{
  nthreads = loglik.ob.nthreads;
  ob.nthreads = nthreads;
  coeff = loglik.coeff;
  basismat = ob.getmat(terms);
  if(!loglik.didnotothess) {
    if(loglik.didfulltothess) coeffcov = inv(loglik.tothess);
    else {
      coeffcov.set_size(coeff.n_elem, coeff.n_elem);
      coeffcov.zeros();
      coeffcov.diag() = loglik.totdiaghess;
    }
  } else {
    coeffcov.set_size(coeff.n_elem, coeff.n_elem);
    coeffcov.zeros();
  }
}

void predr_std::update(const mat& x_) {
  x = x_;
  new (&ob) outerbase(om,x_,false);
  ob.nthreads = nthreads;
  basismat = ob.getmat(terms);
}

vec predr_std::mean() const {
  return (basismat * coeff);
}

vec predr_std::var() const {
  mat adjbasmat = basismat * coeffcov;
  adjbasmat %= basismat;
  vec out = sum(adjbasmat,1);
  out += exp(2*para(0));
  return out;
}