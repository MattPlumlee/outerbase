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
 *************************** logpr_gauss ************************************ 
 **************************************************************************** 
 **************************************************************************** 
 */ 

/* 
 * logpr_gauss::logpr_gauss
 * 
 */ 

logpr_gauss::logpr_gauss(const outermod& om_, umat terms_)
  :  lpdf(), om(om_)
{
  npara = 1;
  terms = terms_;
  
  para0.set_size(1);
  para0[0] = 6;
  paravar.set_size(1);
  paravar[0] = 4;
  
  paranames = {"coeffscale"};
  nterms = terms.n_rows;
  para = para0;
  sca = exp(para(0));
  coeffsd = sqrt(om.getvar(terms));
  coefflvarge = om.getlvar_gradhyp(terms);
}

/* 
 * logpr_gauss::updateom
 * 
 */ 

void logpr_gauss::updateom() {
  coeffsd = sqrt(om.getvar(terms));
  coefflvarge = om.getlvar_gradhyp(terms);
}

/* 
 * logpr_gauss::updatepara
 * 
 */ 

void logpr_gauss::updatepara(vec para_) {
  para = para_;
  sca = exp(para(0));
}

/* 
 * logpr_gauss::updateterms
 * 
 */ 

void logpr_gauss::updateterms(umat terms_) {
  terms = terms_;
  nterms = terms.n_rows;
  coeffsd = sqrt(om.getvar(terms));
  coefflvarge = om.getlvar_gradhyp(terms);
}


/* 
 * logpr_gauss::update
 * 
 */ 

void logpr_gauss::update(const vec& coeff_) {
  coeff = coeff_;
  stdresid = coeff/(coeffsd*sca);
  
  if(compute_val) val = -0.5*sum(square(stdresid)) - sum(log(coeffsd*sca));
  if(compute_gradhyp) gradhyp = (0.5*coefflvarge).t() * (square(stdresid)-1);
  if(compute_gradpara) gradpara = sum(square(stdresid)) - coeffsd.n_elem;
  if(compute_grad) grad = -1.*stdresid/(coeffsd*sca);
}

/* 
 * logpr_gauss::hessmult
 * 
 */ 

vec logpr_gauss::hessmult(const vec& g) {
  return (g / square(coeffsd*sca));
}

/* 
 * logpr_gauss::diaghess
 * 
 */ 

vec logpr_gauss::diaghess() {
  return (1./square(coeffsd*sca));
}

/* 
 * logpr_gauss::diaghessgradhyp
 * 
 */ 

mat logpr_gauss::diaghessgradhyp() {
  mat out = coefflvarge;
  out.each_col() /= square(coeffsd*sca);
  return (-out);
}


/* 
 * logpr_gauss::diaghessgradpara
 * 
 */ 

mat logpr_gauss::diaghessgradpara() {
  return ((-2./square(coeffsd*sca)));
}


/* 
 * logpr_gauss::diaghess
 * 
 */ 

mat logpr_gauss::hess() {
  mat hess(nterms, nterms);
  hess.zeros();
  hess.diag() = 1./square(coeffsd*sca);
  return hess;
}

/* 
 * logpr_gauss::diaghessgradhyp
 * 
 */ 

cube logpr_gauss::hessgradhyp() {
  cube hess_gradhyp(nterms, nterms, coefflvarge.n_cols);
  hess_gradhyp.zeros();
  mat out = coefflvarge;
  out.each_col() /= square(coeffsd*sca);
  for(uword l = 0; l < hess_gradhyp.n_slices; l++) 
    hess_gradhyp.slice(l).diag() = -out.col(l);
  return hess_gradhyp;
}


/* 
 * logpr_gauss::diaghessgradpara
 * 
 */ 

cube logpr_gauss::hessgradpara() {
  cube hess_gradpara(nterms, nterms, 1);
  hess_gradpara.zeros();
  hess_gradpara.slice(0).diag() = ((-2./square(coeffsd*sca)));
  return (hess_gradpara);
}