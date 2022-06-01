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

double covf::lpdf(vec hypp) const {
  double out = 0;
  
  // return infs if needed
  if(hyp.n_elem != hypp.n_elem) return (-datum::inf);
  
  for (unsigned int l = 0; l <  hypp.n_elem; ++l) {
    if(hypub[l] < hypp[l])  return (-datum::inf);
    if(hyplb[l] > hypp[l])  return (-datum::inf);
    out += 5*log(hypub[l]-hypp[l]);
    out += 5*log(hypp[l]-hyplb[l]);
  }
  out -= 0.5*accu(square(hypp-hyp0)/hypvar);
  
  return out;
}


vec covf::lpdf_gradhyp(vec hypp) const {
  vec out;
  out.set_size(hyp.n_elem);
  out.zeros();
  
  // return infs if needed
  if(hyp.n_elem != hypp.n_elem) return out;
  
  for (unsigned int l = 0; l <  hypp.n_elem; ++l) {
    if(hypub[l] < hypp[l])  return out;
    if(hyplb[l] > hypp[l])  return out;
    out[l] -= 5/(hypub[l]-hypp[l]);
    out[l] += 5/(hypp[l]-hyplb[l]);
  }
  out -= (hypp-hyp0)/hypvar;
  
  return out;
}


/* 
 **************************************************************************** 
 **************************************************************************** 
 ********************************* covf_mat25 ******************************* 
 **************************************************************************** 
 **************************************************************************** 
 */ 

/*
 * covf_mat25::covf_mat25
 *
 * constructor for covf_mat25, setting upper and lower bounds
 */

covf_mat25::covf_mat25() 
  : covf()
{
  numhyp = 1;
  
  hypnames = {"scale"};
  
  hyp.set_size(1);
  hyp(0) = 0;
  
  hyplb.set_size(1);
  hyplb(0) = -2.25;
  
  hypub.set_size(1);
  hypub(0) = 1.5;
  
  hyp0.set_size(1);
  hyp0(0) = 0;
  
  hypvar.set_size(1);
  hypvar(0) = 0.1;
  
  lowbnd = 0;
  uppbnd = 1;
}

mat covf_mat25::cov(const vec& x1, const vec& x2) const {
  double expLS = exp(a*hyp(0));
  
  mat h(x1.n_elem,x2.n_elem);
  vec x1t = x1 / expLS;
  vec x2t = x2 / expLS;
  h.each_col() = x1t;
  h.each_row() -= x2t.t();
  
  h = abs(h);
  h = (1 + h + square(h)/3) % exp(-h);
  
  return h;
}

vec covf_mat25::covmdiag(const vec& x) const {
  vec h(x.n_elem);
  h.ones();
  return h;
}

cube covf_mat25::cov_gradhyp(const vec& x1, const vec& x2) const {
  double expLS = exp(a*hyp(0));
  
  vec x1t = x1 / expLS;
  vec x2t = x2 / expLS;
  mat h(x1.n_elem,x2.n_elem);
  cube h_grad(x1.n_elem,x2.n_elem,1);
  
  h.each_col() = x1t;
  h.each_row() -= x2t.as_row();
  
  mat h2 = (h % (1+abs(h)) % exp(-abs(h)));
  
  h_grad.slice(0) = a/3 * (h % h2);
  
  return h_grad;
}

/* 
 **************************************************************************** 
 **************************************************************************** 
 ******************************* covf_mat25pow ****************************** 
 **************************************************************************** 
 **************************************************************************** 
 */ 

/*
 * covf_mat25pow::covf_mat25pow
 *
 * constructor for covf_mat25pow, setting upper and lower bounds
 */

covf_mat25pow::covf_mat25pow() 
  : covf()
{
  numhyp = 2;
  
  hypnames = {"scale","power"};
  
  hyp.set_size(2);
  hyp[0] = 0;
  hyp[1] = 0;
  
  hyplb.set_size(2);
  hyplb[0] = -2.25;
  hyplb[1] = -1.25;
  
  hypub.set_size(2);
  hypub[0] = 1.5;
  hypub[1] = 1.25;
  
  hyp0.set_size(2);
  hyp0[0] = 0;
  hyp0[1] = 0;
  
  hypvar.set_size(2);
  hypvar[0] = 0.1;
  hypvar[1] = 0.01;
  
  lowbnd = 0;
  uppbnd = 1;
}

mat covf_mat25pow::cov(const vec& x1, const vec& x2) const {
  double powv = exp(b*hyp(1));
  double expLS = exp(a*hyp(0)+b*hyp(1));
  
  
  vec x1t = pow(x1, powv)/expLS;
  vec x2t = pow(x2, powv)/expLS;
  mat h(x1.n_elem,x2.n_elem);
  h.each_col() = x1t;
  h.each_row() -= x2t.t();
  
  h = abs(h);
  h = (1 + h + square(h)/3) % exp(-h);
  
  return h;
}

vec covf_mat25pow::covmdiag(const vec& x) const {
  vec h(x.n_elem);
  h.ones();
  return h;
}

cube covf_mat25pow::cov_gradhyp(const vec& x1, const vec& x2) const {
  double powv = exp(b*hyp(1));
  double expLS = exp(a*hyp(0)+b*hyp(1));
  
  vec x1t = pow(x1, powv)/expLS;
  vec x2t = pow(x2, powv)/expLS;
  mat h(x1.n_elem,x2.n_elem);
  cube h_grad(x1.n_elem,x2.n_elem,2);
  h_grad.zeros();
  h.each_col() = x1t;
  h.each_row() -= x2t.as_row();
  
  mat h2 = (h % (1+abs(h)) % exp(-abs(h)));
  
  h_grad.slice(1).each_col() = (log(x1)) % x1t;
  h_grad.slice(1).each_row() -= (log(x2) % x2t).as_row();
  
  h_grad.slice(1) %= -(b * powv/3) * h2;
  h %= h2;
  h_grad.slice(1) += (b/3) * h;
  h_grad.slice(0) = (a/3) * h;
  
  return h_grad;
}

/* 
 **************************************************************************** 
 **************************************************************************** 
 ******************************* covf_mat25ang ****************************** 
 **************************************************************************** 
 **************************************************************************** 
 */ 


covf_mat25ang::covf_mat25ang() 
  : covf() 
{
  numhyp = 2;
  
  hypnames = {"sin.sc","cos.sc"};
  
  hyp.set_size(2);
  hyp[0] = 0;
  hyp[1] = 0;
  
  hyplb.set_size(2);
  hyplb[0] = -2.25;
  hyplb[1] = -2.25;
  
  hypub.set_size(2);
  hypub[0] = 1.5;
  hypub[1] = 1.5;
  
  hyp0.set_size(2);
  hyp0[0] = 0;
  hyp0[1] = 0;
  
  hypvar.set_size(2);
  hypvar[0] = 0.1;
  hypvar[1] = 0.1;
  
  lowbnd = 0;
  uppbnd = 6.283185;
}

mat covf_mat25ang::cov(const vec& x1, const vec& x2) const {
  vec sx1 = sin(x1);
  vec cx1 = cos(x1);
  vec sx2 = sin(x2);
  vec cx2 = cos(x2);
  double expLSs = exp(a*hyp(0));
  double expLSc = exp(a*hyp(1));
  //double crosss = tanh(eta(2));
  
  mat hs(x1.n_elem,x2.n_elem);
  mat hc(x1.n_elem,x2.n_elem);
  sx1 /= expLSs;
  sx2 /= expLSs;
  hs.each_col() = sx1;
  hs.each_row() -= sx2.t();
  cx1 /= expLSc;
  cx2 /= expLSc;
  hc.each_col() = cx1;
  hc.each_row() -= cx2.t();
  
  mat h(x1.n_elem,x2.n_elem);
  h = sqrt(square(hs) + square(hc));
  h = (1 + h + square(h)/3) % exp(-h);
  
  return h;
}

vec covf_mat25ang::covmdiag(const vec& x) const {
  vec h(x.n_elem);
  h.ones();
  return h;
}

cube covf_mat25ang::cov_gradhyp(const vec& x1, const vec& x2) const {
  vec sx1 = sin(x1);
  vec cx1 = cos(x1);
  vec sx2 = sin(x2);
  vec cx2 = cos(x2);
  double expLSs = exp(a*hyp(0));
  double expLSc = exp(a*hyp(1));
  
  mat hs(x1.n_elem,x2.n_elem);
  mat hc(x1.n_elem,x2.n_elem);
  sx1 /= expLSs;
  sx2 /= expLSs;
  hs.each_col() = sx1;
  hs.each_row() -= sx2.t();
  cx1 /= expLSc;
  cx2 /= expLSc;
  hc.each_col() = cx1;
  hc.each_row() -= cx2.t();
  
  mat h(x1.n_elem,x2.n_elem);
  h = sqrt(square(hs) + square(hc));
  
  cube dhde(x1.n_elem,x2.n_elem,2);
  dhde.slice(0) = (a/3)*square(hs);
  dhde.slice(1) = (a/3)*square(hc);
  dhde.slice(0) %= (exp(-h) % (h+1));
  dhde.slice(1) %= (exp(-h) % (h+1));
  
  return dhde;
}

