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

using namespace Rcpp;
using namespace arma;

#include "linalg.h"

/*
 * domult_
 *
 * supporting function for prodmm_
 */

void domult_(vec& out, const vec& a, vec& temp, 
             const umat& terms, const uvec& knotptst, const mat& basemat,
             int num_threads) {
  if(out.n_elem !=  basemat.n_rows) 
    out.set_size(basemat.n_rows); 
  if(temp.n_elem !=  basemat.n_rows) 
    temp.set_size(basemat.n_rows); 
  
  out.zeros();
  uword termsnrows = terms.n_rows;
  uword termsncols = terms.n_cols;
  if (omp_in_parallel()) {
    for (uword k = 0; k < termsnrows; ++k) { 
      temp.fill(a[k]); 
      for (uword l = 0; l < termsncols; ++l)  
        if(terms.at(k,l)>0) temp %= basemat.col(knotptst[l]+terms.at(k,l)); 
      out += temp;  
    }
  } else {
  #pragma omp parallel num_threads(num_threads)
  {
  vec out_ = out;
  vec temp_ = temp;
  out_.zeros();
  #pragma omp for
  for (uword k = 0; k < termsnrows; ++k) { 
    temp_.fill(a[k]); 
    for (uword l = 0; l < termsncols; ++l)  
      if(terms.at(k,l)>0) temp_ %= basemat.col(knotptst[l]+terms.at(k,l)); 
    out_ += temp_;  
  } 
  #pragma omp critical  
  out += out_;  
  }
  }
}

/*
 * prodmm_
 *
 * does matrix multiplication in a memory friendly way that leverages omp.
 */

void prodmm_(vec& out, 
             const umat& terms, 
             const vec& a,
             const mat& basemat, const vec& basescale, const uvec& knotptst,
             bool vertpl, uword chunksize, uword loopsize, int num_threads) {
  if(out.n_elem !=  basemat.n_rows) out.set_size(basemat.n_rows); 
  out.zeros();
  if (vertpl) { // Tall cases 
  #pragma omp parallel num_threads(num_threads)
  { 
  vec out_ = out; //local out 
  out_.zeros();
  vec temp_; //local temp vector 
  #pragma omp for nowait  
  for(uword lcv = 0; lcv < loopsize; lcv+=1) { 
    uword startind = lcv*chunksize; 
    uword endind = std::min((lcv+1)*chunksize-1,basemat.n_rows-1); 
    domult_(out_, a, temp_, terms, knotptst, basemat.rows(startind,endind), 
            num_threads);
    
    #pragma omp critical 
    out.subvec(startind,endind) = out_; //putting it back in the right place 
  }
  } 
  } else { // Wide cases 
    vec temp; //local temp vector 
    domult_(out, a, temp, terms, knotptst, basemat, num_threads);
  } 
  out %= basescale; 
}

/*
 * domultgesub_
 *
 * supporting function for domultge_
 */

void domultgesub_(vec& out, mat& outge, const vec& a, 
                  vec& temp, vec& tempalt,
                  const umat& terms, const uvec& knotptst, const mat& basemat,
                  const mat& basematge, const uvec& gest, const uvec& hypmatch,
                  const uword& k) {
  
  temp.fill(a[k]); 
  
  for (uword l = 0; l < terms.n_cols; ++l)  
    if(terms.at(k,l)>0) temp %= basemat.col(knotptst[l]+terms.at(k,l)); 
    out += temp;  
    for (uword l = 0; l < (gest.n_elem-1); ++l) { 
      if(terms.at(k,hypmatch[l])>0){  
        tempalt.fill(a[k]); 
        
        for (uword m = 0; m < terms.n_cols; ++m) 
          if(terms.at(k,m)>0 && m != hypmatch[l]) 
            tempalt %= basemat.col(knotptst[m]+terms.at(k,m)); 
          
          outge.col(l) += tempalt%basematge.col(gest[l]+ 
            terms.at(k,hypmatch[l]))-temp%basematge.col(gest[l]); 
      }
    }
}

/*
 * domultge_
 *
 * supporting function for prodmmge_
 */

void domultge_(vec& out, mat& outge, const vec& a, 
               vec& temp, vec& tempalt,
               const umat& terms, const uvec& knotptst, const mat& basemat,
               const mat& basematge, const uvec& gest, const uvec& hypmatch,
               int num_threads) {
  if(out.n_elem !=  basemat.n_rows) 
    out.set_size(basemat.n_rows); 
  if(outge.n_rows !=  basemat.n_rows || outge.n_cols !=  gest.n_elem-1) 
    outge.set_size(basemat.n_rows, gest.n_elem-1); 
  if(temp.n_elem !=  basemat.n_rows) 
    temp.set_size(basemat.n_rows); 
  if(tempalt.n_elem !=  basemat.n_rows) 
    tempalt.set_size(basemat.n_rows); 
  
  out.zeros();
  outge.zeros();
  
  uword termsnrows = terms.n_rows;
  if (omp_in_parallel()) {
    for (uword k = 0; k < termsnrows; ++k)
      domultgesub_(out, outge, a, temp, tempalt,
                   terms, knotptst, basemat,
                   basematge, gest, hypmatch,
                   k);
  } else {
  #pragma omp parallel num_threads(num_threads)
  {
  vec out_ = out; //local out 
  mat outge_ = outge; //local out 
  vec temp_ = temp; //local temp vector 
  vec tempalt_ = tempalt; //local temp vector 
  out_.zeros();
  outge_.zeros();
  #pragma omp for nowait  
  for (uword k = 0; k < termsnrows; ++k)
    domultgesub_(out_, outge_, a, temp_, tempalt_,
                 terms, knotptst, basemat,
                 basematge, gest, hypmatch,
                 k);
  #pragma omp critical  
  {
  out += out_;  
  outge += outge_;
  }
  }
  }
}

/*
 * prodmmge_
 *
 * gradient of prodmm_ with respect to the hyper parameters
 */

void prodmmge_(vec& out, mat& outge, 
               const umat& terms, 
               const vec& a,
               const mat& basemat, const vec& basescale, const uvec& knotptst,
               const mat& basematge, const uvec& gest, const uvec& hypmatch,
               bool vertpl, uword chunksize, uword loopsize, int num_threads) {
  if(out.n_elem !=  basemat.n_rows) 
    out.set_size(basemat.n_rows); 
  if(outge.n_rows !=  basemat.n_rows || outge.n_cols !=  gest.n_elem-1) 
    outge.set_size(basemat.n_rows, gest.n_elem-1); 
  out.zeros();
  outge.zeros();
  if (vertpl) { 
  #pragma omp parallel num_threads(num_threads)
  { 
  vec out_(chunksize); 
  mat outge_(chunksize, gest.n_elem-1); 
  out_.zeros();
  outge_.zeros();
  vec temp_; 
  vec tempalt_; 
  uword startind; 
  uword endind; 
  #pragma omp for nowait  
  for(uword lcv = 0; lcv < loopsize; lcv+=1){ 
    startind = lcv*chunksize; 
    endind = std::min((lcv+1)*chunksize-1, basemat.n_rows-1); 
    
    domultge_(out_, outge_, a, temp_, tempalt_,
              terms, knotptst, basemat.rows(startind,endind),
              basematge.rows(startind,endind), gest, hypmatch,
              num_threads);
    
  #pragma omp critical
  {
  out.subvec(startind,endind) = out_; 
  outge.rows(startind,endind) = outge_; 
  }
  } 
  } 
  } else { 
    vec temp_; 
    vec tempalt_; 
    domultge_(out, outge, a, temp_, tempalt_,
              terms, knotptst, basemat,
              basematge, gest, hypmatch,
              num_threads);
  } 
  for (uword l = 0; l < (gest.n_elem-1); ++l) 
    outge.col(l) += basematge.col(gest[l]) % out; 
  out %= basescale;
  outge.each_col() %= basescale;  
}


/*
 * dotmultsub_
 *
 * supporting function for tprodmm_
 */

void dotmultsub_(vec& out, vec& temp, 
                 const mat& basemat,
                 const uvec& knotptst, const umat& terms, const vec& b,
                 const uword& k) {
  temp = b; 
  for (uword l = 0; l < terms.n_cols; ++l) 
    if(terms.at(k,l)>0) temp %= basemat.col(knotptst[l]+terms.at(k,l)); 
  out(k) += sum(temp); 
}

/*
 * tprodmm_
 *
 * does transpose matrix multiplication in a memory friendly way that leverages 
 * omp.
 */

void tprodmm_(vec& out, 
              const umat& terms, 
              const vec& a,
              const mat& basemat, const vec& basescale, const uvec& knotptst,
              bool vertpl, uword chunksize, uword loopsize, int num_threads) {
  if(out.n_elem !=  terms.n_rows) 
    out.set_size(terms.n_rows); 
  
  out.zeros(); //zero it out just in case 
  vec b = basescale % a; 
  if(vertpl){ // Tall cases 
  #pragma omp parallel num_threads(num_threads)
  { 
  vec out_ = out; //local out 
  vec temp_; //local temp vector 
  mat basemat_; 
  vec b_;
  out_.zeros();
  #pragma omp for nowait  
  for(uword lcv = 0; lcv < loopsize; lcv+=1){
    uword startind = lcv*chunksize; 
    uword endind = std::min((lcv+1)*chunksize-1, basemat.n_rows-1); 
    basemat_ = basemat.rows(startind,endind); //extract rows 
    b_ = b.subvec(startind,endind);
    
    uword termsnrows = terms.n_rows;
    for (uword k = 0; k < termsnrows; ++k) 
      dotmultsub_(out_, temp_, basemat_,
                  knotptst, terms, b_,
                  k);
  }
  #pragma omp critical 
  out += out_; 
  } 
  } else { // Wide cases 
  #pragma omp parallel num_threads(num_threads)
  { 
  vec out_ = out; //local out 
  vec temp_; //local temp vector 
  out_.zeros();
  
  uword termsnrows = terms.n_rows;
  #pragma omp for nowait  
  for (uword k = 0; k < termsnrows; ++k) 
    dotmultsub_(out_, temp_, basemat,
                knotptst, terms, b,
                k);
  
  #pragma omp critical 
  out += out_; 
  }
  }
}


/*
 * dotmultsub_
 *
 * supporting function for dotmultgesub_
 */

void dotmultgesub_(vec& out, mat& outge,
                   vec& temp, vec& tempalt, 
                   const mat& basemat, const mat& basematge, 
                   const uvec& gest, const uvec& hypmatch, const uvec& knotptst, 
                   const umat& terms, const vec& b,
                   const uword& k) {
  temp = b;
  
  for (uword l = 0; l < terms.n_cols; ++l) 
    if(terms.at(k,l)>0) temp %= basemat.col(knotptst[l]+terms.at(k,l));
  out(k) += sum(temp);
  
  for (uword l = 0; l < outge.n_cols; ++l) { 
    if (terms.at(k,hypmatch[l]) > 0){ 
      tempalt = b;
      for (uword m = 0; m < terms.n_cols; ++m) 
        if(terms.at(k,m)>0 && m != hypmatch[l]) 
          tempalt %= basemat.col(knotptst[m]+terms.at(k,m)); 
      outge(k,l) += dot(tempalt, basematge.col(gest[l] +
        terms.at(k,hypmatch[l])));
    } else outge(k,l) += dot(temp, basematge.col(gest[l]));
  } 
}

/*
 * tprodmmge_
 *
 * gradient of tprodmm_ with respect to the hyper parameters
 */

void tprodmmge_(vec& out, mat& outge, 
                const umat& terms, 
                const vec& a,
                const mat& basemat, const vec& basescale, const uvec& knotptst,
                const mat& basematge, const uvec& gest, const uvec& hypmatch,
                bool vertpl, uword chunksize, uword loopsize, int num_threads) {
  if(out.n_elem !=  terms.n_rows) 
    out.set_size(terms.n_rows); 
  if(outge.n_rows !=  terms.n_rows || outge.n_cols !=  gest.n_elem-1)  
    outge.set_size(terms.n_rows, gest.n_elem-1);  
  
  vec b = basescale % a; 
  out.zeros();
  outge.zeros();
  if(vertpl) { // Tall cases 
  #pragma omp parallel num_threads(num_threads)
  {
  mat outge_ = outge;
  vec out_ = out;
  
  vec temp;
  vec tempalt;
  
  mat basemat_;
  mat basematge_;
  vec b_;
  out_.zeros();
  outge_.zeros();
  #pragma omp for nowait  
  for (uword lcv = 0; lcv < loopsize; lcv+=1) { 
    uword startind = lcv*chunksize; 
    uword endind = std::min((lcv+1)*chunksize-1,basemat.n_rows-1); 
    basemat_ = basemat.rows(startind,endind); //extract rows 
    basematge_ = basematge.rows(startind,endind); //extract rows 
    b_ = b.subvec(startind,endind);
    uword termsnrows = terms.n_rows;
    for (uword k = 0; k < termsnrows; ++k) 
      dotmultgesub_(out_, outge_,
                    temp, tempalt, 
                    basemat_, 
                    basematge_, 
                    gest, hypmatch, knotptst, 
                    terms, b_, k);
  }
  #pragma omp critical
  {
    out += out_;
    outge += outge_;
  }
  }
  } else {
  #pragma omp parallel num_threads(num_threads)
  { 
  vec temp;
  vec tempalt;
  
  mat outge_ = outge;
  vec out_ = out;
  
  out_.zeros();
  outge_.zeros();
  uword termsnrows = terms.n_rows;
  #pragma omp for nowait  
  for (uword k = 0; k < termsnrows; ++k) 
    dotmultgesub_(out_, outge_,
                  temp, tempalt, 
                  basemat, 
                  basematge, 
                  gest, hypmatch, knotptst, 
                  terms, b, k);
  #pragma omp critical
  {
  out += out_;
  outge += outge_;
  }
  }
  }
}



/*
 * domult_
 *
 * supporting function for prodmm_
 */

void domultm_(mat& out, const mat& a, vec& temp, 
              const umat& terms, const uvec& knotptst, const mat& basemat,
              int num_threads) {
  if(out.n_elem !=  basemat.n_rows || out.n_cols !=  a.n_cols) 
    out.set_size(basemat.n_rows, a.n_cols); 
  if(temp.n_elem !=  basemat.n_rows) 
    temp.set_size(basemat.n_rows); 
  
  out.zeros();
  uword termsnrows = terms.n_rows;
  if (omp_in_parallel()) {
    for (uword k = 0; k < termsnrows; ++k) { 
      temp.ones(); 
      for (uword l = 0; l < terms.n_cols; ++l)  
        if(terms.at(k,l)>0) temp %= basemat.col(knotptst[l]+terms.at(k,l)); 
      out += temp * a.row(k);  
    }
  } else {
  #pragma omp parallel num_threads(num_threads)
  {
  mat out_ = out;
  vec temp_ = temp;
  if(out_.n_elem !=  basemat.n_rows || out_.n_cols !=  a.n_cols) 
    out_.set_size(basemat.n_rows, a.n_cols); 
  if(temp_.n_elem !=  basemat.n_rows) 
    temp_.set_size(basemat.n_rows); 
  out_.zeros();
  #pragma omp for
  for (uword k = 0; k < termsnrows; ++k) { 
    temp_.ones(); 
    for (uword l = 0; l < terms.n_cols; ++l)  
      if(terms.at(k,l)>0) temp_ %= basemat.col(knotptst[l]+terms.at(k,l)); 
    out_ += temp_ * a.row(k);  
  } 
  #pragma omp critical  
  out += out_;  
  }
  }
}

/*
 * prodmm_
 *
 * matrix version of prodmm_
 */

void prodmm_(mat& out, 
             const umat& terms, 
             const mat& a,
             const mat& basemat, const vec& basescale, const uvec& knotptst,
             bool vertpl, uword chunksize, uword loopsize, int num_threads) {
  if(out.n_elem !=  basemat.n_rows || out.n_cols !=  a.n_cols) 
    out.set_size(basemat.n_rows, a.n_cols); 
  out.zeros();
  if (vertpl) { // Tall cases 
  #pragma omp parallel num_threads(num_threads)
  { 
  mat out_ = out; //local out 
  out_.zeros();
  vec temp_; //local temp vector 
  #pragma omp for nowait  
  for(uword lcv = 0; lcv < loopsize; lcv+=1) { 
    uword startind = lcv*chunksize; 
    uword endind = std::min((lcv+1)*chunksize-1,basemat.n_rows-1); 
    domultm_(out_, a, temp_, terms, knotptst, basemat.rows(startind,endind), 
             num_threads);
    
    #pragma omp critical 
    out.rows(startind,endind) = out_; //putting it back in the right place 
  } 
  } 
  } else { // Wide cases 
    vec temp; //local temp vector 
    domultm_(out, a, temp, terms, knotptst, basemat, num_threads);
  } 
  out.each_col() %= basescale; 
}



/*
 * dotmultsub_
 *
 * supporting function for tprodmm_ 
 */

void dotmmultsub_(mat& out, vec& temp, 
                  const mat& basemat,
                  const uvec& knotptst, const umat& terms, const mat& b,
                  const uword& k) {
  temp.ones(); 
  for (uword l = 0; l < terms.n_cols; ++l) 
    if(terms.at(k,l)>0) temp %= basemat.col(knotptst[l]+terms.at(k,l)); 
  out.row(k) += temp.t() * b; 
}

/*
 * tprodmm_
 *
 * matrix version of tprodmm_
 */

void tprodmm_(mat& out, 
              const umat& terms, 
              const mat& a,
              const mat& basemat, const vec& basescale, const uvec& knotptst,
              bool vertpl, uword chunksize, uword loopsize, int num_threads) {
  
  if(out.n_elem !=  terms.n_rows || out.n_cols !=  a.n_cols) 
    out.set_size(terms.n_rows, a.n_cols); 
  
  out.zeros(); //zero it out just in case 
  mat b = a; 
  b.each_col() %= basescale;
  if(vertpl){ // Tall cases 
  #pragma omp parallel num_threads(num_threads)
  { 
  mat out_ = out; //local out 
  vec temp_; //local temp vector 
  mat basemat_; 
  mat b_;
  out_.zeros();
  #pragma omp for nowait  
  for(uword lcv = 0; lcv < loopsize; lcv+=1){
    uword startind = lcv*chunksize; 
    uword endind = std::min((lcv+1)*chunksize-1, basemat.n_rows-1); 
    basemat_ = basemat.rows(startind,endind); //extract rows 
    b_ = b.rows(startind,endind);
    temp_.set_size(b_.n_rows);
    uword termsnrows = terms.n_rows;
    for (uword k = 0; k < termsnrows; ++k) 
      dotmmultsub_(out_, temp_, basemat_,
                   knotptst, terms, b_,
                   k);
  }
  #pragma omp critical 
  out += out_; 
  } 
  } else { // Wide cases 
  #pragma omp parallel num_threads(num_threads)
  { 
  mat out_ = out; //local out 
  vec temp_ ; //local temp vector 
  out_.zeros();
  temp_.set_size(b.n_rows);
  
  #pragma omp for nowait  
  for (uword k = 0; k < terms.n_rows; ++k) 
    dotmmultsub_(out_, temp_, basemat,
                 knotptst, terms, b,
                 k);
  
  #pragma omp critical 
  out += out_; 
  }
  }
}



/*
 * domat_
 *
 * supporting function for prodmm_
 */

void domat_(mat& out, vec& temp, 
            const umat& terms, const uvec& knotptst, const mat& basemat,
            int num_threads) {
  if(out.n_elem !=  basemat.n_rows || out.n_cols !=  terms.n_rows) 
    out.set_size(basemat.n_rows, terms.n_rows); 
  if(temp.n_elem !=  basemat.n_rows) temp.set_size(basemat.n_rows); 
  if (omp_in_parallel()) {
    uword termsnrows = terms.n_rows;
    for (uword k = 0; k < termsnrows; ++k) { 
      temp.ones(); 
      for (uword l = 0; l < terms.n_cols; ++l)  
        if(terms.at(k,l)>0) temp %= basemat.col(knotptst[l]+terms.at(k,l)); 
      out.col(k) = temp;  
    }
  } else {
  #pragma omp parallel num_threads(num_threads)
  {
  vec temp_ = temp;
  if(temp_.n_elem !=  basemat.n_rows) temp_.set_size(basemat.n_rows); 
  uword termsnrows = terms.n_rows;
  #pragma omp for nowait
  for (uword k = 0; k < termsnrows; ++k) { 
    temp_.ones(); 
    for (uword l = 0; l < terms.n_cols; ++l)  
      if(terms.at(k,l)>0) temp_ %= basemat.col(knotptst[l]+terms.at(k,l)); 
    #pragma omp critical  
    out.col(k) = temp_;  
  }
  }
  }
}

/*
 * prodmm_
 *
 * matrix version of prodmm_
 */

void getm_(mat& out, 
           const umat& terms, 
           const mat& basemat, const vec& basescale, const uvec& knotptst,
           bool vertpl, uword chunksize, uword loopsize, int num_threads) {
  if(out.n_elem !=  basemat.n_rows || out.n_cols !=  terms.n_rows) 
    out.set_size(basemat.n_rows, terms.n_rows); 
  out.zeros();
  if (vertpl) { // Tall cases 
  #pragma omp parallel num_threads(num_threads)
  { 
  
  mat out_ = out; //local out 
  out_.zeros();
  vec temp_; //local temp vector 
  #pragma omp for nowait  
  for(uword lcv = 0; lcv < loopsize; lcv+=1) { 
    uword startind = lcv*chunksize; 
    uword endind = std::min((lcv+1)*chunksize-1,basemat.n_rows-1); 
    domat_(out_, temp_, terms, knotptst, basemat.rows(startind,endind), 
           num_threads);
    
    #pragma omp critical 
    out.rows(startind,endind) = out_; //putting it back in the right place 
  } 
  } 
  } else { // Wide cases 
    vec temp; //local temp vector 
    domat_(out, temp, terms, knotptst, basemat, num_threads);
  } 
  out.each_col() %= basescale; 
}


/*
 * domultge_
 *
 * supporting function for prodmmge_
 */

void dogetmge_(cube& outge, 
               vec& temp, vec& tempalt,
               const umat& terms, const uvec& knotptst, const mat& basemat,
               const mat& basematge, const uvec& gest, const uvec& hypmatch,
               int num_threads) {
  if(outge.n_rows !=  basemat.n_rows || 
     outge.n_cols !=  terms.n_rows ||
     outge.n_slices !=  gest.n_elem-1) 
    outge.set_size(basemat.n_rows,terms.n_rows, gest.n_elem-1); 
  if(temp.n_elem !=  basemat.n_rows) 
    temp.set_size(basemat.n_rows); 
  if(tempalt.n_elem !=  basemat.n_rows) 
    tempalt.set_size(basemat.n_rows); 
  
  outge.zeros();
  
  if (omp_in_parallel()) {
    for (uword k = 0; k < terms.n_rows; ++k){
      for (uword l = 0; l < (gest.n_elem-1); ++l) { 
        tempalt.ones(); 
        for (uword m = 0; m < terms.n_cols; ++m)
          if(terms.at(k,m)>0 && m != hypmatch[l])
            tempalt %= basemat.col(knotptst[m]+terms.at(k,m)); 
        outge.slice(l).col(k) = tempalt % basematge.col(gest[l] + 
          terms.at(k,hypmatch[l])); 
      }
    }
  } else {
  #pragma omp parallel num_threads(num_threads)
  {
  vec temp_ = temp; //local temp vector 
  vec tempalt_ = tempalt; //local temp vector 
  #pragma omp for nowait  
  for (uword k = 0; k < terms.n_rows; ++k){
    for (uword l = 0; l < (gest.n_elem-1); ++l) { 
      tempalt_.ones(); 
      for (uword m = 0; m < terms.n_cols; ++m)
        if(terms.at(k,m)>0 && m != hypmatch[l])
          tempalt_ %= basemat.col(knotptst[m]+terms.at(k,m)); 
      #pragma omp critical  
      outge.slice(l).col(k) = tempalt_%basematge.col(gest[l]+
        terms.at(k,hypmatch[l])); 
    }
  }
  }
  }
}

/*
 * prodmmge_
 *
 * gradient of prodmm_ with respect to the hyper parameters
 */

void getmge_(cube& outge, 
             const umat& terms, 
             const mat& basemat, const vec& basescale, const uvec& knotptst,
             const mat& basematge, const uvec& gest, const uvec& hypmatch,
             bool vertpl, uword chunksize, uword loopsize, int num_threads) {
  if(outge.n_rows !=  basemat.n_rows || 
     outge.n_cols !=  terms.n_rows ||
     outge.n_slices !=  gest.n_elem-1) 
    outge.set_size(basemat.n_rows,terms.n_rows, gest.n_elem-1); 
  outge.zeros();
  if (vertpl) { 
  #pragma omp parallel num_threads(num_threads)
  { 
  vec out_(chunksize); 
  mat outge_(chunksize, gest.n_elem-1); 
  out_.zeros();
  outge_.zeros();
  vec temp_; 
  vec tempalt_; 
  uword startind; 
  uword endind; 
  #pragma omp for nowait  
  for(uword lcv = 0; lcv < loopsize; lcv+=1){ 
    startind = lcv*chunksize; 
    endind = std::min((lcv+1)*chunksize-1, basemat.n_rows-1); 
    dogetmge_(outge, temp_, tempalt_,
              terms, knotptst, basemat.rows(startind,endind),
              basematge.rows(startind,endind), gest, hypmatch,
              num_threads);
    #pragma omp critical 
    outge.rows(startind,endind) = outge_; 
  }
  } 
  } else { 
    vec temp_; 
    vec tempalt_; 
    dogetmge_(outge, temp_, tempalt_,
              terms, knotptst, basemat,
              basematge, gest, hypmatch,
              num_threads);
  } 
  for(uword k = 0; k < outge.n_slices; k+=1){ 
    outge.slice(k).each_col() %= basescale;  
  }
}
