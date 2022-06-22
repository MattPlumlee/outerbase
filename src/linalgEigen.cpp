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
#include <RcppEigen.h>


// No omp specification needed, let Eigen handle it


using namespace Rcpp;
using namespace arma;

#include "linalgEigen.h"

/*
 * prodmmE_
 *
 * does matrix multiplication in a memory friendly way that leverages Eigen.
 */

void prodmmE_(vec& out, 
             const umat& terms, 
             const vec& a,  
             const mat& basemat, const vec& basescale, const uvec& knotptst,
             bool vertpl, uword chunksize, uword loopsize, int num_threads) {
  Eigen::setNbThreads(4);  
  
  if(out.n_elem !=  basemat.n_rows) out.set_size(basemat.n_rows);   
    
  const Eigen::ArrayXXd eigen_B = Eigen::Map<const Eigen::ArrayXXd>(
    basemat.memptr(),basemat.n_rows,basemat.n_cols);
  const Eigen::ArrayXXf eigen_Bf = eigen_B.cast <float> ();
  Eigen::ArrayXf eigen_outf(basemat.n_rows);
  Eigen::ArrayXf eigen_tempf(basemat.n_rows);
  
  const Eigen::ArrayXd eigen_a = Eigen::Map<const Eigen::ArrayXd>(
    a.memptr(),a.n_elem);
  const Eigen::ArrayXXf eigen_af = eigen_a.cast <float> ();
  
  eigen_outf.setZero();
  unsigned int termsnrows = terms.n_rows;
  unsigned int termsncols = terms.n_cols;
  for (unsigned int k = 0; k < termsnrows; ++k) { 
    eigen_tempf.setConstant(eigen_af(k));
    for (unsigned int l = 0; l < termsncols; ++l)  
      if(terms.at(k,l)>0) eigen_tempf *= eigen_Bf.col(knotptst[l]+terms.at(k,l)); 
    eigen_outf += eigen_tempf;
  }
  
  Eigen::ArrayXd eigen_out = eigen_outf.cast <double> ();
  out = arma::vec(eigen_out.data(), basemat.n_rows, false, false);
  
  out %= basescale; 
  
}

