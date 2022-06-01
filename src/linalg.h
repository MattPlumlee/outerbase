#ifndef __LINALG__
#define __LINALG__

#ifndef __COVFUNCS__
#define __COVFUNCS__
#include "covfuncs.h"
#endif

void prodmm_(vec& out, 
             const umat& terms, 
             const vec& a,
             const mat& basemat, const vec& basescale, const uvec& knotptst,
             bool vertpl, uword chunksize, uword loopsize, int num_threads) ;

void prodmm_(mat& out, 
             const umat& terms, 
             const mat& a,
             const mat& basemat, const vec& basescale, const uvec& knotptst,
             bool vertpl, uword chunksize, uword loopsize, int num_threads) ;

void tprodmm_(vec& out, 
              const umat& terms, 
              const vec& a,
              const mat& basemat, const vec& basescale, const uvec& knotptst,
              bool vertpl, uword chunksize, uword loopsize, int num_threads) ;


void tprodmm_(mat& out, 
              const umat& terms, 
              const mat& a,
              const mat& basemat, const vec& basescale, const uvec& knotptst,
              bool vertpl, uword chunksize, uword loopsize, int num_threads) ;

void prodmmge_(vec& out, mat& outge, 
               const umat& terms, 
               const vec& a,
               const mat& basemat, const vec& basescale, const uvec& knotptst,
               const mat& basematge, const uvec& gest, const uvec& hypmatch,
               bool vertpl, uword chunksize, uword loopsize, int num_threads) ;

void tprodmmge_(vec& out, mat& outge, 
                  const umat& terms, 
                  const vec& a,
                  const mat& basemat, const vec& basescale, const uvec& knotptst,
                  const mat& basematge, const uvec& gest, const uvec& hypmatch,
                  bool vertpl, uword chunksize, uword loopsize, int num_threads) ;


void getm_(mat& out,
            const umat& terms, 
            const mat& basemat, const vec& basescale, const uvec& knotptst,
            bool vertpl, uword chunksize, uword loopsize, int num_threads) ;

void getmge_(cube& out,
           const umat& terms, 
           const mat& basemat, const vec& basescale, const uvec& knotptst,
           const mat& basematge, const uvec& gest, const uvec& hypmatch,
           bool vertpl, uword chunksize, uword loopsize, int num_threads) ;
  
#endif