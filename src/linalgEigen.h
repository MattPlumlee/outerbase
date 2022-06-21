#ifndef __LINALGEIGEN__
#define __LINALGEIGEN__

void prodmmE_(vec& out, 
             const umat& terms, 
             const vec& a,
             const mat& basemat, const vec& basescale, const uvec& knotptst,
             bool vertpl, uword chunksize, uword loopsize, int num_threads) ;
  
#endif