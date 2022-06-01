#ifndef __MODANDBASE__
#define __MODANDBASE__

#ifndef __COVFUNCS__
#define __COVFUNCS__
#include "covfuncs.h"
#endif

class outermod {
public:
  uword d; //dim of input
  vec basisvar; //measures the log-size of each basis
  
  bool setcovfs = false;
  std::vector<covf*> covflist;
  
  vec hyp; //vector of hyperparameters for the basis
  
  uvec knotptst; //where knot points start on the vector for each dim
  uvec hypmatch; //which dim does each parameter belong to
  uvec hypst; //where the vectors start for each dim
  uvec gest; //
  uvec knotptstge; //eqv. of rasmat for rasmatge
  vec knotpt; //vector of knot points for all dims
  
  bool setknots = false;
  
  ivec maxlevel; //maximum allowable level for each basis
  
  outermod(){}
  void buildob(mat&, const mat&, const uword&) const;
  void buildob(mat&, cube&, const mat&, const uword&) const;
  vec totvar(const mat&) const;
  
  vec getvar(const umat& terms) const;
  mat getlvar_gradhyp(const umat& terms) const;
  umat selectterms(const unsigned int) const;
  
  double hyplpdf(vec) const;
  vec hyplpdf_grad(vec) const;
  
  void hyp_init();
  void hyp_set(vec);
  
  void build();
private:
  mat rotmat; //modifier matrix to make nearly orthogonal
  
  mat rotmat_gradhyp; // gives the gradient of rasmat from the hyperparameters
  vec logbasisvar_gradhyp; // gives the gradient of basisvar from the 
  // hyperparameters
  
  void setsizes_();
};


class outerbase {
public:
  const outermod& om;
  const mat xp;
  mat basemat;
  
  uword d;
  uword n_row;
  uword n_hyp;
  bool dograd;
  uvec hypst;
  
  uword nthreads;//used for omp
  bool vertpl = false;//used for omp
  
  outerbase(const outermod&, mat);
  outerbase(const outermod&, mat, bool);
  void build();
  
  mat getbase(const uword) const;
  mat getmat(const umat& terms) const;  
  cube getmat_gradhyp(const umat& terms) const;  
  
  void mm(vec& out, const umat& terms, const vec& a) const; 
  void tmm(vec& out, const umat& terms, const vec& a) const; 
  void mm_gradhyp(vec& out, mat& out_gradhyp, 
                  const umat& terms, const vec& a) const;  
  void tmm_gradhyp(vec& out, mat& out_gradhyp, 
                   const umat& terms, const vec& a) const;  
  
  
  vec mm_out(const umat& terms, const vec& a) const; 
  vec tmm_out(const umat& terms, const vec& a) const; 
  mat mm_gradhyp_out(const umat& terms, const vec& a) const;  
  mat tmm_gradhyp_out(const umat& terms, const vec& a) const;
  
  vec sqcolsums(const umat& terms) const;
  mat sqcolsums_gradhyp(const umat& terms) const;
  
  vec sqmm(const umat& terms, const vec& a) const;
  mat sqmm_gradhyp(const umat& terms, const vec& a) const;
  
  vec sqtmm(const umat& terms, const vec& a) const;
  mat sqtmmm(const umat& terms, const mat& a) const;
  mat sqtmm_gradhyp(const umat& terms, const vec& a) const;
  
  vec residvar(const umat& terms) const;
  mat residvar_gradhyp(const umat& terms) const;
private:
  uword n_base;
  uvec knotptst;
  vec basescale;
  mat basescalemat;
  mat basemat_gradhyp;
  mat basematsq;
  vec basescalesq;
  mat basematsq_gradhyp;
  
  uvec gest;
  uvec hypmatch;
  
  uword loopsize; //used for omp
  uword chunksize = 200; //used for omp
  
  void setvals_();
  void setsizes_();
  mat sqbasemat_gradhyp() const;
};

#endif


