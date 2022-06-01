#ifndef __FIT__
#define __FIT__

#ifndef __COVFUNCS__
#define __COVFUNCS__
#include "covfuncs.h"
#endif


class predf {
public:
  
  predf(){}
  virtual ~predf() {}
  virtual void update(const mat&){}
  virtual vec mean() const { return {};}
  virtual vec var() const { return {};}
};  


class lpdf {
public:
  double val;
  vec grad;
  vec gradhyp;
  vec gradpara;
  vec para;
  
  umat terms;
  vec coeff;
  vec totdiaghess;
  mat tothess;
  
  bool didfulltothess = false;
  bool didnotothess = true;
  bool fullhess = false;
  
  std::vector<std::string> paranames;
  
  bool compute_val = true;
  bool compute_grad = true;
  bool compute_gradhyp = false;
  bool compute_gradpara = false;
  unsigned int npara = 0;
  unsigned int nterms = 0;
  
  lpdf(){}
  virtual ~lpdf() {}
  
  virtual predf* pred() const {
    throw std::invalid_argument("cannot produce a predictor from this obj.");
    return {};
  }
  
  virtual void setnthreads(int) { return;}
  
  vec para0;
  vec paravar;
  virtual double paralpdf(vec) const;
  virtual vec paralpdf_grad(vec) const;
  
  virtual void optcg(double, unsigned int);
  virtual void optnewton();
  
  virtual void updateom(){}
  virtual void updatepara(vec){}
  virtual void updateterms(umat){}
  virtual void update(const vec&){}
  
  virtual vec hessmult(const vec&) {return {};}
  
  virtual vec diaghess() {return {};}
  virtual mat diaghessgradhyp() {return {};}
  virtual mat diaghessgradpara() {return {};}
  virtual void settotdiaghess(vec diaghess) {
    totdiaghess = diaghess;
    didfulltothess = false;
    didnotothess = false;
  }
  virtual void settothess(mat hess) {
    tothess = hess;
    didfulltothess = true;
    didnotothess = false;
  }
  virtual mat hess() {return {};}
  virtual cube hessgradhyp() {return {};}
  virtual cube hessgradpara() {return {};}
};


class lpdfvec: public lpdf {
public:
  double val_margadj;
  vec gradhyp_margadj;
  vec gradpara_margadj;
  bool domargadj = true;
  
  lpdfvec(lpdf&, lpdf&);
  
  void setnthreads(int);
  
  void updateom();
  void updatepara(vec);
  void updateterms(umat);
  void update(const vec&);
  
  vec hessmult(const vec&);
  
  vec diaghess();
  mat diaghessgradhyp();
  mat diaghessgradpara();
  void settotdiaghess(vec diaghess);
  void settothess(mat hess);
  
  mat hess();
  cube hessgradhyp();
  cube hessgradpara();
  
  double paralpdf(vec) const;
  vec paralpdf_grad(vec) const;
private:
  vec diaghessv;
  mat diaghessgradhypv;
  mat diaghessgradparav;
  
  mat hessv;
  cube hessgradhypv;
  cube hessgradparav;
  
  bool redohess = true;
  std::vector<std::reference_wrapper<lpdf>> lpdflist;
  
  uvec parasrt;
  uvec paraend;
  
  void buildhess();
  void margadj();
  
  vec diaghess_();
  mat diaghessgradhyp_();
  mat diaghessgradpara_();
  
  mat hess_();
  cube hessgradhyp_();
  cube hessgradpara_();
};


class logpr_gauss: public lpdf {
public:
  const outermod& om;
  
  vec coeffsd;
  
  logpr_gauss(const outermod&, umat);
  
  
  void updateom();
  void updatepara(vec);
  void updateterms(umat);
  void update(const vec& coeff);
  
  vec hessmult(const vec& g);
  
  vec diaghess();
  mat diaghessgradhyp();
  mat diaghessgradpara();
  
  mat hess();
  cube hessgradhyp();
  cube hessgradpara();
private:
  mat coefflvarge;
  mat coefflvargetemp;
  vec stdresid;
  vec residtemp;
  double sca;
};

class loglik_std: public lpdf {
public:
  const outermod& om;
  outerbase ob;
  mat basismat;
  cube basismat_gradhyp;
  
  vec y; 
  vec yhat;
  mat x;
  
  loglik_std(const outermod&, umat, vec, mat);
  
  void updateom();
  void updatepara(vec);
  void updateterms(umat);
  void update(const vec&);
  
  vec hessmult(const vec&);
  
  vec diaghess();
  mat diaghessgradhyp();
  mat diaghessgradpara();
  
  mat hess();
  cube hessgradhyp();
  cube hessgradpara();
  
  predf* pred() const; 
};


class predr_std: public predf {
public:
  const outermod& om;
  const vec para;
  const umat terms;
  mat basismat; 
  
  mat x;
  vec coeff;
  outerbase ob;  
  mat coeffcov;
  bool domargadj = true;
  
  predr_std(const loglik_std&);
  void update(const mat&);
  vec mean() const;
  vec var() const;
};

class loglik_gauss: public lpdf {
public:
  const outermod& om;
  outerbase ob;
  vec y;
  mat x;
  
  loglik_gauss(const outermod&, umat, vec, mat);
  
  void setnthreads(int);
  void updateom();
  void updatepara(vec);
  void updateterms(umat);
  void update(const vec& coeff);
  
  vec hessmult(const vec&);
  
  vec diaghess();
  mat diaghessgradhyp();
  mat diaghessgradpara();
  
  predf* pred() const; 
private:
  vec yhat;
  vec obsvar;
  vec lobsvar;
  vec obssd;
  mat yhatge;
  vec gradtemp;
  mat yhatgetemp;
  mat residtge;
  vec yhattemp;
  vec residtemp;
  vec residtemp2;
};

class pred_gauss : public predf {
public:
  const outermod& om;
  const vec para;
  const umat terms;
  
  mat x;
  vec coeff;
  outerbase ob;
  vec coeffvar;
  bool domargadj = true;
  
  pred_gauss(const loglik_gauss&);
  void update(const mat&);
  vec mean() const;
  vec var() const;
};

class loglik_gda: public lpdf {
public:
  const outermod& om;
  outerbase ob;
  vec y;
  mat x;
  
  bool doda;   
  
  loglik_gda(const outermod&, umat, vec, mat);
  
  void setnthreads(int);
  void updateom();
  void updatepara(vec);
  void updateterms(umat);
  void update(const vec&);
  
  vec hessmult(const vec&);
  
  vec diaghess();
  mat diaghessgradhyp();
  mat diaghessgradpara();
  
  predf* pred() const; 
private:
  vec yhat;
  
  bool redostd;
  vec obssd;
  mat obssd_gradpara;
  mat obssd_gradhyp;
  mat yhatge;
  vec gradtemp;
  mat yhatgetemp;
  mat residtge;
  vec yhattemp;
  vec residtemp;
  vec residtemp2;
  
  void buildstd();
};

class pred_gda : public predf {
public:
  const outermod& om;
  const vec para;
  const umat terms;
  
  mat x;
  vec coeff;
  outerbase ob;
  bool doda;
  bool domargadj = true;
  vec coeffvar;
  
  pred_gda(const loglik_gda&);
  void update(const mat& x_);
  vec mean() const;
  vec var() const;
};


class predictor {
public:
  predf* pred;
  
  predictor(const lpdf& logpdf) {pred = logpdf.pred();}
  void update(const mat& x){pred->update(x);}
  vec mean() const {return pred->mean();}
  vec var() const {return pred->var();}
};

#endif