#ifndef __COVFUNCS__
#define __COVFUNCS__

class covf {
public:
  vec hyp;
  vec hypub;
  vec hyplb;
  vec hyp0;
  vec hypvar;
  double lowbnd;
  double uppbnd;
  
  unsigned int numhyp;
  std::vector<std::string> hypnames;
  
  covf(){}
  virtual ~covf() {}
  
  double lpdf(vec) const;
  vec lpdf_gradhyp(vec) const;
  
  virtual bool inputcheck(const vec& x) const {
    if (x.min() < lowbnd) return false;
    if (x.max() > uppbnd) return false;
    return true;
    }
  
  virtual mat cov(const vec&, const vec&) const {return {};}
  virtual vec covmdiag(const vec&) const {return {};}
  virtual cube cov_gradhyp(const vec&, const vec&) const {return {};}
};

class covf_mat25: public covf {
public:
  covf_mat25();
  
  mat cov(const vec&, const vec&) const;
  vec covmdiag(const vec&) const;
  cube cov_gradhyp(const vec&, const vec&) const;
private:
  double a = 2.;
};

class covf_mat25pow: public covf {
public:
  covf_mat25pow();
  
  mat cov(const vec&, const vec&) const;
  vec covmdiag(const vec&) const;
  cube cov_gradhyp(const vec&, const vec&) const;
private:
  double a = 2.;
  double b = 0.25;
};


class covf_mat25ang: public covf {
public:
  covf_mat25ang();
  
  mat cov(const vec&, const vec&) const;
  vec covmdiag(const vec&) const;
  cube cov_gradhyp(const vec&, const vec&) const;
private:
  double a = 2.;
};


#endif