


/* 
 **************************************************************************** 
 **************************************************************************** 
 *************************** loglik_gauss *********************************** 
 **************************************************************************** 
 **************************************************************************** 
 */ 

/*
 * loglik_gauss::loglik_gauss
 *
 */

loglik_gauss::loglik_gauss(const outermod& om_, umat terms_, vec y_, mat x_)
  :  lpdf(), om(om_), ob(om,x_,true),  y(y_), x(x_)
{
  terms = terms_;
  npara = 1;
  
  para0.set_size(1);
  para0[0] = log(0.01*var(y));
  paravar.set_size(1);
  paravar[0] = 1;
  
  paranames = {"noisescale"};
  para = para0;
  obssd.resize(y.n_elem);
  obssd.fill(exp(para(0)));
  obsvar.resize(y.n_elem);
  obsvar.fill(exp(2.*para(0)));
  nterms = terms.n_rows;
  lobsvar = log(obsvar);
}


/*
 * loglik_gauss::updateom
 *
 */

void loglik_gauss::updateom() {
  ob.build();
}

/* 
 * loglik_gda::setnthreads
 * 
 */ 

void loglik_gauss::setnthreads(int k) {
  ob.nthreads = k;
}


/*
 * loglik_gauss::updatepara
 *
 */

void loglik_gauss::updatepara(vec para_) {
  para = para_;
  obssd.fill(exp(para(0)));
  obsvar.fill(exp(2.*para(0)));
  lobsvar = log(obsvar);
}

/*
 * loglik_gauss::updateterms
 *
 */

void loglik_gauss::updateterms(umat terms_) {
  terms = terms_;
  nterms = terms.n_rows;
}


/*
 * loglik_gauss::update
 *
 */

void loglik_gauss::update(const vec& coeff_) {
  coeff = coeff_;
  yhat.copy_size(y);
  residtemp.copy_size(y);
  
  if(compute_gradhyp){
    yhatge.set_size(y.n_elem, ob.n_hyp);
    ob.mm_gradhyp(yhat, yhatge, terms, coeff);
  } else ob.mm(yhat, terms, coeff);
  
  residtemp = (yhat-y)/obssd;
  residtemp2 = square(residtemp);
  if(compute_val) val = -0.5*sum(square(residtemp)) - sum(log(obssd));
  if(compute_grad){
    grad.copy_size(coeff);
    residtemp = -1.*(residtemp/obssd);
    ob.tmm(grad, terms, residtemp);
    if(compute_gradhyp) gradhyp = (residtemp.t() * yhatge).as_col();
    if(compute_gradpara) gradpara = sum(residtemp2) - y.n_elem;
  }
}

/*
 * loglik_gauss::hessmult
 *
 */

vec loglik_gauss::hessmult(const vec& g) {
  yhattemp.copy_size(y);
  ob.mm(yhattemp, terms, g);
  yhattemp /= obssd;
  yhattemp /= obssd;
  gradtemp.copy_size(grad);
  ob.tmm(gradtemp, terms, yhattemp);
  return gradtemp;
}


/*
 * loglik_gauss::diaghess
 *
 */


vec loglik_gauss::diaghess() {
  vec lh = ob.sqcolsums(terms);
  return ((exp(-2*para(0)))*lh);
}


/*
 * loglik_gauss::diaghessgradhyp
 *
 */

mat loglik_gauss::diaghessgradhyp() {
  mat h = ob.sqcolsums_gradhyp(terms);
  return ((exp(-2*para(0)))*h);
}


/*
 * loglik_gauss::diaghessgradpara
 *
 */

mat loglik_gauss::diaghessgradpara() {
  vec lh = ob.sqcolsums(terms);
  return ((-2*exp(-2*para(0)))*lh);
}

predf* loglik_gauss::pred() const {
  return (new pred_gauss(*this));
}


/* 
 **************************************************************************** 
 **************************************************************************** 
 *************************** pred_gauss ************************************* 
 **************************************************************************** 
 **************************************************************************** 
 */ 



pred_gauss::pred_gauss(const loglik_gauss& loglik) : 
om(loglik.om), para(loglik.para), terms(loglik.terms),
x(loglik.x), ob(om,loglik.x,false) //private
{
  coeff = loglik.coeff;
  
  if(!loglik.didnotothess) {
    if(loglik.didfulltothess) coeffvar = 1/loglik.totdiaghess.diag();
    else {
      coeffvar = 1/loglik.totdiaghess;
    }
  } else {
    coeffvar = 0*coeff;
  }
}

void pred_gauss::update(const mat& x_) {
  x = x_;
  new (&ob) outerbase(om,x_,false);
}

vec pred_gauss::mean() const {
  return ob.mm_out(terms, coeff);
}
vec pred_gauss::var() const {
  vec out = ob.sqmm(terms, coeffvar);
  out += exp(2*para(0));
  return out;
}