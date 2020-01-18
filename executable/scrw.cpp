// switching correlated random walk

#include <TMB.hpp>
#include <iostream>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(y); // observed location (col(0): easting, col(1): northing)
  DATA_IVECTOR(obs_t); // observed time steps
  DATA_INTEGER(timeSteps); // maximum of time steps
  DATA_SCALAR(minDF); // observed time steps
  DATA_INTEGER(mode); //0: without mode, 1: with mode
  DATA_VECTOR(gamma_fix); //-1: estimate, 0~1: fixed at the value
  DATA_VECTOR(theta_fix); //-1: estimate otherwise: fixed at the value
  DATA_INTEGER(obs_error_type); //0: student-t, 1: normal
  DATA_IVECTOR(sigma_key); // c(0,0) or c(0,1)
  DATA_INTEGER(mixture); // 0 or 1
  
  PARAMETER_VECTOR(logit_alpha);
  PARAMETER(log_phi); // process error for mode transition
  PARAMETER_VECTOR(logit_gamma);
  PARAMETER_VECTOR(trans_theta); // mean angle for each mode
  PARAMETER_VECTOR(log_sigma); // process error for longitude and latitude
  PARAMETER(trans_rho);
  PARAMETER(trans_df);
  PARAMETER_VECTOR(logit_p); // probability of migration mode (latent variable)
  PARAMETER_MATRIX(x); // true location (latent variable)
  
  // int timeSteps = max(obs_t)+1;
  vector<Type> alpha = Type(1.0)/(Type(1.0)+exp(-logit_alpha));
  vector<Type> p = Type(1.0)/(Type(1.0)+exp(-logit_p));
  Type phi = exp(log_phi);
  vector<Type> gamma = Type(1.0)/(Type(1.0)+exp(-logit_gamma));
  for (int i=0; i<gamma.size(); i++) {
    if (gamma_fix(i) > Type(-1.0)) gamma(i) = gamma_fix(i);
  }
  vector<Type> sigma = exp(log_sigma);
  Type rho = (exp(trans_rho)-Type(1.0))/(exp(trans_rho)+Type(1.0)); // -1<rho<1
  Type pi = Type(3.141592653589);
  vector<Type> theta = pi*(exp(trans_theta)-Type(1.0))/(exp(trans_theta)+Type(1.0));
  for (int i=0; i<theta.size(); i++) {
    if (theta_fix(i) > Type(-1.0)) theta(i) = theta_fix(i);
  }
  Type df = minDF + exp(trans_df);
  Type omega = exp(trans_df);
  
  Type nll=0;
  
  // transition of mode
  Type mu = Type(0.0);
  if (mode==1) {
    for (int t=1; t<timeSteps; t++) { // since t=1
      mu = alpha(0)*p(t-1)+(1-alpha(1))*(1-p(t-1));
      nll -= dnorm(logit_p(t), log(mu)-log(Type(1.0)-mu), phi, true);
    }
  }
  
  // latent variable (z = 0 or 1)
  vector<Type> z(p.size());
  z.fill(0.0);
  if (mode == 1) {
    for (int t=0; t<timeSteps; t++) {
      z(t) = CppAD::CondExpLt(p(t),Type(0.5),Type(0.0),Type(1.0));
    }
  }

  // defining covariance matrix
  matrix<Type> sigma_var(2,2);
  for (int i=0; i<2; i++) {
    for (int j=0; j<2; j++) {
      if(i==j) {
        sigma_var(i,j) = sigma(sigma_key(i))*sigma(sigma_key(j));
      } else {
        sigma_var(i,j) = rho*sigma(sigma_key(i))*sigma(sigma_key(j));
      }
    }
  }
  // correlated random walk
  using namespace density;
  MVNORM_t<Type> neg_log_densityF(sigma_var);
  vector<Type> dist1(2);
  vector<Type> dist0(2);
  vector<Type> pred_dist1(2);
  for (int t=2; t<timeSteps; t++) { // since t=2
    dist1 = x.row(t)-x.row(t-1);
    dist0 = x.row(t-1)-x.row(t-2);
    if (mode==0) {
      pred_dist1(0) = gamma(0)*(CppAD::cos(theta(0))*dist0(0)-CppAD::sin(theta(0))*dist0(1));
      pred_dist1(1) = gamma(0)*(CppAD::sin(theta(0))*dist0(0)+CppAD::cos(theta(0))*dist0(1));
    } else {
      if (mixture==1) {
        // // expected from migration mode
        pred_dist1(0) = p(t-1)*gamma(0)*(CppAD::cos(theta(0))*dist0(0)-CppAD::sin(theta(0))*dist0(1));
        pred_dist1(1) = p(t-1)*gamma(0)*(CppAD::sin(theta(0))*dist0(0)+CppAD::cos(theta(0))*dist0(1));
        // // expected from staying mode
        pred_dist1(0) += (Type(1.0)-p(t-1))*gamma(1)*(CppAD::cos(theta(1))*dist0(0)-CppAD::sin(theta(1))*dist0(1));
        pred_dist1(1) += (Type(1.0)-p(t-1))*gamma(1)*(CppAD::sin(theta(1))*dist0(0)+CppAD::cos(theta(1))*dist0(1));
      } else {
        pred_dist1(0) = z(t-1)*gamma(0)*(CppAD::cos(theta(0))*dist0(0)-CppAD::sin(theta(0))*dist0(1));
        pred_dist1(1) = z(t-1)*gamma(0)*(CppAD::sin(theta(0))*dist0(0)+CppAD::cos(theta(0))*dist0(1));
        pred_dist1(0) = z(t-1)*gamma(0)*(CppAD::cos(theta(0))*dist0(0)-CppAD::sin(theta(0))*dist0(1));
        pred_dist1(1) = z(t-1)*gamma(0)*(CppAD::sin(theta(0))*dist0(0)+CppAD::cos(theta(0))*dist0(1));
      }
    }
    nll += neg_log_densityF(dist1-pred_dist1);
  }

  // observation
  Type V_lon = Type(0.0);
  Type V_lat = Type(0.0);
  for (int i=0; i<obs_t.size(); i++) { // since t=0
    V_lon += pow(y(i,0)-x(obs_t(i),0), Type(2.0));
    V_lat += pow(y(i,1)-x(obs_t(i),1), Type(2.0));
  }
  V_lon /= obs_t.size();
  V_lat /= obs_t.size();
  Type tau_lon = pow(V_lon*(df-Type(2.0))/df, Type(0.5));
  Type tau_lat = pow(V_lat*(df-Type(2.0))/df, Type(0.5));
  if (obs_error_type == 0) { // student-t
    for (int i=0; i<obs_t.size(); i++) { // since t=0
      nll -= dt((y(i,0)-x(obs_t(i),0))/tau_lon, df, true);
      nll -= dt((y(i,1)-x(obs_t(i),1))/tau_lat, df, true);
    }
  } else {
    for (int i=0; i<obs_t.size(); i++) { // since t=0
      nll -= dnorm(y(i,0),x(obs_t(i),0), omega, true);
      nll -= dnorm(y(i,1),x(obs_t(i),1), omega, true);
    }
  }

  ADREPORT(alpha);
  // ADREPORT(p);
  ADREPORT(phi);
  ADREPORT(gamma);
  ADREPORT(sigma);
  ADREPORT(rho);
  ADREPORT(theta);
  ADREPORT(df);
  ADREPORT(omega);
  // ADREPORT(x);
  
  return nll;
}
