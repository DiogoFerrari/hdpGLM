// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "src_hmc_binomial.h"
#include "src_hmc_multinomial.h"

using namespace Rcpp;
using namespace arma;


// these depend on family
// ----------------------
double U(colvec theta, List fix)
{
  String family = fix["family"];
  double u_value = 0.0;
  
  if (family == "binomial"){
    u_value = U_bin(theta, fix) ;
  }
  if (family == "multinomial"){
    u_value = U_multi(theta, fix) ;
  }
  return(u_value);
} 
colvec grad_U(colvec theta, List fix)
{
  String family = fix["family"];
  int D = theta.n_rows;
  arma::colvec grad = arma::zeros(D);
  
  if (family == "binomial"){
    grad = grad_U_bin(theta, fix) ;
  }
  if (family == "multinomial"){
    grad = grad_U_multi(theta, fix) ;
  }
  return(grad);
}
colvec q(colvec theta_t, List fix)
{
  String family = fix["family"];
  int D = theta_t.n_rows;
  arma::colvec q_values = arma::zeros(D);
  
  if (family == "binomial"){
    q_values = q_bin(theta_t, fix) ;
  }
  if (family == "multinomial"){
    q_values = q_multi(theta_t, fix) ;
  }
  return(q_values);

}
mat G(colvec theta)
{
  String family = "binomial"; // they are all te same
  int D = theta.n_rows;
  arma::mat g_values = arma::zeros(D, D);
  
  if (family == "binomial"){
    g_values = G_bin(theta) ;
  }
  if (family == "multinomial"){
    g_values = G_multi(theta) ;
  }
  return(g_values);
}

// these does not depend on family
// -------------------------------
double Kinectic(colvec v, colvec theta)
{
  return(  as_scalar( ( v.t() * G(theta) * v)/2.0 ) );
}
double H(double U , double K)
{
  return(U + K);
}
colvec hmc_update(colvec theta_t, double epsilon, int L, List fix)
{

  // Variables
  // ---------
  double u;
  double alpha;
  int D = theta_t.n_rows;
  colvec v_current(D), theta(D), v(D), theta_return(D);

  // Algorithm
  // ---------
  v_current = q(theta_t, fix);

  v = v_current;
  theta = theta_t;

  // leapfrog method together with Modified Euller Method
  v = v - (epsilon/2.0)*grad_U(theta, fix);
  for(int l = 0; l < (L-1); l++){
    theta = theta + epsilon * (G(theta).i() * v);
    v     = v     - epsilon * grad_U(theta, fix);
  }
  theta = theta +   epsilon     * (G(theta).i() * v);
  v     = v     - (epsilon/2.0) * grad_U(theta, fix);
  v     = - v;

  u     = as_scalar(arma::randu(1));
  alpha = exp( - H(U(theta, fix), Kinectic(v, theta)) + H(U(theta_t, fix), Kinectic(v_current, theta_t)) );

  // std::cout << "theta       = " << theta.t() ;
  // std::cout << "theta_t     = " << theta_t.t() ;
  // std::cout << "v           = " << v.t() ;
  // std::cout << "v_current   = " << v_current.t() ;
  // std::cout << "U(theta)    = " << U(theta, fix) << "\n";
  // std::cout << "U(theta_t)  = " << U(theta_t, fix) << "\n";
  // std::cout << "K(v)        = " << Kinectic(v, theta) << "\n";
  // std::cout << "K(v_current)= " << Kinectic(v_current, theta_t) << "\n";
  // std::cout << "-H + H      = " << - H(U(theta, fix), Kinectic(v, theta)) + H(U(theta_t, fix), Kinectic(v_current, theta_t)) << "\n" ;
  // std::cout << "exp(-H + H) = " << exp(- H(U(theta, fix), Kinectic(v, theta)) + H(U(theta_t, fix), Kinectic(v_current, theta_t))) << "\n" ;
  // std::cout << "u           = " << u << "\n";
  // std::cout << "alpha       = " << alpha << "\n";

  if (u <=  alpha) {
    theta_return = theta;
  }else{
    theta_return = theta_t;
  }
  return(theta_return);
}



