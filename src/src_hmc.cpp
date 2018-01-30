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
  
  if (family == "binomial"){
    return ( U_bin(theta, fix) );
  }
  if (family == "multinomial"){
    return ( U_multi(theta, fix) );
  }
} 
colvec grad_U(colvec theta, List fix)
{
  String family = fix["family"];
  
  if (family == "binomial"){
    return ( grad_U_bin(theta, fix) );
  }
  if (family == "multinomial"){
    return ( grad_U_multi(theta, fix) );
  }
}
colvec q(colvec theta_t, List fix)
{
  String family = fix["family"];
  
  if (family == "binomial"){
    return ( q_bin(theta_t, fix) );
  }
  if (family == "multinomial"){
    return ( q_multi(theta_t, fix) );
  }

}
mat G(colvec theta)
{
  String family = "binomial"; // they are all te same
  
  if (family == "binomial"){
    return ( G_bin(theta) );
  }
  if (family == "multinomial"){
    return ( G_multi(theta) );
  }
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
  colvec v_current(D), theta(D), v(D);

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
    return (theta);
  }else{
    return (theta_t);
  }
}



