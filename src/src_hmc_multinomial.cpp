// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;




// to be completed 




double U_multi(colvec theta, List fix)
{
  return(1.0);
} 
colvec grad_U_multi(colvec theta, List fix)
{
  return ( theta ) ;
}
mat G_multi(colvec theta)
{
  int D = theta.n_rows;
  vec diag(D);
  diag.fill(1.0);
  mat I = diagmat(diag);
  return(I);
}
colvec q_multi(colvec theta_t, List fix)
{
  mat Sigma_beta = fix["Sigma_beta"];
  vec mu_beta = fix["mu_beta"];
  int ncols = Sigma_beta.n_cols;
  arma::mat Y = arma::randn(1, ncols);
  mat sample = arma::repmat(mu_beta, 1, 1).t() + Y * arma::chol(G_multi(theta_t));
  return (sample.row(0).t());
}


arma::mat dpGLM_update_theta_multinomial (arma::colvec y,arma::mat X,arma::colvec Z, int K, arma::mat theta, List fix, double epsilon, int leapFrog, int hmc_iter, String family)
{
  return(theta);
}

arma::mat hdpGLM_update_theta_multinomial (arma::colvec y,arma::mat X,arma::colvec Z, int K, arma::mat theta, List fix, double epsilon, int leapFrog, int hmc_iter, String family)
{
  return(theta);
}
