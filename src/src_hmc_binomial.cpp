// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "src_utils.h"
#include "src_hmc.h"


using namespace Rcpp;
using namespace arma;


// DONE
// ----
double U_bin(colvec betak, List fix)
{
  arma::mat Sigma_beta  = Rcpp::as<Rcpp::List>(fix["fix"])["Sigma_beta"];
  arma::colvec mu_beta  = Rcpp::as<Rcpp::List>(fix["fix"])["mu_beta"];
  arma::mat Xk		= fix["Xk"];
  arma::colvec yk	= fix["yk"];
  int d			= Xk.n_cols - 1;

  double log_p = as_scalar( (-(d+1)/2)*log(2*3.141593) - (1/2)*log(det(Sigma_beta)) - ( (1/2) * ((betak - mu_beta).t() * Sigma_beta.i() * (betak- mu_beta)) ) - yk.t()*log(1+exp(- Xk * betak)) - (1-yk).t()*log(1+exp( Xk * betak))  );
  return ( - log_p);
  
} 
colvec grad_U_bin(colvec betak, List fix)
{
  arma::mat Sigma_beta  = Rcpp::as<Rcpp::List>(fix["fix"])["Sigma_beta"];
  arma::colvec mu_beta  = Rcpp::as<Rcpp::List>(fix["fix"])["mu_beta"];
  arma::mat Xk		= fix["Xk"];
  arma::colvec yk	= fix["yk"];
  int d			= Xk.n_cols - 1;
  int n                 = Xk.n_rows;
  
  arma::colvec h1(d+1), h2(d+1);

  for(int j = 0; j < d+1; j++){
    h1(j) = sum(   yk   % Xk.col(j) % (1/(1+exp( Xk * betak))) );
    h2(j) = sum( (1-yk) % Xk.col(j) % (1/(1+exp(-Xk * betak))) );
  }
  arma::colvec grad_log_p = - ((betak - mu_beta).t() * Sigma_beta.i()).t() + h1 - h2;

  return( - grad_log_p);
}
mat G_bin(colvec theta)
{
  int D = theta.n_rows;
  vec diag(D);
  diag.fill(1.0);
  mat I = diagmat(diag);
  return(I);
}
colvec q_bin(colvec theta_t, List fix)
{
  arma::mat Sigma_beta  = Rcpp::as<Rcpp::List>(fix["fix"])["Sigma_beta"];
  arma::colvec mu_beta = Rcpp::as<Rcpp::List>(fix["fix"])["mu_beta"];
  int ncols = Sigma_beta.n_cols;
  arma::mat Y = arma::randn(1, ncols);
  mat sample = arma::repmat(mu_beta, 1, 1).t() + Y * arma::chol(G_bin(theta_t));
  return (sample.row(0).t());
}


// DONE
// ----
arma::mat dpGLM_update_theta_binomial (arma::colvec y,arma::mat X,arma::colvec Z, int K, arma::mat theta, List fix, double epsilon, int leapFrog, int hmc_iter, String family)
{
  // constants
  // ---------
  int d = X.n_cols - 1;

  // indexes
  arma::uvec idx_col_X = arma::linspace<arma::uvec>(0, d, d+1);
  arma::uvec idx_col_betas = arma::linspace<arma::uvec>(1, d+1, d+1);


  // update beta for Zstar
  // ---------------------
  arma::colvec Zstar = unique(Z);
  for(int i = 0; i < Zstar.size(); i++){
    double           k = Zstar[i];
    arma::uvec   idx_k = find(Z == k);

    arma::mat       Xk = X.submat(idx_k, idx_col_X);
    arma::colvec    yk = y.elem( find (Z==k) );
    int             Nk = yk.size();
    Rcpp::List fix_expanded = Rcpp::List::create(Rcpp::Named("fix") = fix,
						 Rcpp::Named("family") = family,
						 Rcpp::Named("Xk") = Xk,
						 Rcpp::Named("yk") = yk,
						 Rcpp::Named("Nk") = Nk);
    for(int i = 0; i < hmc_iter; i++){
      arma::colvec betak = theta(k-1, span(1,d+1)).t();
      arma::colvec betak_new   = hmc_update(betak, epsilon, leapFrog, fix_expanded);
      theta(k-1, span(1,(d+1)) ) = betak_new.t();

      // update acceptance rate:
      dpGLM_MCMC_TRIAL +=1.0;
      if(any(betak != betak_new)) dpGLM_ACCEPTANCE_COUNT += 1.0;
    }
  }

  // update beta for Zstar complement
  // --------------------------------
  arma::colvec Z_all = arma::linspace<arma::colvec>(1, K, K);
  arma::colvec Zstar_complementar = set_diff(Z_all, Zstar);
  for(int i = 0; i < Zstar_complementar.size(); i++){
    int k = Zstar_complementar[i];

    // update beta
    // -----------
    arma::rowvec beta_new = rmvnormArma(1, fix["mu_beta"], fix["Sigma_beta"]);
    theta(k-1, span(1,d+1)) = beta_new;
  }

  return(theta);
}
