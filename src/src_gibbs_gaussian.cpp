// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include "src_utils.h"

using namespace Rcpp;
using namespace arma;

arma::mat dpGLM_update_theta_gaussian(arma::colvec y,arma::mat X,arma::colvec Z, int K, arma::mat theta, List fix)
{
  
  dpGLM_ACCEPTANCE_COUNT +=1;
  dpGLM_MCMC_TRIAL +=1;
  dpGLM_ACCEPTANCE_RATE_AVERAGE = 1;

  // constants
  // ---------
  int d = X.n_cols - 1;

  // indexes
  arma::uvec idx_col_X = arma::linspace<arma::uvec>(0, d, d+1);
  arma::uvec idx_col_betas = arma::linspace<arma::uvec>(1, d+1, d+1);
  arma::uvec idx_col_sigma = arma::linspace<arma::uvec>(d+1 +1, d+1 +1, 1);

  // update beta for Zstar
  // ---------------------
  arma::colvec Zstar = unique(Z);
  for(int i = 0; i < Zstar.size(); i++){
    double           k = Zstar[i];
    arma::uvec   idx_k = find(Z == k);

    arma::mat       Xk = X.submat(idx_k, idx_col_X);
    arma::colvec    yk = y.elem( find (Z==k) );
    int             Nk = yk.size();
    double      sigmak = theta(k-1, d+1 +1); // k-1 b/c the index starts on 0
    arma::colvec betak = theta(k-1, span(1,d+1)).t();

    // update beta
    // -----------
    arma::mat Sigma_beta = fix["Sigma_beta"];
    // arma::mat Sk = (Sigma_beta + Xk.t()*Xk).i();                    # original
    // arma::mat Sigma_betak = Sk * (sigmak*2);                        # original
    arma::mat Sk          = (Sigma_beta.i()*pow(sigmak,2) + Xk.t()*Xk).i();
    arma::mat Sigma_betak = Sk * pow(sigmak,2);                       
    arma::colvec mu_betak = Sk * Xk.t() * yk;
    arma::rowvec beta_new = rmvnormArma(1, mu_betak, Sigma_betak);
    theta(k-1, span(1,d+1)) = beta_new;
    
    // update sigma
    // ------------
    double nu = fix["df_sigma"];
    double s2 = fix["s2_sigma"];
    double s2khat = as_scalar( (1.0/Nk) * ( (yk - Xk*betak).t() * (yk - Xk*betak) ) );
    double df = nu + Nk;
    double scale = as_scalar( (nu*s2 + Nk*s2khat)/( nu + Nk) );
    if(scale>1){scale=1;}
    theta(k-1, d+1 +1) = as_scalar( inv_scaled_chisq(1, df, scale) );
    
  }

  // update beta for Zstar complement
  // --------------------------------
  arma::colvec Z_all = arma::linspace<arma::colvec>(1, K, K);
  arma::colvec Zstar_complementar = set_diff(Z_all, Zstar);
  for(int i = 0; i < Zstar_complementar.size(); i++){
    int              k = Zstar_complementar[i];

    // update beta
    // -----------
    arma::rowvec beta_new = rmvnormArma(1, fix["mu_beta"], fix["Sigma_beta"]);
    theta(k-1, span(1,d+1)) = beta_new;

    // update sigma
    // ------------
    double df    = fix["df_sigma"];
    double scale = fix["s2_sigma"];
    theta(k-1, d+1 +1) = as_scalar( inv_scaled_chisq(1, df, scale) );
  }

  return(theta);
}
