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

arma::mat hdpGLM_update_theta_gaussian(arma::colvec y,arma::mat X, arma::mat W, arma::colvec C, arma::colvec Z, int K, arma::mat tau, arma::mat theta, List fix)
{
  hdpGLM_ACCEPTANCE_COUNT +=1;
  hdpGLM_MCMC_TRIAL +=1;
  hdpGLM_ACCEPTANCE_RATE_AVERAGE = 1;

  // constants
  // ---------
  int J  = W.n_rows;
  int d  = X.n_cols - 1;
  int Dw = W.n_cols - 1;

  // indexes
  arma::uvec idx_col_X = arma::linspace<arma::uvec>(0, d, d+1);
  arma::uvec idx_col_betas = arma::linspace<arma::uvec>(1, d+1, d+1);
  arma::uvec idx_col_sigma = arma::linspace<arma::uvec>(d+1 +1, d+1 +1, 1);

  // update beta for Zjstar
  // ----------------------
  for(int j = 1; j <= J; j++){
    arma::uvec    idx_j = find(C == j);
    arma::colvec Zjstar = unique(Z.elem(idx_j));

    for(int i = 0; i < Zjstar.size(); i++){
      int               k = Zjstar[i];
      int              jk = K*(j-1) + k;
      arma::uvec    idx_k = find(Z == k);
      arma::uvec   idx_jk = find(Z == k && C ==j );
      
      arma::colvec     Wj = W.row(j-1).t();
      arma::mat        Xk = X.submat(idx_k, idx_col_X);
      arma::colvec     yk = y.elem( find (Z==k) );
      int             Njk = conv_to<colvec>::from(find(Z == k)).size();
      arma::mat       Xjk = X.submat(idx_jk, idx_col_X);
      arma::colvec    yjk = y.elem( find (Z==k && C ==j) );
      double      sigmajk = theta(jk-1, d+2 +1);           // jk-1 b/c the index starts on 0, d+2 because d(# of covars) +2 (column for C and Z in theta) +1 (intercept)
      arma::colvec betajk = theta(jk-1, span(2,d+2)).t();

      // update beta
      // -----------
      arma::mat Sigma_beta  = fix["Sigma_beta"];
      arma::mat Sjk         = ( Sigma_beta.i()*pow(sigmajk,2) + Xjk.t()*Xjk ).i();
      arma::mat Sigma_betajk = Sjk * pow(sigmajk,2);                       
      arma::colvec mu_betajk = Sjk * (Sigma_beta.i() * (Wj.t() * tau).t() + (Xjk.t() * yjk)/pow(sigmajk,2) ) * pow(sigmajk,2);
      arma::rowvec beta_new = rmvnormArma(1, mu_betajk, Sigma_betajk);
      theta(jk-1, span(2,d+2)) = beta_new;
    
      // update sigma
      // ------------
      double nu = fix["df_sigma"];
      double s2 = fix["s2_sigma"];
      double s2khat = as_scalar( (1.0/Njk) * ( (yjk - Xjk*betajk).t() * (yjk - Xjk*betajk) ) );
      double df = nu + Njk;
      double scale = as_scalar( (nu*s2 + Njk*s2khat)/( nu + Njk) );
      if(scale>1){scale=1;}
      theta(jk-1, d+3) = as_scalar( inv_scaled_chisq(1, df, scale) ); // d is the number of covars, +3 for column in thetat with k, j, and intercept of X

    }
  }
    
  // update beta for Zjstar complement
  // ---------------------------------
  arma::colvec Z_all = arma::linspace<arma::colvec>(1, K, K);
  for(int j = 1; j <= J; j++){
    arma::uvec                 idx_j = find(C == j);
    arma::colvec              Zjstar = unique(Z.elem(idx_j));
    arma::colvec Zjstar_complementar = set_diff(Z_all, Zjstar);

    for(int i = 0; i < Zjstar_complementar.size(); i++){
      int              k = Zjstar_complementar[i];
      int             jk = K*(j-1) + k;

      // update beta
      // -----------
      arma::rowvec beta_new = rmvnormArma(1, fix["mu_beta"], fix["Sigma_beta"]);
      theta(jk-1, span(2,d+2)) = beta_new;

      // update sigma
      // ------------
      double df    = fix["df_sigma"];
      double scale = fix["s2_sigma"];
      theta(jk-1, d+3) = as_scalar( inv_scaled_chisq(1, df, scale) );
    }
  }
  return(theta);
}
