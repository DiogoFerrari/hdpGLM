// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "src_utils.h"
#include "src_hmc_binomial.h"
#include "src_hmc_multinomial.h"
#include "src_gibbs_gaussian.h"

using namespace Rcpp;
using namespace arma;


// Global
// ------
double dpGLM_ACCEPTANCE_COUNT = 0.0;
double dpGLM_ACCEPTANCE_RATE_AVERAGE = 1.0;
double dpGLM_MCMC_TRIAL = 0.0;

// {{{ ancillary }}}

void dpGLM_display_message(String family, int burn_in, int n_iter, int iter, int K, int max_active_cluster_at_a_iter, int active_clusters_at_iter, arma::colvec Z)
{
  // compute table
  vec Zstar = unique(Z);
  uvec Zstar_count = hist(Z,Zstar);
  mat A(Zstar.n_rows, 2);
  A.col(0) = Zstar;
  A.col(1) = conv_to<vec>::from(Zstar_count);
  A.col(1) = 100*A.col(1)/sum(A.col(1));
  arma::uvec idx_larger_clusters = find(A.col(1) > 5);
  arma::uvec idx_col = arma::linspace<arma::uvec>(0, 1, 2);
  arma::mat A_subset = A(idx_larger_clusters, idx_col);

  dpGLM_ACCEPTANCE_RATE_AVERAGE = (dpGLM_ACCEPTANCE_RATE_AVERAGE + dpGLM_ACCEPTANCE_COUNT/dpGLM_MCMC_TRIAL)/2.0;

  Rcpp::Rcout << std::endl << std::endl;
  Rcpp::Rcout << "-----------------------------------------------------" <<  std::endl;
  Rcpp::Rcout << "MCMC in progress ...." << std::endl;
  Rcpp::Rcout << std::endl;
  Rcpp::Rcout << "Family of the link function of the mixture components: " << family.get_cstring();
  Rcpp::Rcout << std::endl;
  Rcpp::Rcout << "Burn-in: " << burn_in << std::endl;
  Rcpp::Rcout << "Number of MCMC samples: " << n_iter << std::endl;
  Rcpp::Rcout << std::endl;
  Rcpp::Rcout << "Iteration: " << iter+1 << std::endl;
  Rcpp::Rcout << std::endl;
  Rcpp::Rcout << "Acceptance rate for beta         : " << dpGLM_ACCEPTANCE_COUNT/dpGLM_MCMC_TRIAL  << std::endl;
  Rcpp::Rcout << "Average acceptance rate for beta : " << dpGLM_ACCEPTANCE_RATE_AVERAGE << std::endl;
  Rcpp::Rcout << std::endl;
  Rcpp::Rcout << "Maximum Number of cluster allowed (K): " << K << std::endl;
  Rcpp::Rcout << "Maximum Number of cluster activated  : " << max_active_cluster_at_a_iter  << std::endl;
  Rcpp::Rcout << "Current number of active clusters    : " << active_clusters_at_iter  << std::endl;
  Rcpp::Rcout << std::endl;
  Rcpp::Rcout << "Percentage of data classified in each clusters k at current iteraction (displaying only clusters with more than 5% of the data)" << std::endl;
  Rcpp::Rcout << A_subset.t() << std::endl;
  
}
arma::mat dpGLM_update_countZik(arma::mat countZik, arma::mat Z)
{
  for(int i = 0; i < Z.n_rows; i++){
    countZik(i, Z(i,0)-1) += 1;
  }
  return(countZik);
}
arma::mat dpGLM_get_pik(arma::mat countZik)
{
  arma::mat pik = arma::zeros(countZik.n_rows, countZik.n_cols);
  for(int i = 0; i < countZik.n_rows; i++){
    pik.row(i) = countZik.row(i)/sum(countZik.row(i));
  }
  return(pik);
}
arma::mat dpGLM_get_theta_active(arma::mat theta, arma::colvec Z)
{
  arma::colvec Zstar = unique(Z);
  arma::mat theta_new(Zstar.n_rows,theta.n_cols);
  for(int i = 0; i < Zstar.n_rows; i++){
    theta_new.row(i) = theta.row(Zstar(i)-1);
  }
  return(theta_new);
}

// }}}



// {{{ constants and inits }}}

arma::mat dpGLM_get_inits(int K,int d, String family, List fix)
{
  arma::mat theta(K, 1); 

  // initializing k
  for(int k = 0; k < K; k++){
    theta(k,0) = k+1;
  }

  // initializing beta
  arma::mat beta = rmvnormArma(K, fix["mu_beta"], fix["Sigma_beta"]);
  theta.resize(theta.n_rows,theta.n_cols + beta.n_cols);
  for(int d = 0; d < beta.n_cols; d++){
    theta.col(d+1) = beta.col(d);
  }


  // initialize sigma if family is gaussian
  if(family == "gaussian"){
    arma::colvec  sigma(K);
    sigma = inv_scaled_chisq(K,fix["df_sigma"],fix["s2_sigma"]);

    theta.resize(theta.n_rows, theta.n_cols + 1);
    theta.col(theta.n_cols-1) = sigma;
  }
  return(theta);
}

// }}}



// {{{ update parameters (theta, pi, Z) }}}


// pi and Z
// --------
arma::vec dpGLM_update_pi(arma::mat Z, int K, Rcpp::List fix)
{
  arma::colvec  V = zeros<colvec>(K);
  arma::colvec pi = zeros<colvec>(K);
  arma::vec     N = zeros<colvec>(K);
  double    alpha = Rcpp::as<double>(fix["alpha"]);
  int       sumNl = 0;
  double     prod = 1.0;
  
  // computing Nk
  for(int k = 0; k < K; k++) {
    for(int i = 0; i < Z.n_rows; i++) {
      if (Z[i] == k+1) N[k] += 1;
    }
  }

  // computing pi[0] (=p[1] in R)
  sumNl = 0; for(int l = 1; l < K; l++){sumNl += N[l];}
  V[0]  = R::rbeta( 1 + N[0], alpha + sumNl);
  pi[0] = V[0];
  // computing pi[1:K-2] (=pi[2:K-1] in R)
  for(int k = 1; k <= K-2; k++){
    sumNl=0; for(int l = k+1; l < K; l++){sumNl += N[l];}
    V[k]  = R::rbeta(1 + N[k], alpha + sumNl);
    prod *= (1 - V[k-1]);
    pi[k] = V[k] * prod;
  }
  V[K-1] = 1;            // b/c index starts on 0, so the last element is K-1
  prod *= (1 - V[K-2]);
  pi[K-1] = V[K-1] * prod;

  return(pi);
}
arma::colvec dpGLM_update_Z(arma::colvec y, arma::mat X, arma::colvec pi, int K, arma::mat theta, String family)
{
  int                n = X.n_rows;
  arma::colvec       Z = zeros(n);
  arma::mat        phi = zeros(n, K);
  int                D = X.n_cols - 1;
  IntegerVector oneToK = seq(1,K);
  
  for(int k = 0; k < K; k++){
    int          idx_k = conv_to<int>::from( find(theta.col(0)==k+1) );
    arma::colvec betak = theta(idx_k, span(1,D+1)).t();   // -1 b/c cpp indexes starts at zero
    arma::vec      pyk = zeros<vec>(n);

    // updating p(y|X,betak)
    // ---------------------
    if ( family == "gaussian"){
      double  sigmak = theta(idx_k, D+2);
      arma::colvec shatk = y - X*betak;
      pyk = normpdf(shatk, 0.0, sigmak);
    }
    if ( family == "binomial"){
      arma::colvec nuk = X * betak;
      // for(int i = 0; i < n; i++){
      // 	if (y[i] == 1) {pyk[i] = 1/(1+exp(-nuk[i]));}
      // 	if (y[i] == 0) {pyk[i] = 1 - 1/(1+exp(-nuk[i]));}
      // }
      pyk  = (y==1)%(1 / (1+exp(-nuk))) + (y==0)%(1 - 1 / (1+exp(-nuk)));
    }

    // computing phik
    phi.col(k) = pi[k]*pyk;
  }

  // sampling Zi
  // -----------
  for(int i = 0; i < n; i++){
    arma::rowvec phii = phi.row(i)/sum(phi.row(i));
    Z(i) = as<double>(sample(oneToK, 1, false, wrap(phii)));
  }
  return(Z);
}
// Theta
// -----
arma::mat dpGLM_update_theta(arma::colvec y, arma::mat X, arma::colvec Z, int K, arma::mat theta, List fix, String family, double epsilon, int leapFrog, int hmc_iter)
{
  if(family == "gaussian"){
    theta = dpGLM_update_theta_gaussian(y, X, Z, K, theta, fix);
  }
  if(family == "binomial"){
    theta = dpGLM_update_theta_binomial(y, X, Z, K, theta, fix, epsilon, leapFrog, hmc_iter, family);
  }
  if(family == "multinomial"){
    theta = dpGLM_update_theta_multinomial(y, X, Z, K, theta, fix, epsilon, leapFrog, hmc_iter, family);
  }
  return(theta);
}


// }}}


// [[Rcpp::export]]
List dpGLM_mcmc(arma::colvec y, arma::mat X, arma::colvec weights, int K, List fix, String family, List mcmc, double epsilon, int leapFrog, int n_display, int hmc_iter)
{
  // meta
  // ----
  int active_clusters_at_iter = 1;
  int max_active_cluster_at_a_iter = 1;
  int n_display_count = 0;

  // Constants
  // ---------
  int d = X.n_cols - 1;
  int n = X.n_rows;
  int n_iter = mcmc["n.iter"] ;
  int burn_in = mcmc["burn.in"] ;
  int N = burn_in + n_iter;

  // initialization
  // --------------
  arma::colvec Z = arma::ones(n);
  arma::mat theta = dpGLM_get_inits(K, d, family, fix);
  arma::mat countZik = arma::zeros(n, K);
  arma::colvec pi(n);

  arma::mat pik(n,K);
  countZik.col(0) = arma::ones(n);

  // MCMC
  // ----
  int n_parameters = 1 + d+1;   // 1=Z, d+1="d betas for the covars + intercept" 
  if( family == "gaussian"){n_parameters+=1;} // for sigma
  arma::mat samples(0, n_parameters);

  // MCMC iterations
  // ---------------
  for(int iter = 0; iter < N; iter++){
    // sample parameters
    // -----------------
    theta    = dpGLM_update_theta(y, X, Z, K, theta,  fix, family, epsilon, leapFrog, hmc_iter);
    pi	     = dpGLM_update_pi(Z, K, fix);
    Z	     = dpGLM_update_Z(y, X, pi, K, theta, family);
    // saving samples
    // --------------
    if(iter+1 > burn_in){
      arma::mat theta_new = dpGLM_get_theta_active(theta, Z);
      samples.resize(samples.n_rows + theta_new.n_rows, samples.n_cols);
      for(int i = 0; i < theta_new.n_rows; i++){
	samples.row(samples.n_rows - theta_new.n_rows +i) = theta_new.row(i);
      }
    }

    // update countZik and pik
    // -----------------------
    countZik = dpGLM_update_countZik(countZik, Z);
    pik	     = dpGLM_get_pik(countZik);

    // meta
    // ----
    arma::colvec Zstar = unique(Z);
    active_clusters_at_iter =  Zstar.size();
    if (active_clusters_at_iter > max_active_cluster_at_a_iter){max_active_cluster_at_a_iter = active_clusters_at_iter;};
    
    // display information
    // -------------------
    n_display_count+=1;
    if(n_display_count == n_display){
      dpGLM_display_message(family, burn_in, n_iter, iter, K, max_active_cluster_at_a_iter, active_clusters_at_iter, Z);
      n_display_count=0;
    }

    progress_bar(iter,N);
  } // end of MCMC iterations

  Rcpp::List results = Rcpp::List::create(Rcpp::Named("samples") = samples,
					  Rcpp::Named("pik") = pik,
					  Rcpp::Named("max_active") = max_active_cluster_at_a_iter,
					  Rcpp::Named("n.iter") = n_iter,
					  Rcpp::Named("burn.in") = burn_in);

  dpGLM_ACCEPTANCE_COUNT  = 0;
  dpGLM_ACCEPTANCE_RATE_AVERAGE = 0.0;
  dpGLM_MCMC_TRIAL = 0;

  return(results);
}
