// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "src_utils.h"
#include "src_hmc_binomial.h"
#include "src_hmc_multinomial.h"
#include "src_gibbs_gaussian.h"
#include "src_pi_z_update.h"

using namespace Rcpp;
using namespace arma;


// Global
// ------
double hdpGLM_ACCEPTANCE_COUNT = 0.0;
double hdpGLM_ACCEPTANCE_RATE_AVERAGE = 1.0;
double hdpGLM_MCMC_TRIAL = 0.0;


// {{{ ancillary }}}

void hdpGLM_display_message(String family, int burn_in, int n_iter, int iter, int K, int J, int max_active_cluster_at_a_iter, int active_clusters_at_iter, arma::colvec Z, arma::colvec C)
{

  hdpGLM_ACCEPTANCE_RATE_AVERAGE = (hdpGLM_ACCEPTANCE_RATE_AVERAGE + hdpGLM_ACCEPTANCE_COUNT/hdpGLM_MCMC_TRIAL)/2.0;

  Rcpp::Rcout << std::endl << std::endl;
  Rcpp::Rcout << "-----------------------------------------------------" <<  std::endl;
  Rcpp::Rcout << "MCMC in progress ...." << std::endl;
  Rcpp::Rcout << std::endl;
  Rcpp::Rcout << "Family of the distribution of the outcome variable of the mixture components: " << family.get_cstring() << std::endl;
  Rcpp::Rcout << std::endl;
  Rcpp::Rcout << "Burn-in: " << burn_in << std::endl;
  Rcpp::Rcout << "Number of MCMC samples: " << n_iter << std::endl;
  Rcpp::Rcout << std::endl;
  Rcpp::Rcout << "Iteration: " << iter+1 << std::endl;
  Rcpp::Rcout << std::endl;
  Rcpp::Rcout << "Acceptance rate for beta         : " << hdpGLM_ACCEPTANCE_COUNT/hdpGLM_MCMC_TRIAL  << std::endl;
  Rcpp::Rcout << "Average acceptance rate for beta : " << hdpGLM_ACCEPTANCE_RATE_AVERAGE << std::endl;
  Rcpp::Rcout << std::endl;
  Rcpp::Rcout << "Maximum Number of cluster allowed (K)                      : " << K << std::endl;
  Rcpp::Rcout << "Maximum Number of cluster activated among all contexts     : " << max_active_cluster_at_a_iter  << std::endl;
  Rcpp::Rcout << "Maximum Number of clusters active in the current iteration : " << active_clusters_at_iter  << std::endl;
  Rcpp::Rcout << std::endl;
  Rcpp::Rcout << "(displaying only clusters with more than 5% of the data)" << std::endl;
  
  for(int j = 1; j <= J; j++){
    arma::uvec idx_j = find(C == j); 
    arma::colvec Zj  = Z(idx_j);
    arma::colvec Zjstar = unique(Zj);
    uvec Zjstar_count = hist(Zj,Zjstar);
    arma::mat tabj(Zjstar.n_rows, 2);
    tabj.col(0) = Zjstar;
    tabj.col(1) = conv_to<vec>::from(Zjstar_count);
    tabj.col(1) = 100*tabj.col(1)/sum(tabj.col(1));

    arma::uvec idx_larger_clusters = find(tabj.col(1) > 5);
    arma::uvec idx_col = arma::linspace<arma::uvec>(0, 1, 2);
    arma::mat tabj_subset = tabj(idx_larger_clusters, idx_col);

    Rcpp::Rcout << "Clusters in context with index " << j << std::endl;
    Rcpp::Rcout << tabj.t() << std::endl;
  }
}
arma::mat hdpGLM_update_countZik(arma::mat countZik, arma::mat Z)
{
  for(int i = 0; i < Z.n_rows; i++){
    countZik(i, Z(i,0)-1) += 1;
  }
  return(countZik);
}
arma::mat hdpGLM_get_pik(arma::mat countZik)
{
  arma::mat pik = arma::zeros(countZik.n_rows, countZik.n_cols);
  for(int i = 0; i < countZik.n_rows; i++){
    pik.row(i) = countZik.row(i)/sum(countZik.row(i));
  }
  return(pik);
}
arma::mat hdpGLM_get_theta_active(arma::mat theta, arma::colvec Z, arma::colvec C)
{
  // Z is a n X 1 vector sampled at the current iteration
  // C is a n X 1 vectos from the data, indicating the context of observation i
  // theta is a J*K x nCov+2 matrix
  // theta.col(0) is the cluster of the betas
  // theta.col(1) is the context of the betas
  // we only save the betas in active clusters
  int J = max(theta.col(1));			// number of contexts
  int n_active=0;				// to count the total number of active clusters across contexts
  for(int j = 1; j <= J; j++){	// I count n_active b/c it is faster that resizing the matrix theta_new below
    arma::uvec            idx_j = find(C == j);   
    arma::colvec  k_active_in_j = unique(Z(idx_j));
    n_active += k_active_in_j.size();
  }
  arma::mat theta_new(n_active,theta.n_cols);	// matrix to receive te active clusters

  int count=0;
  for(int j = 1; j <= J; j++){
    arma::uvec     idx_j = find(C == j);   
    arma::colvec  Zjstar = unique(Z(idx_j));
    for(int i = 0; i < Zjstar.n_rows; i++){
      int idx_k_active_in_j = conv_to<int>::from( find(theta.col(0) == Zjstar(i) && theta.col(1) == j) ); 
      theta_new.row(count) = theta.row(idx_k_active_in_j);
      count+=1;
    }
  }

  return(theta_new);
}

// }}}


// {{{ constants and inits }}}

arma::mat hdpGLM_get_inits_theta(int J, int K,int d, String family, List fix)
{
  arma::mat theta(J*K, 2); 

  // initializing k
  for(int j = 1; j <= J; j++){
    for(int k = 1; k <= K; k++){
      int jk = K*(j-1) + k;
      theta(jk-1,0) = k;   	// we have to subtract 1 b/c indexes in C++ starts at zero
      theta(jk-1,1) = j;
    }
  }

  // initializing beta
  arma::mat beta = rmvnormArma(K*J, fix["mu_beta"], fix["Sigma_beta"]);
  theta.resize(theta.n_rows, theta.n_cols + beta.n_cols);
  for(int d = 0; d < beta.n_cols; d++){
    theta.col(d+2) = beta.col(d);
  }


  // initialize sigma if family is gaussian
  if(family == "gaussian"){
    arma::colvec sigmak = inv_scaled_chisq(K,fix["df_sigma"],fix["s2_sigma"]);
    arma::colvec sigma(K*J);
      for(int j = 1; j <= J; j++){
	for(int k = 1; k <= K; k++){
	  int jk = K*(j-1) + k;
	  sigma(jk-1) = sigmak(k-1); 
	}
      }

    theta.resize(theta.n_rows, theta.n_cols + 1);
    theta.col(theta.n_cols-1) = sigma;
  }
  return(theta);
}

arma::mat hdpGLM_get_inits_tau(int d, int Dw, String family, List fix)
{
  arma::mat tau = rmvnormArma(d+1, fix["mu_tau"], fix["Sigma_tau"]);
  tau = tau.t();
  return(tau);
}

// }}}


// {{{ update theta and tau }}}

// Theta
// -----
arma::mat hdpGLM_update_theta(arma::colvec y, arma::mat X, arma::mat W, arma::colvec C, arma::colvec Z, int K, arma::mat tau, arma::mat theta, List fix, String family, double epsilon, int leapFrog, int hmc_iter)
{
  if(family == "gaussian"){
    theta = hdpGLM_update_theta_gaussian(y, X, W, C, Z, K, tau, theta, fix);
  }
  if(family == "binomial"){
    theta = hdpGLM_update_theta_binomial(y, X, Z, K, theta, fix, epsilon, leapFrog, hmc_iter, family);
  }
  if(family == "multinomial"){
    theta = hdpGLM_update_theta_multinomial(y, X, Z, K, theta, fix, epsilon, leapFrog, hmc_iter, family);
  }
  return(theta);
}

// tau
// ---
arma::mat hdpGLM_update_tau(arma::mat W, int K, int Dx, arma::mat theta, List fix)
{
  int Dw = W.n_cols -1;
  arma::mat Sigma_tau  = fix["Sigma_tau"];
  arma::mat Sigma_beta = fix["Sigma_beta"];
  double    sigma_beta = Sigma_beta(0,0);
  arma::mat tau = arma::zeros(Dw+1, Dx+1);
  arma::mat muAk(Dw+1, K);

  arma::mat SA = ( Sigma_tau.i() * pow(sigma_beta,2) + W.t()*W ).i();
  arma::mat SigmaA = SA * pow(sigma_beta, 2);
  arma::mat Sigma_bar_tau_d = (1.0/K) * SigmaA;

  for(int d = 0; d < Dx+1; d++){
    arma::colvec betadk;
    arma::colvec mu_tau_d(Dw+1);
    for(int k = 1; k <= K; k++){
      arma::uvec idx_k = find(theta.col(0) == k);
      arma::colvec betad = theta.col(d+2);
      betadk = betad(idx_k);
      muAk.col(k-1) = SA * W.t() * betadk;
    }
    mu_tau_d   = (1.0/K) * sum(muAk, 1);
    tau.col(d) = rmvnormArma(1, mu_tau_d, Sigma_bar_tau_d).t();
  }
  return(tau);
}

// }}}


// [[Rcpp::export]]
List hdpGLM_mcmc(arma::colvec y, arma::mat X, arma::mat W, arma::colvec C, arma::colvec weights, int K, List fix, String family, List mcmc, double epsilon, int leapFrog, int n_display, int hmc_iter)
{
  // meta
  // ----
  int active_clusters_at_iter = 1;
  int max_active_cluster_at_a_iter = 1;
  int n_display_count = 0;

  // Constants
  // ---------
  int d  = X.n_cols - 1;
  int Dw = W.n_cols - 1;
  int J = W.n_rows;
  int n = X.n_rows;
  int n_iter = mcmc["n.iter"] ;
  int burn_in = mcmc["burn.in"] ;
  int N = burn_in + n_iter;

  // initialization
  // --------------
  arma::colvec Z     = arma::ones(n);
  arma::mat theta    = hdpGLM_get_inits_theta(J, K, d, family, fix);
  arma::mat tau	     = hdpGLM_get_inits_tau(d, Dw, family, fix);
  arma::mat countZik = arma::zeros(n, K);
  arma::mat pi(K, J);

  arma::mat pik(n,K);
  countZik.col(0) = arma::ones(n);

  // MCMC
  // ----
  int n_parameters = 1+1+d+1;                  // d is number of covars, +1 for interctpt +1 for column with k +1 for column with j
  if( family == "gaussian"){n_parameters+=1;} // +1 for sigma
  arma::mat samples(0, n_parameters);
  arma::mat samples_tau(n_iter, (Dw+1)*(d+1));

  // MCMC iterations
  // ---------------
  for(int iter = 0; iter < N; iter++){
    if (iter % 80 == 0) Rcpp::checkUserInterrupt();

    // sample parameters
    // -----------------
    pi	     = hdpGLM_update_pi(Z, C, K, fix);
    Z	     = hdpGLM_update_Z(y, X, W, C, pi, K, theta, family);
    theta    = hdpGLM_update_theta(y, X, W, C, Z, K, tau, theta,  fix, family, epsilon, leapFrog, hmc_iter);
    tau      = hdpGLM_update_tau(W, K, d, theta, fix);


    // saving samples
    // --------------
    if(iter+1 > burn_in){
      arma::mat theta_new = hdpGLM_get_theta_active(theta, Z, C);
      samples.resize(samples.n_rows + theta_new.n_rows, samples.n_cols);
      for(int i = 0; i < theta_new.n_rows; i++){
	samples.row(samples.n_rows - theta_new.n_rows +i) = theta_new.row(i);
      }
      samples_tau.row(iter - burn_in)= vectorise(tau).t();
    }

    // update countZik and pik
    // -----------------------
    countZik = hdpGLM_update_countZik(countZik, Z);
    pik	     = hdpGLM_get_pik(countZik);

    // meta
    // ----
    active_clusters_at_iter = 1;
    for(int j = 0; j < J; j++){
      arma::uvec    idx_j = find(C == j); 
      arma::colvec Zjstar = unique(Z(idx_j));
      int active_clusters_j_at_iter =  Zjstar.size();
      if(active_clusters_j_at_iter > active_clusters_at_iter){active_clusters_at_iter = active_clusters_j_at_iter;};
    }
    if (active_clusters_at_iter > max_active_cluster_at_a_iter){max_active_cluster_at_a_iter = active_clusters_at_iter;};
    
    // display information
    // -------------------
    n_display_count+=1;
    if(n_display_count == n_display){
      hdpGLM_display_message(family, burn_in, n_iter, iter, K, J, max_active_cluster_at_a_iter, active_clusters_at_iter, Z, C);
      n_display_count=0;
    }

    progress_bar(iter,N);
  } // end of MCMC iterations

  Rcpp::List results = Rcpp::List::create(Rcpp::Named("samples")    = samples,
					  Rcpp::Named("tau")        = samples_tau,
					  Rcpp::Named("pik")        = pik,
					  Rcpp::Named("max_active") = max_active_cluster_at_a_iter,
					  Rcpp::Named("n.iter")     = n_iter,
					  Rcpp::Named("burn.in")    = burn_in);

  hdpGLM_ACCEPTANCE_COUNT  = 0;
  hdpGLM_ACCEPTANCE_RATE_AVERAGE = 0.0;
  hdpGLM_MCMC_TRIAL = 0;

  return(results);
}
