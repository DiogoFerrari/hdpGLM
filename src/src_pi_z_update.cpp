// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>

using namespace Rcpp;
using namespace arma;


// dpGLM
// -----
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

// hdpGLM
// ------
arma::mat hdpGLM_update_pi(arma::colvec Z, arma::colvec C, int K, Rcpp::List fix)
{
  arma::colvec idxj = unique(C) ;
  int             J = idxj.size();

  arma::mat        V	= zeros<mat>(K, J);              // for the stick -breaking construction
  arma::mat       pi	= zeros<mat>(K, J);              // probability of sampling case in cluster k (row) in context j(column)
  arma::mat        N	= zeros<mat>(K, J);              // matrix received the total number of cases classified in cluster k (row) in the context j (column)
  arma::colvec sumNl	= zeros<colvec>(J);              // cumulative: total number of cases in other clusters after cluster k in context j, i.e. l=k+1, k+2,..., K and Nl(j) = N(k+1,j) + .... + N(K,j)
  arma::colvec prod	= ones<colvec>(J);
  double       alpha	= Rcpp::as<double>(fix["alpha"]);
  
  // computing Njk (total number of cases in cluster k in context j)
  for(int j = 0; j < J; j++){
    for(int k = 0; k < K; k++) {
      for(int i = 0; i < Z.n_rows; i++) {
	if (Z[i] == k+1 && C[i] ==j+1) N(k, j) += 1;
      }
    }
  }

  // for each context...
  for(int j = 0; j < J; j++){
    // ...compute pi[0] (=pi[1] in R) 
    sumNl(j) = 0; for(int l = 1; l < K; l++){sumNl(j) += N(l,j);} // for k=0, cumulative sum from l=1 to K
    V(0, j)  = R::rbeta( 1 + N(0,j), alpha + sumNl(j));
    pi(0, j) = V(0,j);
    // computing pi[1:K-2] (=pi[2:K-1] in R)
    for(int k = 1; k <= K-2; k++){
      sumNl(j)=0; for(int l = k+1; l < K; l++){sumNl(j) += N(l, j);}
      V(k, j)  = R::rbeta(1 + N(k, j), alpha + sumNl(j));
      prod(j) *= (1 - V(k-1, j));
      pi(k, j) = V(k, j) * prod(j);
    }
    V(K-1, j) = 1;            // b/c index starts on 0, so the last element is K-1
    prod(j) *= (1 - V(K-2, j));
    pi(K-1, j) = V(K-1, j) * prod(j);
  }
  return(pi);
}
arma::colvec hdpGLM_update_Z(arma::colvec y, arma::mat X, arma::mat W, arma::colvec C, arma::mat pi, int K, arma::mat theta, String family)
{
  int                J = conv_to<colvec>::from( unique(theta.col(1)) ).n_rows;
  int                n = X.n_rows;
  arma::colvec       Z = zeros(n);
  int                D = X.n_cols - 1;
  IntegerVector oneToK = seq(1,K);
  
  for(int i = 0; i < n; i++){
    arma::vec phii = zeros(K);
    for(int k = 1; k <= K; k++){
      arma::vec pyi;
      int         j = C(i);
      int    idx_jk = conv_to<int>::from( find(theta.col(0)==k && theta.col(1)==j) );
      
      // beta_jk
      arma::colvec betajk = theta(idx_jk, span(2,D+2)).t();   // -1 b/c cpp indexes starts at zero

      if ( family == "gaussian"){
	double       sigmak = theta(idx_jk, D+3);
	arma::rowvec shati  = y.row(i) - X.row(i) * betajk;
	pyi = normpdf(shati, 0.0, sigmak);
      }
      if ( family == "binomial"){
	arma::vec nuijk = X.row(i) * betajk;
	pyi = (y.row(i)==1)%(1 / (1+exp(-nuijk))) + (y.row(i)==0)%(1 - 1 / (1+exp(-nuijk)));
      }
      phii(k-1) = conv_to<double>::from(pi(k-1, j-1) * pyi);
    }
    // sampling Zi
    // -----------
    arma::vec pik = phii/sum(phii);
    Z(i) = as<double>(sample(oneToK, 1, false, wrap(pik)));
  }
  return(Z);
}
