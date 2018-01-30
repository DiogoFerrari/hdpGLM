//=================================
// include guard
// -------------
#ifndef __HMC_MULTINOMIAL_H_INCLUDED__     // if .h hasn't been included yet...
#define __HMC_MULTINOMIAL_H_INCLUDED__     //   #define this so the compiler knows it has been included
//=================================

// forward declared dependencies
// ----------------------------
// class <class-to-include>

// Include dependencies
// --------------------
// #include <Rcpp.h>

double       U_multi(arma::colvec theta, Rcpp::List fix);
arma::colvec grad_U_multi(arma::colvec theta, Rcpp::List fix);
arma::mat    G_multi(arma::colvec theta);
arma::colvec q_multi(arma::colvec theta_t, Rcpp::List fix);
arma::mat dpGLM_update_theta_multinomial (arma::colvec y,arma::mat X,arma::colvec Z, int K, arma::mat theta, Rcpp::List fix, double epsilon, int leapFrog, int hmc_iter, Rcpp::String family);

//=================================
#endif
//=================================
