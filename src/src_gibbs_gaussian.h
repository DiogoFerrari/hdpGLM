//=================================
// include guard
// -------------
#ifndef __GIBBS_GAUSSIAN_H_INCLUDED__     // if .h hasn't been included yet...
#define __GIBBS_GAUSSIAN_H_INCLUDED__     //   #define this so the compiler knows it has been included
//=================================



// forward declared dependencies
// ----------------------------
// class <class-to-include>

// Include dependencies
// --------------------
// #include <Rcpp.h>


arma::mat dpGLM_update_theta_gaussian(arma::colvec y,arma::mat X,arma::colvec Z, int K, arma::mat theta, Rcpp::List fix);
arma::mat hdpGLM_update_theta_gaussian(arma::colvec y,arma::mat X, arma::mat W, arma::colvec C, arma::colvec Z, int K, arma::mat tau, arma::mat theta, Rcpp::List fix);

//=================================
#endif
//=================================
