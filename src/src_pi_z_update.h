//=================================
// include guard
// -------------
#ifndef __PI_Z_UPDATE_H_INCLUDED__     // if .h hasn't been included yet...
#define __PI_Z_UPDATE_H_INCLUDED__     //   #define this so the compiler knows it has been included
//=================================



// forward declared dependencies
// ----------------------------
// class <class-to-include>

// Include dependencies
// --------------------
// #include <Rcpp.h>

arma::vec dpGLM_update_pi(arma::mat Z, int K, Rcpp::List fix);
arma::colvec dpGLM_update_Z(arma::colvec y, arma::mat X, arma::colvec pi, int K, arma::mat theta, Rcpp::String family);

arma::vec hdpGLM_update_pi(arma::colvec Z, arma::colvec C, int K, Rcpp::List fix);
arma::colvec hdpGLM_update_Z(arma::colvec y, arma::mat X, arma::mat W, arma::colvec C, arma::colvec pi, int K, arma::mat theta, Rcpp::String family);

//=================================
#endif
//=================================
