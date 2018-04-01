//=================================
// include guard
// -------------
#ifndef __UTILS_H_INCLUDED__     // if .h hasn't been included yet...
#define __UTILS_H_INCLUDED__     //   #define this so the compiler knows it has been included
//=================================



// forward declared dependencies
// ----------------------------
// class <class-to-include>

// Include dependencies
// --------------------
// #include <Rcpp.h>

// Global
// ------
extern double dpGLM_ACCEPTANCE_COUNT;
extern double dpGLM_ACCEPTANCE_RATE_AVERAGE;
extern double dpGLM_MCMC_TRIAL;

extern double hdpGLM_ACCEPTANCE_COUNT;
extern double hdpGLM_ACCEPTANCE_RATE_AVERAGE;
extern double hdpGLM_MCMC_TRIAL;

void progress_bar(int t, int T);
arma::colvec set_diff(arma::colvec& v1, arma::colvec& v2);
arma::mat rmvnormArma(int n, arma::vec mu, arma::mat sigma);
arma::colvec inv_scaled_chisq(int n, double df, double scale);



//=================================
#endif
//=================================
