// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;


void progress_bar(int t, int T)
{
  // int barWidth = 70;
  // double pos = (t+2)*(double(barWidth/double(T)));
  // std::cout <<  "[";
  // for(int i = 0; i <= barWidth; i++){
  //   if(i <= int(pos)){
  //     std::cout << "="; 
  //   }else{
  //     std::cout << " ";
  //   }
  // }
  // std::cout << "] " << int((pos/barWidth)*100.0)<<" %\r";
  // std::cout.flush();
  int barWidth = 70;
  double pos = (t+2)*(double(barWidth/double(T)));
  Rcpp::Rcout <<  "[";
  for(int i = 0; i <= barWidth; i++){
    if(i <= int(pos)){
      Rcpp::Rcout << "="; 
    }else{
      Rcpp::Rcout << " ";
    }
  }
  Rcpp::Rcout << "] " << int(fmin((pos/barWidth),1.0)*100.0)<<" %\r";
  Rcpp::Rcout.flush();

}
arma::colvec set_diff(arma::colvec& v1, arma::colvec& v2)
{

  std::vector<double> a = arma::conv_to< std::vector<double> >::from(arma::sort(v1));
  std::vector<double> b = arma::conv_to< std::vector<double> >::from(arma::sort(v2));
  std::vector<double> out;

  std::set_difference(a.begin(), a.end(), b.begin(), b.end(), std::inserter(out, out.end()));

  return arma::conv_to<arma::colvec>::from(out);
}
arma::mat rmvnormArma(int n, arma::vec mu, arma::mat sigma)
{

  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  arma::mat betaK = arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
  return betaK;
}
arma::colvec inv_scaled_chisq(int n, double df, double scale)
{
  arma::colvec z = arma::ones(n);
  arma::colvec x = arma::ones(n);
  for(int i = 0; i < n; i++){
    z[i] = R::rchisq(df);
    if(z[i]==0.0) z[i]=1e-100;
    x[i] = (df*scale)/z[i];
  }
  return(x);
}
