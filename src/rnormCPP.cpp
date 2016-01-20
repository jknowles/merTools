// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//' @title
//' Faster multivariate normal sampling
//' @description
//' Speed up sampling from multivariate normal using C++
//' @param n number of samples to draw
//' @param mu vector of means to use for the distribution
//' @param sigma matrix of sigma values for distribution
//'
//' @export
// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}
