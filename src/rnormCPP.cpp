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
  arma::mat R = arma::randn(n, ncols);
  try {
    R = chol(sigma);
  } catch ( ... ) {
    sigma += arma::eye(sigma.n_rows,sigma.n_rows) * 1e-6;
    R = chol(sigma);
  }
  return arma::repmat(mu, 1, n).t() + Y * R;
}
