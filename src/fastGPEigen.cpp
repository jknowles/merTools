// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

//' @title
//' find the Cholesky decomposition of a matrix with RcppEigen
//' @description
//' Speed up Cholesky SPD matrix decomposition
//' @param A a matrix to decompose
//' @export
//[[Rcpp::export]]
Eigen::MatrixXd rcppeigen_get_chol(const Eigen::Map<Eigen::MatrixXd> & A){
  Eigen::MatrixXd chol = A.ldlt().matrixL();
  return chol;
}
