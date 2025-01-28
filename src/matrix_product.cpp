#include <Rcpp.h>
#include <RcppEigen.h>
//[[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
// A matrix multiplication in C++ which has proven faster than the R default %*% (on one device - not extensive peformance comparison)

// [[Rcpp::export]]
NumericMatrix matrix_product(Eigen::MatrixXd& tm, Eigen::MatrixXd& tm2) {
    // Convert R matrices to Eigen matrices

    // Perform matrix multiplication
    Eigen::MatrixXd prod = tm * tm2;

    // Return the result as a NumericMatrix
    return wrap(prod);
  }


