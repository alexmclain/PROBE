// [[Rcpp::depends(RcppArmadillo)]]

#include <fstream>
#include <iostream>
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
List inv_cpp(const arma::mat Sigma_y) {
    
  int n = Sigma_y.n_rows;
  arma::vec v(n, arma::fill::ones);
  
  if(Sigma_y.is_sympd()) {
    arma::mat Sigma_y_inv = arma::inv_sympd(Sigma_y);
    return List::create(Named("inverse") = Sigma_y_inv);
  } else {
    arma::mat Sigma_y_inv = arma::diagmat(v);
    return List::create(Named("inverse") = Sigma_y_inv);
  }

}







   
  
