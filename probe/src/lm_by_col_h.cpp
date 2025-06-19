// [[Rcpp::depends(RcppArmadillo)]]

#include <fstream>
#include <iostream>
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
List LM_by_col_h(const arma::vec y, const arma::mat X, double sigma2, const arma::mat Sigma_y_inv) {
  
  int n = X.n_rows, d = X.n_cols;
  
  arma::mat coef_mat(d,2);
  arma::mat se(d,2);
  arma::mat X1(n,1,arma::fill::ones);
  arma::mat StY = Sigma_y_inv*y;

  for(int col=0;col<d; ++col){
    
    arma::mat X2=arma::join_rows(X1,X.col(col));
    arma::mat X2SX2 = arma::trans(X2)*Sigma_y_inv*X2;

    // Checking if (X'X) is invertable and calculating it's inverse
    if(X2SX2.is_sympd()){
      arma::mat X2SX2_inv = arma::inv_sympd(X2SX2);

    // Updating beta
    arma::colvec coef = X2SX2_inv*arma::trans(X2)* StY;
    arma::colvec resid = y - X2*coef; 

    double sig2 = sigma2;
    arma::colvec stderrest = 
      arma::sqrt(arma::diagvec( X2SX2_inv ));

    coef_mat.row(col) = arma::trans(coef);
    se.row(col) = arma::trans(stderrest);
    }

  }
  
  return List::create(Named("Coefficients") = coef_mat,Named("StdErr") = se);
}








