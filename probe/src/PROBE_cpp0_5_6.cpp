// [[Rcpp::depends(RcppArmadillo)]]

#include <fstream>
#include <iostream>
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
List PROBE_cpp0_5_6(const arma::vec y, const arma::mat Z, const arma::colvec Wt, 
                    const arma::colvec W_var, const arma::colvec delta, 
                    const arma::colvec beta_vec, const arma::mat Z2, double sigma2) {
  
  // Getting the dimensions and initializing outputs
  int n = Z.n_rows, d = Z.n_cols;
  arma::mat coef_mat(d,3);
  arma::mat se(d,3,arma::fill::ones);
  arma::mat W1(n,1,arma::fill::ones);
  
  for(int col=0;col<d; ++col){
    // Getting Z_m and p_m
    arma::colvec curr_Z = Z.col(col); 
    arma::colvec curr_Z2 = Z2.col(col); 
    arma::colvec curr_delta = delta.row(col); 
    
    // Taking hypothesis m off of the expected mean and variance.
    arma::colvec t_Wt  = Wt    - curr_Z*curr_delta(0)*beta_vec.row(col);
    arma::colvec t_Wt2 = W_var - curr_Z2*curr_delta(0)*(arma::square(beta_vec.row(col))*(1-curr_delta(0)));
    
    // Calculating the second moment (variance + mean^2) and it's sum.
    t_Wt2 = t_Wt2 + arma::square(t_Wt);
    arma::colvec t_Wt2_sum = arma::trans(W1)*t_Wt2;
    
    // Making the X, the (X'X) expectation matrix.
    arma::mat bXt_1  = arma::join_rows(W1,curr_Z);
    arma::mat bXt  = arma::join_rows(bXt_1,t_Wt);
    arma::mat bXXt = arma::trans(bXt)*bXt;
    arma::mat bXXt2= bXXt;
    bXXt2(2,2) = t_Wt2_sum(0);
    
    // Checking if (X'X) is invertable and calculating it's inverse
    if(bXXt2.is_sympd()){
      arma::mat XXt_inv = arma::inv_sympd(bXXt2);
      
      // Updating beta
      arma::colvec t_beta = XXt_inv*arma::trans(bXt)*y;
  
      // Estimating the covariance of beta and their SE's
      arma::mat Vbt = sigma2*arma::diagvec(arma::trans(XXt_inv)*(bXXt2)*XXt_inv);
      arma::colvec stderrest = arma::sqrt(Vbt);
      
      // Exporting results.
      coef_mat.row(col) = arma::trans(t_beta);
      se.row(col) = arma::trans(stderrest);
    }
  }
  
  // Test statistics
  arma::mat T_vals = coef_mat/se;
  
  return List::create(Named("Coefficients") = coef_mat,Named("StdErr") = se, 
                      Named("T_statistics") = T_vals);
}








