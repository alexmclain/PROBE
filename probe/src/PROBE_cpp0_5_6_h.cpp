// [[Rcpp::depends(RcppArmadillo)]]

#include <fstream>
#include <iostream>
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
List PROBE_cpp0_5_6_h(const arma::vec y, const arma::mat Z, const arma::colvec Wt, 
                    const arma::colvec W_var, const arma::colvec delta, 
                    const arma::colvec beta_vec, const arma::mat Z2, double sigma2, 
		    const arma::mat Sigma_y_inv) {
 
  // Matrix multiplication with Z and Sigma_y_inv outside of the loop
  // Creating elements needed for the loop
  int n = Z.n_rows, d = Z.n_cols;
  arma::mat intercept(n,1,arma::fill::ones);
  arma::mat Z_int  = arma::join_rows(intercept,Z);
  arma::colvec weights = arma::diagvec(Sigma_y_inv);
  //arma::sp_mat Sigma_y_inv_sparse = arma::sp_mat(Sigma_y_inv);
  //arma::mat ZtSZ = arma::trans(Z_int)*Sigma_y_inv_sparse*Z_int;
  arma::mat ZtSZ = (Z_int%Z_int);
  ZtSZ.each_col()%= weights;
  arma::colvec Z2S = arma::trans(sum(ZtSZ,0));
  arma::mat ZS = Z_int.each_col()%weights;
  arma::colvec int2S = arma::trans(sum(ZS,0));
  arma::mat StY = weights%y;

  //arma::colvec Z2S = arma::diagvec(ZtSZ);
  //arma::colvec int2S = ZtSZ.col(0);

  // Getting the dimensions and initializing outputs
  arma::mat coef_mat(d,3);
  arma::mat se(d,3,arma::fill::ones);
  
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

    // Making the X, the (X'X) expectation matrix.
    arma::colvec XtSX_11 = Z2S.row(0);
    arma::colvec XtSX_12 = int2S.row(col+1);
    arma::colvec XtSX_13 = arma::trans(weights)*t_Wt;
    arma::colvec XtSX_22 = Z2S.row(col+1);
    arma::colvec XtSX_23 = arma::trans(weights)*(curr_Z%t_Wt);
    arma::colvec XtSX_33 = arma::trans(weights)*t_Wt2;

    arma::mat XtSX_1 = arma::join_rows(XtSX_11, XtSX_12, XtSX_13);
    arma::mat XtSX_2 = arma::join_rows(XtSX_12, XtSX_22, XtSX_23);
    arma::mat XtSX_3 = arma::join_rows(XtSX_13, XtSX_23, XtSX_33);
    arma::mat XtSX  = arma::join_cols(XtSX_1, XtSX_2, XtSX_3);

    arma::mat bXt_1  = arma::join_rows(intercept,curr_Z);
    arma::mat bXt  = arma::join_rows(bXt_1,t_Wt);
    
    // Checking if (X'X) is invertable and calculating it's inverse
    if(XtSX.is_sympd()){
      arma::mat XXt_inv = arma::inv_sympd(XtSX);
      
      // Updating beta
      arma::colvec t_beta = XXt_inv*arma::trans(bXt)*StY;

      // Estimating the covariance of beta and their SE's
      arma::mat Vbt = arma::diagvec(arma::trans(XXt_inv)*(XtSX)*XXt_inv);
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








