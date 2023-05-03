// [[Rcpp::depends(RcppArmadillo)]]

#include <fstream>
#include <iostream>
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
List PROBE_one_cpp(const arma::vec y, const arma::mat Z, const arma::colvec Wt, 
                   const arma::colvec W_var, const arma::colvec delta, 
                   const arma::colvec beta_vec, const arma::mat Z2, double sigma2,
                   const arma::colvec update_order) {
  
  // Getting the dimensions and initializing outputs
  int n = Z.n_rows, d = Z.n_cols;
  
  arma::mat coef_mat(d,1);
  coef_mat.col(0) = beta_vec;
  arma::mat coef_mat_alpha(d+1,1);
  arma::mat se(d,1,arma::fill::ones);
  arma::mat W1(n,1,arma::fill::ones);
  arma::colvec t_Wtn = 0*Wt;
  arma::colvec t_W_varn = 0*W_var;
  
  // Updating sigma2 and initial remapping
  arma::colvec Wm_start = Wt;
  arma::colvec Wm_var_start = W_var;
  arma::colvec t_Wt2 = Wm_var_start + arma::square(Wm_start);
  arma::colvec t_Wt2_sum = arma::trans(W1)*t_Wt2;
  arma::colvec tXY = arma::trans(Wm_start)*y;
  // double alpha_start = tXY(0)/t_Wt2_sum(0);
  // arma::colvec t_Wto = alpha_start*Wm_start;
  // arma::colvec t_W_varo = alpha_start*alpha_start*W_var;
  arma::colvec t_Wto = Wm_start;
  arma::colvec t_W_varo = W_var;
  arma::colvec resid = y - t_Wto;
  double new_sigma2 = std::pow(arma::norm(resid,2),2)/(n-1);
   coef_mat_alpha.row(0) = 1;
  // coef_mat *= alpha_start;
  
  // Rcout << " t_Wto: " << t_Wto(0) << " " << t_Wto(1) << " " << t_Wto(2) << "\n";
  // 
  for(int col_ind=0;col_ind<d; ++col_ind){
    
    int col = update_order(col_ind);
    
    // Getting Z_m and p_m
    arma::colvec curr_Z = Z.col(col); 
    arma::colvec curr_Z2 = Z2.col(col); 
    arma::colvec curr_delta = delta.row(col); 
    
    // Taking hypothesis m off of the expected mean and variance.
    t_Wto -= curr_Z*curr_delta(0)*coef_mat.row(col);
    t_W_varo -= curr_Z2*curr_delta(0)*(arma::square(coef_mat.row(col))*(1-curr_delta(0)));
    
    arma::colvec t_Wt  =  t_Wto ;
    arma::colvec t_W_var = t_W_varo;
    
    // if(col_ind == 0){
    //   Rcout << " t_Wto: " << t_Wt(0) << " " << t_Wt(1) << " " << t_Wt(2) << "\n";
    //   Rcout << " t_W_var: " << t_W_var(0) << " " << t_W_var(1) << " " << t_W_var(2) << "\n";
    // 
    // }
    // Calculating the second moment (variance + mean^2) and it's sum.
    arma::colvec t_Wt2 = t_W_var + arma::square(t_Wt);
    arma::colvec t_Wt2_sum = arma::trans(W1)*t_Wt2;
    
    if(t_Wt2_sum(0) > 1e-20){
    // Making the X, the (X'X) expectation matrix.
    arma::mat bXt  = arma::join_rows(curr_Z,t_Wt);
    arma::mat bXXt2 = arma::trans(bXt)*bXt;
    double D11 = bXXt2(1,0);
    bXXt2(1,0) *= curr_delta(0);
    arma::colvec y_Wt = y - t_Wtn;
    arma::colvec tXY = arma::trans(bXt)*y_Wt;
    bXXt2(1,1) = t_Wt2_sum(0);
    
    // Checking if (X'X) is invertable and calculating it's inverse
    
    // Rcout << "Col1 double before: " << col << " Wt: " << t_Wt2_sum(0) << "\n";
    arma::mat XXt_inv = arma::inv(bXXt2);
    // Rcout << "Col1 double: " << col << " Wt: " << t_Wt2_sum(0) << "\n";
    // Updating beta
    arma::colvec t_beta = XXt_inv*tXY;
    
    // Estimating the variance of beta 
    // The derivative of beta wrt W_+
    double tH_inv = 1/det(bXXt2);
    arma::colvec th_i = 2*(t_Wt*bXXt2(0,0) - curr_Z*bXXt2(1,0));
    arma::colvec td_i = 2*t_Wt*tXY(0) - y*D11 - curr_Z*tXY(1);
    arma::colvec b1=arma::square(tH_inv*(td_i - th_i*t_beta(0))) ;
    // The derivative of beta wrt W_-
    arma::mat d_beta1 = XXt_inv*arma::trans(bXt); // 2 x n
    arma::colvec b2=arma::trans(arma::square( d_beta1.row(0) ));
    // The expected variance
    bXXt2(1,0) = D11;
    arma::mat XXt_inv2 = arma::inv(bXXt2);
    arma::mat Vbt = sigma2*arma::diagvec(XXt_inv2);
    // The variance of beta
    arma::colvec beta_var = Vbt(0,0) + arma::trans(b1)*t_W_varo  + arma::trans(b2)*t_W_varn;
    arma::colvec stderrest = arma::sqrt(beta_var);
    
    // Reverse mapping step
    t_Wto *= t_beta(1);
    t_W_varo *= std::pow(t_beta(1), 2);
    t_Wtn += curr_Z*curr_delta(0)*t_beta(0);
    t_W_varn += curr_Z2*curr_delta(0)*(std::pow(t_beta(0),2)*(1-curr_delta(0)));
    
    if(col_ind < d - 1){
      arma::colvec ids = update_order(arma::span(col_ind+1, d-1));
      arma::uvec ids2 = arma::conv_to<arma::uvec>::from(arma::colvec(ids));
      coef_mat.elem(ids2) *= t_beta(1);
    }
    
    // Exporting results.
    coef_mat_alpha.row(col+1) = t_beta(1);
    coef_mat.row(col) = t_beta(0);
    se.row(col) = stderrest(0);
    }else{
      arma::mat bXXt2 = arma::trans(curr_Z)*curr_Z;
      arma::colvec y_Wt = y - t_Wtn;
      arma::colvec tXY = arma::trans(curr_Z)*y_Wt;
      
      // Rcout << "Col1 single: " << col << "Wt" << t_Wt2_sum(0) << "\n";
      arma::mat XXt_inv = 1/bXXt2;
      // Updating beta
      arma::colvec t_beta = XXt_inv*tXY;
      
      // Estimating the variance of beta 
      // The expected variance
      arma::mat Vbt = sigma2*arma::diagvec(XXt_inv);
      // The derivative of beta wrt W_-
      arma::mat d_beta1 = XXt_inv*arma::trans(curr_Z); // 1 x n
      arma::colvec b2 = arma::trans(arma::square( d_beta1.row(0) ));
      // The variance of beta
      arma::colvec beta_var = Vbt(0,0) + arma::trans(b2)*t_W_varn;
      arma::colvec stderrest = arma::sqrt(beta_var);
      
      t_Wtn += curr_Z*curr_delta(0)*t_beta(0);
      t_W_varn += curr_Z2*curr_delta(0)*(std::pow(t_beta(0),2)*(1-curr_delta(0)));
      
      // Exporting results.
      coef_mat_alpha.row(col+1) = 0;
      coef_mat.row(col) = t_beta(0);
      se.row(col) = stderrest(0);
    }
    
    // Rcout << "Col: " << col << " t_beta: " << t_beta(0) << " " << t_beta(1) << " Std: " << stderrest(0) << "\n";
  }
  
  // Rcout << "The value of b1 : " << coef_mat.col(0) << "\n";
  // Test statistics
  
  return List::create(Named("Coefficients") = coef_mat,
                      Named("Coefficients_alpha") = coef_mat_alpha,
                      Named("StdErr") = se,
                      Named("Sigma2") = new_sigma2);
}








