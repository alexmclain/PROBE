m_update_func <- function(Z,Z_2,beta_tilde, delta, beta_tilde_var=0){
  Z_delta <- MVM(Z,delta*beta_tilde)$Res 
  W_ast<- c(Row_sum(as.matrix(Z_delta))$Rowsum)
  W_ast[is.nan(W_ast) | is.na(W_ast)] <- 0
  W_ast_var <- NULL
  if(!is.null(Z_2)){
    Z_delta2 <- MVM(Z_2,beta_tilde^2*delta*(1-delta) + delta*beta_tilde_var)$Res 
    W_ast_var <- c(Row_sum(as.matrix(Z_delta2))$Rowsum)
    W_ast_var[is.nan(W_ast_var) | is.na(W_ast_var)] <- 0
  }
  return(list(W_ast=W_ast,W_ast_var=W_ast_var))
}