
lasso <- function(Y,X){
  cv.out <- cv.glmnet(X,Y,alpha = 1,nfolds = 10,lambda.min.ratio=0.001)
  # while(cv.out$lambda.min==min(cv.out$lambda)){
  #   ### Further decreasing the lambda if the smallest lambda was the best.
  #   cv.out <- cv.glmnet(X,Y,alpha = 1,nfolds = 10,lambda=exp(c(log(cv.out$lambda.1se),seq(log(min(cv.out$lambda)),log(min(cv.out$lambda)*0.001),length.out = 100))))
  # }
  cv.out
}

adap_lasso <- function(Y,X){
  ridge1_cv <- cv.glmnet(X, Y, alpha = 0, nfolds = 10)
  best_ridge_coef <- as.numeric(coef(ridge1_cv, s = ridge1_cv$lambda.min))[-1]
  cv.out <- cv.glmnet(X, Y, alpha = 1, nfolds = 10, lambda.min.ratio = 0.001,
                      penalty.factor = 1 / abs(best_ridge_coef))
  
  # while(cv.out$lambda.min==min(cv.out$lambda)){
  #   ### Further decreasing the lambda if the smallest lambda was the best.
  #   cv.out <- cv.glmnet(X, Y, alpha = 1, nfolds = 10,
  #                       penalty.factor = 1 / abs(best_ridge_coef),
  #                       lambda=exp(c(log(cv.out$lambda.1se), 
  #                                    seq(log(min(cv.out$lambda)),
  #                                        log(min(cv.out$lambda)*0.001),
  #                                        length.out = 100))) )
  # }
  list(cv.out = cv.out, 
       best_ridge_coef = best_ridge_coef)
}

scad_func <- function(Y,X,N,M){
  cv.out2 <- cv.ncvreg(matrix(as.numeric(X),N,M),Y,family = "gaussian",penalty = "SCAD",lambda.min=0.01, max.iter = 50000)
  # while(cv.out2$lambda.min==min(cv.out2$lambda)){
  #   ### Further decreasing the lambda if the smallest lambda was the best.
  #   cv.out2 <- cv.ncvreg(matrix(as.numeric(X),N,M),Y,family = "gaussian",penalty = "SCAD",lambda=exp(seq(log(min(cv.out2$lambda)),log(min(cv.out2$lambda)/100),length.out = 20)), max.iter = 50000)
  # }
  cv.out2
}

MCP_func <- function(Y,X,N,M){
  cv.out3 <- cv.ncvreg(matrix(as.numeric(X),N,M),Y,family = "gaussian",penalty = "MCP",lambda.min=0.01, max.iter = 50000)
  # while(cv.out3$lambda.min==min(cv.out3$lambda)){
  #   ### Further decreasing the lambda if the smallest lambda was the best.
  #   cv.out3 <- cv.ncvreg(matrix(as.numeric(X),N,M),Y,family = "gaussian",penalty = "MCP",lambda=exp(seq(log(min(cv.out3$lambda)),log(min(cv.out3$lambda)/100),length.out = 20)), max.iter = 50000)
  # }
  cv.out3
}


MTR <- function(E_step, alpha, signal = NULL) {
  
  p_vals <- E_step$p_vals
  lfdr_val <- E_step$lfdr
  p_hat <- E_step$pi0
  M <- length(p_vals)
  alpha_hat <- alpha/p_hat
  
  T_R <- p.adjust(p_vals, method = "BY")
  R_BH <- p.adjust(p_vals, method = "BH")
  
  threshold <- 0
  lfdr_val[is.na(lfdr_val)] <- 1
  if (min(lfdr_val) < alpha) {
    threshold <- max(sort(lfdr_val)[cumsum(sort(lfdr_val)) < alpha*(1:M)])
  }
  
  BH_res <- data.frame(BY = 1 * I(T_R < alpha), BH = 1 * I(R_BH < alpha), 
                       LFDR = 1 * I(lfdr_val <= threshold))
  
  R_BY <- p.adjust(p_vals, method = "bonferroni")
  R2_BY <- p.adjust(p_vals, method = "holm")
  
  Bonf_res <- data.frame(Holm = 1 * I(R2_BY < alpha), Bonf = 1 * I(R_BY < alpha))
  
  BH_sum <- apply(BH_res,2,sum)
  Bonf_sum <- apply(Bonf_res,2,sum)
  if (!is.null(signal)) {
    n_signal <- !(1:M %in% signal)
    BH_sum <- c(sum(BH_res$BH), sum(BH_res$BH[signal]), sum(BH_res$BH[n_signal]))
    LFDR_res <- c(sum(BH_res$LFDR), sum(BH_res$LFDR[signal]), sum(BH_res$LFDR[n_signal]))
    BY_res <- c(sum(BH_res$BY), sum(BH_res$BY[signal]), sum(BH_res$BY[n_signal]))
    BH_sum <- data.frame(BY_sum = BY_res, LFDR_sum = LFDR_res, BH_sum = BH_sum)
    
    Bonf_sum <- c(sum(Bonf_res$Bonf), sum(Bonf_res$Bonf[signal]), sum(Bonf_res$Bonf[n_signal]))
    Holm_res2 <- c(sum(Bonf_res$Holm), sum(Bonf_res$Holm[signal]), 
                   sum(Bonf_res$Holm[n_signal]))
    
    Bonf_sum <- data.frame(Holm_sum = Holm_res2, Bonf_sum = Bonf_sum)
  }
  
  
  return(list(BH_res = BH_res, Bonf_res = Bonf_res, Bonf_sum = Bonf_sum, 
              BH_sum = BH_sum))
}

