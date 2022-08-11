#' A function providing quantities related to the E-step after applying various corrections
#'
#' @param E_step The object returned by the E_step_func function
#' @param alpha Type I error; significance level
#' @param signal An optional argument, a vector of true non-null coefficients to calculate the true and false discovery rates
#' @return A list of results adjusted via the FDR correction, the Benjamini-Hochberg-Yekutieli BY correction, and the Bonferroni correction 
#' @examples
#' #not run
#' #MTR_res <- mtr_func(E_step, alpha, signal)
mtr_func <- function(E_step, alpha, signal = NULL) {
  
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
    threshold <- max(sort(lfdr_val)[cumsum(sort(lfdr_val)) < alpha])
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
    row.names(BH_sum) <- c("Total", "Correct", "Errors")
    row.names(Bonf_sum) <- c("Total", "Correct", "Errors")
  }
  
  
  return(list(BH_res = BH_res, Bonf_res = Bonf_res, Bonf_sum = Bonf_sum, 
              BH_sum = BH_sum))
}