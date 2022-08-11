#' A function to plot various results from the probe algorithm
#'
#' @param full_res The output from the probe function
#' @param test_plot A logical value (true/false) indicating whether to include test data when plotting
#' @param alpha Target false discovery rate; only used if test_plot = FALSE.
#' @return Plots a figure but does not return an object
#' @examples
#' #not run
#' #plot_probe_func(full_res, test_plot = !is.null(Y_test), alpha = alpha)
plot_probe_func <- function(full_res, test_plot, alpha) {
  
  plot_dat <- full_res$plot_dat
  
  a <- quantile(plot_dat$Pred_err, probs = 0.9)
  b <- min(plot_dat$Pred_err)
  L_pred <- plot_dat$Pred_err[length(plot_dat$Pred_err)]
  a <- max(c(a, b + (L_pred-b)*5/3 ))
  
  par(mar = c(4.2, 4.5, 0.5, 4.5), mfrow = c(1, 1))
  ylab_val <- expression(paste("Signal  ", MSE[t]))
  if (test_plot) {
    ylab_val <- expression(paste("Test  ", MSPE[t]))
  }
  
  plot(plot_dat$Iter, plot_dat$Pred_err, type = "l", ylim = c(b, a), 
       xlim = range(plot_dat$Iter), xlab = "Iteration", ylab = ylab_val, 
       lwd = 2, cex.lab = 1.4, cex.axis = 1.2, las = 1)
  
  p_vec <- plot_dat$Total_Disc
  b9 <- min(p_vec)
  a9 <- quantile(p_vec, probs = 0.9)
  
  axis(4, las = 1, at = seq(b, a, length.out = 4), 
       labels = round(seq(b9, a9, length.out = 4), 0), cex.axis = 1.2, las = 2)
  mtext("Number of rejections", side = 4, line = 2.9, cex = 1.2)
  trans_crit <- (p_vec - b9)/quantile(p_vec - b9, probs = 0.9) * (a - 
                                                                    b) + b
  col_vec <- rep(1, length(p_vec))
  if (!test_plot) {
    denom <- plot_dat$Total_Disc
    denom[denom == 0] <- 1
    col_vec <- ifelse(plot_dat$FP/denom < alpha, "grey60", 1)
  }
  points(plot_dat$Iter, trans_crit, col = col_vec, pch = 19, cex = 0.5)
  
  if (test_plot) { legend("topright", legend = c("Number of rejections \n using lfdr",
                                                 "Test MSPE"), lwd = 2, lty = c(0,1),
                          pch = c(20,-1), pt.lwd = c(2,0), cex = 1.3)
  }else{ legend("topright", legend = c("","Number of rejections \n using lfdr with:", 
                                       expression(FDR>alpha), expression(FDR<= alpha),
                                       "MSE of signal"),
                lwd = c(0,0,0,0,2), lty = c(1,1,3,3,1), col = c("white","white",1,"grey60",1),
                pch = c(-1,-1,20,20,-1), pt.lwd = c(0,0,2,2,0), cex = 1.3)}
}