report_func <- function(count, E_step, MTR_res, Xt_conv1, signal_track, 
                        report_pred) {
  
  CC_round <- signif(Xt_conv1, 2)
  if (is.null(signal_track)) {
    disc <- sum(MTR_res$BH_res$LFDR)
    if (!is.null(report_pred)) {
      report_pred <- signif(report_pred, 2)
      cat("Iteration=", count, "Number of discoveries (using lfdr)=", 
          disc, "Sum(delta)=", round(sum(E_step$delta), 1), " MSE(test)=", 
          report_pred, "Convergence Crit=", CC_round, "\n")
    } else {
      cat("Iteration=", count, "Number of discoveries (using lfdr)=", 
          disc, "Sum(delta)=", round(sum(E_step$delta), 1), "Convergence Crit=", 
          CC_round, "\n")
    }
  } else {
    if (report_pred <= 1) {
      report_pred <- signif(report_pred, 3)
    }
    if (report_pred > 1) {
      report_pred <- round(report_pred, 1)
    }
    cat("Iteration=", count, " Hyp testing (using lfdr) TP=", signal_track[3], 
        " FP=", signal_track[1], " TN=", signal_track[2]-signal_track[4], " FN=", signal_track[4], 
        " MSE(signal)=", report_pred, " Convergence Crit=", CC_round, 
        "\n", sep = "")
  }
}