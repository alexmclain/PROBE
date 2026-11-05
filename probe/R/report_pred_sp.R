report_pred_sp <- function(no_improve, L, signal_track, MTR_res, report_pred, count, E_step,
                           label = "Patience") {

  patience_str <- paste0(label, "=", no_improve, "/", L)

  if (is.null(signal_track)) {
    disc <- sum(MTR_res$BH_res$LFDR)
    if (!is.null(report_pred)) {
      if (report_pred <= 1) {
        report_pred <- signif(report_pred, 3)
      }
      if (report_pred > 1) {
        report_pred <- round(report_pred, 1)
      }
      message("Iteration=", count, " Number of discoveries (using lfdr)=",
              disc, " Sum(gamma)=", round(sum(E_step$gamma), 1), " MSPE(test)=",
              report_pred, " ", patience_str, "\n")
    } else {
      message("Iteration=", count, " Number of discoveries (using lfdr)=",
              disc, " Sum(gamma)=", round(sum(E_step$gamma), 1), " ", patience_str, "\n")
    }
  } else {
    if (report_pred <= 1) {
      report_pred <- signif(report_pred, 3)
    }
    if (report_pred > 1) {
      report_pred <- round(report_pred, 1)
    }
    message("Iteration=", count, " Hyp testing (using lfdr) TP=", signal_track[3],
            " FP=", signal_track[1], " TN=", signal_track[2]-signal_track[4], " FN=", signal_track[4],
            " MSE(signal)=", report_pred, " ", patience_str,
            "\n", sep = "")
  }
}
