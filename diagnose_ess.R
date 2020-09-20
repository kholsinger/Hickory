diagnose_bulk_ess <- function(fit, threshold = 400) {
  fit_df <- as.data.frame(fit)
  for (name in colnames(fit_df)) {
    tmp <- ess_bulk(fit_df[[name]])
    if (tmp < threshold) {
      cat(name, ": ", tmp, "\n", sep = "")
    }
  }
}

diagnose_tail_ess <- function(fit, threshold = 400) {
  fit_df <- as.data.frame(fit)
  for (name in colnames(fit_df)) {
    tmp <- ess_tail(fit_df[[name]])
    if (tmp < threshold) {
      cat(name, ": ", tmp, "\n", sep = "")
    }
  }
}
