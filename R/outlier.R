#' Report locus or population outliers
#'
#' @export
#' @param fit The model results to check for outliers
#' @param genos The original genotype data
#' @param alpha Use 1 - alpha central credible interval to identify outliers
#' @param locus If TRUE report locus outliers
#' @param pop If TRUE report population outliers
#'
#' Locus or population specific estimates of theta are reported if the credible interval for the difference
#' between those estimates and the overall estimate of theta does not overlap 0. The report lists the
#' locus or population identified as an outlier, the mean difference between the specific and overall
#' estimate and the (alpha/2, 1 - alpha/2) credible interval for the difference.
#'
report_outliers <- function(fit,
                            genos,
                            alpha = 0.05,
                            locus = TRUE,
                            pop = TRUE)
{
  if (locus) {
    labels <- colnames(genos$N)
    theta <- rstan::extract(fit, pars = "theta")$theta
    theta_i <- rstan::extract(fit, pars = "theta_i")$theta_i
    report_outliers_(labels, theta, theta_i, alpha)
  }
  if (pop) {
    labels <- rownames(genos$N)
    theta <- rstan::extract(fit, pars = "theta")$theta
    theta_j <- rstan::extract(fit, pars = "theta_j")$theta_j
    report_outliers_(labels, theta, theta_j, alpha)
  }
}

#' Report outliers. Not intended for direct use. Use report_outliers().
#'
#' @export
#' @param labels Labels for each item in the report
#' @param theta Central estimate
#' @param theta_i Individual estimate
#' @param alpha 1 - alpha central credible interval
#'
report_outliers_ <- function(labels, theta, theta_i, alpha) {
  n_labels <- length(labels)
  n_sample <- length(theta)
  for (i in 1:n_labels) {
    diff <- theta_i[, i] - theta
    interval <- quantile(diff, c(alpha/2.0, 1.0 - alpha/2.0))
    if ((interval[1] > 0) || (interval[2] < 0)) {
      cat(labels[i], ": ", round(mean(diff), 3), " (",
          round(interval[1], 3), ",", round(interval[2], 3), ")\n", sep = "")
    }
  }
}
