#' Summarize locus-, population-, or both specific estimates of theta
#'
#' @export
#' @param fit The model results to summarize
#' @param genos The original genotype data
#' @param locus If TRUE report locus-specific effects
#' @param pop If TRUE report population-specific effects
#' @param alpha 1 - alpha central credible interval
#'
report <- function(fit,
                      genos,
                      locus = FALSE,
                      pop = FALSE,
                      alpha = 0.05)
{
  if (locus) {
    labels <- colnames(genos$N)
    theta <- rstan::extract(fit, pars = "theta_l")$theta_l
    report_(labels, theta, alpha)
  }
  if (pop) {
    labels <- rownames(genos$N)
    theta <- rstan::extract(fit, pars = "theta_p")$theta_p
    report_(labels, theta, alpha)
  }
}

#' Summarize estimates. Not intended for direct use. Use summarize()
#'
#' @export
#' @param labels Labels for each item in the summary
#' @param theta Individual estimates of theta for each item
#' @param alpha 1 - alpha central credible interval
#'
report_ <- function(labels, theta, alpha = 0.05) {
  n_labels <- length(labels)
  for (i in 1:n_labels) {
    xbar <- mean(theta[, i])
    quant <- quantile(theta[ ,i], c(alpha/2.0, 1.0 - alpha/2.0))
    cat(labels[i], ": ", round(xbar, 3), " (",
        round(quant[1], 3), ",",
        round(quant[2], 3), ")\n", sep = "")
  }
}
