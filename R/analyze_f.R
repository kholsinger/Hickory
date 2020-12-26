#' Estimate inbreeding coefficients at each population/locus combination
#' in the data set
#'
#' IMPORTANT NOTE: The underlying Stan code, while relatively straightforward
#' has not been validated. Use with caution.
#'
#' @description The summary printed before the function exits identifies the
#' loci for which the central (1-`prob`) credible interval for f does not
#' include 0.
#'
#' @export
#' @param genos A list in the format retured by `read_marker_data()`
#' @param prob Sets the credible interval to flag population locus combinations,
#' i.e., those where the central 1-`prob` credible interval does not overlap 0
#' @param ... Arguments passed to `rstan::sampling` (e.g., iter, chains)
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
analyze_f <- function(genos,
                      prob = 0.05,
                      ...)
{
  stan_data <- list(N_loci = genos$N_loci,
                    N_pops = genos$N_pops,
                    n = genos$n)
  fit <- rstan::sampling(stanmodels$analyze_f,
                         data = stan_data,
                         ...)
  return(fit)
}

#' Summarize the results of `analyze_f()`
#'
#' @export
#' @param fit A stanfit object returned by `analyze_f()`
#' @param genos The genotype counts that were analyzed
#' @param prob Sets the credible interval to flag population locus combinations,
#' i.e., those where the central 1-`prob` credible interval does not overlap 0
#' @return No values returned
#'
summarize_f <- function(fit, genos, prob = 0.05) {
  f <- rstan::extract(fit, pars = "f")$f
  for (i in 1:genos$N_pops) {
    for (j in 1:genos$N_loci) {
      interval <- quantile(f[, i, j], c(prob/2.0, 1 - prob/2.0))
      if ((interval[1] > 0.0) || (interval[2] < 0.0)) {
        cat("f[",
            rownames(genos$N)[i], ",",
            colnames(genos$N)[j], "]: ",
            round(mean(f[, i, j]), 3), " (",
            round(interval[1], 3), ",",
            round(interval[2], 3), ")\n",
            sep = "")
      }
    }
  }
}
