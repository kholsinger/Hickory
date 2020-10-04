#' Leave-one-out cross validation for Hickory models
#'
#' Use leave-one-out cross validation to estimate the expected log pointwise
#' predictive density (`elpd_loo`), the effective number of parameters,
#' (`p_loo`), and the LOO information criterion (`looic`). See the
#' documentation for `loo()` in the `loo` package for details.
#'
#' @export
#' @param fit A stanfit object returned from `analyze_codominant()` or
#' `analyze_dominant()`
#' @param ... Optional parameters passed to `loo::loo()`
#' @return A named list from `loo::loo()`.
#'
loo <- function(fit, ...) {
  log_lik_array <- loo::extract_log_lik(fit, merge_chains = FALSE)
  rel_n_eff <- loo::relative_eff(x = exp(log_lik_array))
  fit_loo <- loo::loo.array(log_lik_array,
                            r_eff = rel_n_eff,
                            ...)
  return(fit_loo)
}

#' Compare two models using already calculated estimates of posterior
#' log likelihoods
#'
#' @export
#' @param model_1 An object of class "`loo`" for the first model
#' @param model_2 An object of class "`loo`" for the second model
#' @return A matrix with class "`compare.loo`"
#'
compare <- function(model_1, model_2) {
  cmp <- loo::loo_compare(model_1, model_2)
  return(cmp)
}
