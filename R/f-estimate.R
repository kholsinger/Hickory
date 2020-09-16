#' Estimate inbreeding coefficient at one locus with two alleles
#'
#' @export
#' @param n Numeric vector of genotype counts
#' @param ... Arguments passed to `rstan::sampling` (e.g., iter, chains)
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
f_stan <- function(n, ...) {
  stan_data <- list(n = n)
  fit <- rstan::sampling(stanmodels$f-estimate,
                         data = stan_data,
                         ...)
  return(fit)
}

