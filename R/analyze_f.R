#' Estimate inbreeding coefficients at each population/locus combination
#' in the data set
#'
#' @export
#' @param genos A list in the format retured by `read_marker_data()`
#' @param ... Arguments passed to `rstan::sampling` (e.g., iter, chains)
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
analyze_f <- function(genos,
                      ...)
{
  stan_data <- list(N_loci = genos$N_loci,
                    N_pops = genos$N_pops,
                    n = genos$n)
  fit <- rstan::sampling(stanmodels$analyze_f,
                         data = stan_data,
                         refresh = 0,
                         ...)
  return(fit)
}

