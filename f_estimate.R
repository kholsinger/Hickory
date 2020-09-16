library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#' Estimate inbreeding coefficient at one locus with two alleles
#'
#' @export
#' @param n Numeric vector of genotype counts
#' @param ... Arguments passed to `rstan::sampling` (e.g., iter, chains)
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
f_estimate <- function(n) {
  stan_data <- list(N_loci = nrow(n),
                    n = n)
  stan_pars <- c("p",
                 "f",
                 "w")
  fit <- stan(file = "f_estimate.stan",
              data = stan_data,
              pars = stan_pars,
              refresh = 0,
              control = list(adapt_delta = 0.999,
                             max_treedepth = 20))
  return(fit)
}

