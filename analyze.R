library(rstan)
library(bayesplot)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

color_scheme_set("brightblue")

analyze_codominant <- function(genos, ...) {
  stan_data <- list(N_loci = genos$N_loci,
                    N_pops = genos$N_pop,
                    n = genos$n)
  stan_pars <- c("f",
                 "theta",
                 "p",
                 "pi")
  fit <- stan(file = "theta_estimate_codominant.stan",
              data = stan_data,
              pars = stan_pars,
              ...)
  print(fit, pars = c("f", "theta"), digits_summary = 3)
  p <- mcmc_intervals(fit, pars = c("f", "theta"))
  print(p)
  return(fit)
}

analyze_dominant <- function(genos, ...) {
  stan_data <- list(N_loci = genos$N_loci,
                    N_pops = genos$N_pop,
                    n = genos$n[, , 2],
                    N = genos$N)
  stan_pars <- c("f",
                 "theta",
                 "p",
                 "pi")
  fit <- stan(file = "theta_estimate_dominant.stan",
              data = stan_data,
              pars = stan_pars,
              control = list(adapt_delta = 0.9,
                             max_treedepth = 20),
              ...)
  print(fit, pars = c("f", "theta"), digits_summary = 3)
  p <- mcmc_intervals(fit, pars = c("f", "theta"))
  print(p)
  return(fit)
}

