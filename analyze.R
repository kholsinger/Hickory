library(rstan)
library(bayesplot)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

color_scheme_set("brightblue")

logit <- function(p) {
  return(log(p/(1-p)))
}

logit_prior <- function(prior) {
  mean <- logit(prior$mean)
  upper <- logit(prior$upper)
  ## set standard deviation as half of distance from mean to upper bound
  ##
  sd <- (upper - mean)/2.0
  return(c(mean, sd))
}

## put priors in global namespace so that they are accessible from
## initialize_chains()
##
set_priors <- function(prior_pi, prior_f, prior_theta, N_loci) {
  prior_pi <<- prior_pi
  prior_f <<- prior_f
  prior_theta <<- prior_theta
  N_loci <<- N_loci
}

initialize_chains <- function() {
  logit_prior_pi <- logit_prior(prior_pi)
  logit_prior_f <- logit_prior(prior_f)
  logit_prior_theta <- logit_prior(prior_theta)

  logit_pi <- rnorm(N_loci, logit_prior_pi[1], logit_prior_pi[2])
  logit_f <- rnorm(1, logit_prior_f[1], logit_prior_f[2])
  logit_theta <- rnorm(1, logit_prior_theta[1], logit_prior_theta[2])

  list(logit_pi = logit_pi,
       logit_f = logit_f,
       logit_theta = logit_theta)
}

analyze_codominant <- function(genos,
                               prior_pi = list(mean = 0.5, upper = 0.9),
                               prior_f = list(mean = 0.1, upper = 0.2),
                               prior_theta = list(mean = 0.1, upper = 0.2),
                               ...)
{
  set_priors(prior_pi = prior_pi,
             prior_f = prior_f,
             prior_theta = prior_theta,
             N_loci = genos$N_loci)
  logit_prior_pi <- logit_prior(prior_pi)
  logit_prior_f <- logit_prior(prior_f)
  logit_prior_theta <- logit_prior(prior_theta)
  stan_data <- list(N_loci = genos$N_loci,
                    N_pops = genos$N_pop,
                    n = genos$n,
                    mu_pi = logit_prior_pi[1],
                    sd_pi = logit_prior_pi[2],
                    mu_f = logit_prior_f[1],
                    sd_f = logit_prior_f[2],
                    mu_theta = logit_prior_theta[1],
                    sd_theta = logit_prior_theta[2])
  stan_pars <- c("f",
                 "theta",
                 "p",
                 "pi")
  fit <- stan(file = "theta_estimate_codominant.stan",
              data = stan_data,
              pars = stan_pars,
              init = initialize_chains,
              ...)
  print(fit, pars = c("f", "theta", "lp__"), digits_summary = 3)
  suppressMessages(
    p <- mcmc_intervals(fit, pars = c("f", "theta")) +
      scale_y_discrete(labels = c("f", expression(theta))) +
      xaxis_text(size = 16, family = "sans") +
      yaxis_text(size = 16, family = "sans")
  )
  print(p)
  return(fit)
}

analyze_dominant <- function(genos,
                             prior_pi = list(mean = 0.5, upper = 0.9),
                             prior_f = list(mean = 0.1, upper = 0.2),
                             prior_theta = list(mean = 0.1, upper = 0.2),
                             ...)
{
  set_priors(prior_pi = prior_pi,
             prior_f = prior_f,
             prior_theta = prior_theta,
             N_loci = genos$N_loci)
  logit_prior_pi <- logit_prior(prior_pi)
  logit_prior_f <- logit_prior(prior_f)
  logit_prior_theta <- logit_prior(prior_theta)
  stan_data <- list(N_loci = genos$N_loci,
                    N_pops = genos$N_pop,
                    n = genos$n[, , 2],
                    N = genos$N,
                    mu_pi = logit_prior_pi[1],
                    sd_pi = logit_prior_pi[2],
                    mu_f = logit_prior_f[1],
                    sd_f = logit_prior_f[2],
                    mu_theta = logit_prior_theta[1],
                    sd_theta = logit_prior_theta[2])
  stan_pars <- c("f",
                 "theta",
                 "p",
                 "pi")
  fit <- stan(file = "theta_estimate_dominant.stan",
              data = stan_data,
              pars = stan_pars,
              init = initialize_chains,
              ...)
  print(fit, pars = c("f", "theta", "lp__"), digits_summary = 3)
  suppressMessages(
    p <- mcmc_intervals(fit, pars = c("f", "theta")) +
      scale_y_discrete(labels = c("f", expression(theta))) +
      xaxis_text(size = 16, family = "sans") +
      yaxis_text(size = 16, family = "sans")
  )
  print(p)
  return(fit)
}

