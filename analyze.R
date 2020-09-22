library(rstan)
library(bayesplot)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

color_scheme_set("brightblue")

logit <- function(p) {
  return(log(p/(1-p)))
}

logit_prior <- function(prior) {
  lo <- logit(prior$lo)
  hi <- logit(prior$hi)
  ## set standard deviation as half of distance from mean to upper bound
  ##
  mean <- (hi + lo)/2.0
  sd <- (hi - lo)/4.0
  return(list(mean = mean, sd = sd))
}

initialize_chains <- function() {
  logit_prior_pi <- logit_prior(prior_pi)
  logit_prior_f <- logit_prior(prior_f)
  logit_prior_theta <- logit_prior(prior_theta)

  logit_pi <- rnorm(analyze_N_loci, logit_prior_pi$mean, logit_prior_pi$sd)
  logit_f <- rnorm(1, logit_prior_f$mean, logit_prior_f$sd)
  logit_theta <- rnorm(1, logit_prior_theta$mean, logit_prior_theta$sd)

  list(logit_pi = logit_pi,
       logit_f = logit_f,
       logit_theta = logit_theta)
}

analyze_codominant <- function(genos,
                               prior_pi = list(lo = 0.1, hi = 0.9),
                               prior_f = list(lo = 0.01, hi = 0.2),
                               prior_theta = list(lo = 0.01, hi = 0.2),
                               ...)
{
  logit_prior_pi <- logit_prior(prior_pi)
  logit_prior_f <- logit_prior(prior_f)
  logit_prior_theta <- logit_prior(prior_theta)
  stan_data <- list(N_loci = genos$N_loci,
                    N_pops = genos$N_pop,
                    n = genos$n,
                    mu_pi = logit_prior_pi$mean,
                    sd_pi = logit_prior_pi$sd,
                    mu_f = logit_prior_f$mean,
                    sd_f = logit_prior_f$sd,
                    mu_theta = logit_prior_theta$mean,
                    sd_theta = logit_prior_theta$sd)
  stan_pars <- c("f",
                 "theta",
                 "p",
                 "pi")
  fit <- stan(file = "theta_estimate_codominant.stan",
              data = stan_data,
              pars = stan_pars,
#              init = initialize_chains,
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
                             prior_pi = list(lo = 0.1, hi = 0.9),
                             prior_f = list(lo = 0.01, hi = 0.2),
                             prior_theta = list(lo = 0.01, hi = 0.2),
                             ...)
{
  logit_prior_pi <- logit_prior(prior_pi)
  logit_prior_f <- logit_prior(prior_f)
  logit_prior_theta <- logit_prior(prior_theta)
  stan_data <- list(N_loci = genos$N_loci,
                    N_pops = genos$N_pop,
                    n = genos$n[, , 2],
                    N = genos$N,
                    mu_pi = logit_prior_pi$mean,
                    sd_pi = logit_prior_pi$sd,
                    mu_f = logit_prior_f$mean,
                    sd_f = logit_prior_f$sd,
                    mu_theta = logit_prior_theta$mean,
                    sd_theta = logit_prior_theta$sd)
  stan_pars <- c("f",
                 "theta",
                 "p",
                 "pi")
  fit <- stan(file = "theta_estimate_dominant.stan",
              data = stan_data,
              pars = stan_pars,
#              init = initialize_chains,
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

