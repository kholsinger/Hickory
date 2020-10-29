#' Estimate f and theta from dominant marker data
#'
#' @export
#' @param genos A list in the format returned by read_marker_data()
#' @param prior_pi A vector specifying lower and upper limits on pi
#' @param prior_f A vector specifying lower and upper limits on f
#' @param prior_theta A vector specifying lower and upper limits on theta
#' @param f_zero TRUE for f = 0 model
#' @param f_one TRUE for f = 1 model
#' @param theta_zero TRUE for theta = 0 model
#' @param ... Optional arguments passed to `rstan::sampling()`
#' @return An object of class `stanfit` returned by `rstan::sampling()`
#'
analyze_dominant <- function(genos,
                             prior_pi = list(lower = 0.1, upper = 0.9),
                             prior_f = list(lower = 0.01, upper = 0.2),
                             prior_theta = list(lower = 0.01, upper = 0.2),
                             f_zero = FALSE,
                             f_one = FALSE,
                             theta_zero = FALSE,
                             ...)
{
  if (f_zero && f_one) {
    stop("f_zero and f_one cannot both be TRUE")
  }
  set_priors(prior_pi = prior_pi,
             prior_f = prior_f,
             prior_theta = prior_theta,
             N_loci = genos$N_loci)
  if (theta_zero) {
    return(analyze_dominant_t0(genos,
                               prior_pi,
                               prior_f,
                               ...))
  }
  logit_prior_pi <- logit_prior(prior_pi)
  logit_prior_f <- logit_prior(prior_f)
  logit_prior_theta <- logit_prior(prior_theta)
  stan_data <- list(N_loci = genos$N_loci,
                    N_pops = genos$N_pop,
                    n = genos$n[, , 2],
                    N = genos$N,
                    mu_pi = logit_prior_pi$mu,
                    sd_pi = logit_prior_pi$sd,
                    mu_f = logit_prior_f$mu,
                    sd_f = logit_prior_f$sd,
                    mu_theta = logit_prior_theta$mu,
                    sd_theta = logit_prior_theta$sd,
                    f_zero = f_zero,
                    f_one = f_one)
  fit <- rstan::sampling(stanmodels$analyze_dominant,
                         data = stan_data,
                         init = initialize_chains,
                         ...)
  print(fit, pars = c("f", "theta", "lp__"), digits_summary = 3)
  color_scheme_set("brightblue")
  suppressMessages(
    p <- bayesplot::mcmc_intervals(fit, pars = c("f", "theta")) +
      ggplot2::scale_y_discrete(labels = c("f", expression(theta))) +
      xaxis_text(size = 16, family = "sans") +
      yaxis_text(size = 16, family = "sans")
  )
  print(p)
  return(fit)
}

#' Estimate f from dominant marker data assuming theta = 0, i.e.,
#' assuming that all populations have the same allele frequency at each locus
#'
#' @export
#' @param genos A list in the format returned by read_marker_data()
#' @param prior_pi A vector specifying lower and upper limits on pi
#' @param prior_f A vector specifying lower and upper limits on f
#' @param ... Optional arguments passed to `rstan::sampling()`
#' @return An object of class `stanfit` returned by `rstan::sampling()`
#'
analyze_dominant_t0 <- function(genos,
                                prior_pi,
                                prior_f,
                                ...)
{
  logit_prior_pi <- logit_prior(prior_pi)
  logit_prior_f <- logit_prior(prior_f)
  stan_data <- list(N_loci = genos$N_loci,
                    N_pops = genos$N_pop,
                    n = genos$n[, , 2],
                    N = genos$N,
                    mu_pi = logit_prior_pi$mu,
                    sd_pi = logit_prior_pi$sd,
                    mu_f = logit_prior_f$mu,
                    sd_f = logit_prior_f$sd)
  fit <- rstan::sampling(stanmodels$analyze_dominant_t0,
                         data = stan_data,
                         init = initialize_chains,
                         ...)
  print(fit, pars = c("f", "lp__"), digits_summary = 3)
  color_scheme_set("brightblue")
  suppressMessages(
    p <- bayesplot::mcmc_intervals(fit, pars = c("f")) +
      ggplot2::scale_y_discrete(labels = c("f", expression(theta))) +
      xaxis_text(size = 16, family = "sans") +
      yaxis_text(size = 16, family = "sans")
  )
  print(p)
  return(fit)
}
