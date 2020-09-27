#' Estimate f and theta from codominant marker data
#'
#' @export
#' @param genos A list in the format returned by read_marker_data()
#' @param prior_pi A vector specifying lower and upper limits on pi
#' @param prior_f A vector specifying lower and upper limits on f
#' @param prior_theta A vector specifying lower and upper limits on theta
#' @param ... Optional arguments passed to `rstan::sampling()`
#' @return An object of class `stanfit` returned by `rstan::sampling()`
#'
analyze_codominant <- function(genos,
                               prior_pi = list(lo = 0.1, upper = 0.9),
                               prior_f = list(lo = 0.01, upper = 0.2),
                               prior_theta = list(lo = 0.01, upper = 0.2),
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
                    N_pops = genos$N_pops,
                    n = genos$n,
                    mu_pi = logit_prior_pi$mu,
                    sd_pi = logit_prior_pi$sd,
                    mu_f = logit_prior_f$mu,
                    sd_f = logit_prior_f$sd,
                    mu_theta = logit_prior_theta$mu,
                    sd_theta = logit_prior_theta$sd)
  stan_pars <- c("f",
                 "theta",
                 "p",
                 "pi",
                 "log_lik")
  fit <- rstan::sampling(stanmodels$analyze_codominant,
                         data = stan_data,
                         pars = stan_pars,
                         init = initialize_chains,
                         ...)
  print(fit, pars = c("f", "theta", "lp__"), digits_summary = 3)
  color_scheme_set("brightblue")
  suppressMessages(
    p <- bayesplot::mcmc_intervals(fit, pars = c("f", "theta")) +
      scale_y_discrete(labels = c("f", expression(theta))) +
      xaxis_text(size = 16, family = "sans") +
      yaxis_text(size = 16, family = "sans")
  )
  print(p)
  return(fit)
}

