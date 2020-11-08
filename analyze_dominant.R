analyze_dominant <- function(genos,
                             prior_pi = list(lower = 0.1, upper = 0.9),
                             prior_f = list(lower = 0.01, upper = 0.2),
                             prior_theta = list(lower = 0.01, upper = 0.2),
                             f_zero = FALSE,
                             f_one = FALSE,
                             theta_zero = FALSE,
                             theta_i = FALSE,
                             alpha = 0.1,
                             ...)
{
  if (f_zero && f_one) {
    stop("f_zero and f_one cannot both be TRUE")
  }
  set_priors(prior_pi = prior_pi,
             prior_f = prior_f,
             prior_theta = prior_theta,
             N_loci = genos$N_loci,
             N_pops = genos$N_pops)
  if (theta_zero) {
    return(analyze_dominant_t0(genos,
                               prior_pi,
                               prior_f,
                               ...))
  }
  logit_prior_pi <- logit_prior(prior_pi)
  logit_prior_f <- logit_prior(prior_f)
  logit_prior_theta <- logit_prior(prior_theta)
  if (theta_i) {
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
                      alpha = alpha)
    fit <- rstan::sampling(stanmodels$analyze_dominant_locus_pop,
                           data = stan_data,
                           init = initialize_chains,
                           ...)
    print(fit, pars = c("f", "theta", "lp__"), digits_summary = 3)
  } else {
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
    fit <- rstan::stan(file = "analyze_dominant_locus_pop_mixture.stan",
                           data = stan_data,
                           init = initialize_chains,
                           ...)
    print(fit, pars = c("f", "theta", "lp__"), digits_summary = 3)
  }
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
