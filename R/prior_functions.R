#' Calculate the logit of p
#'
#' @export
#' @param p A frequency in (0, 1)
#' @return log(p/(1-p))
#'
logit <- function(p) {
  return(log(p/(1-p)))
}

#' Calculate the inverse logit of x
#'
#' @export
#' @param x log(p/(1-p))
#' @return p
#'
inv_logit <- function(x) {
  return(exp(x)/(1 + exp(x)))
}

#' Set mean and standard deviation on the logit scale given lower and upper
#' bound
#'
#' mean = (upper + lower)/2
#' sd = (upper - lower)/4
#'
#' @export
#' @param prior A list with two components: mean and upper
#' @return A vector with element 1 = mean and element 2 = sd
#'
logit_prior <- function(prior) {
  lower <- logit(prior$lower)
  upper <- logit(prior$upper)
  ## set standard deviation as half of distance from mean to upper bound
  ##
  mean <- (upper + lower)/2.0
  sd <- (upper - lower)/4.0
  return(list(mu = mean, sd = sd))
}

#' Put priors in global namespace so that they are available to
#' `initialize_chains()`.
#'
#' @export
#' @param prior_pi A vector with the mean and upper bound on pi
#' @param prior_f A vector with the mean and upper bound on f
#' @param prior_theta A vector with the mean and upper bound on theta
#' @param N_loci The number of loci in the sample
#' @param N_pops The number of populations in the sample
#' @return None
#'
set_priors <- function(prior_pi, prior_f, prior_theta, N_loci, N_pops) {
  assign("prior_pi", prior_pi, envir = HickoryEnv)
  assign("prior_f", prior_f, envir = HickoryEnv)
  assign("prior_theta", prior_theta, envir = HickoryEnv)
  assign("N_loci", N_loci, envir = HickoryEnv)
  assign("N_pops", N_pops, envir = HickoryEnv)
}

#' Put alpha ("tightness" of locus-specific theta) in global namespace
#' so that it is available to `initialize_loc_pop()`.
#'
#' @export
#' @param alpha The value of alpha to store
#' @return None
#'
set_alpha <- function(alpha) {
  assign("prior_alpha", (1.0 - alpha)/alpha, envir = HickoryEnv)
}

#' Initialize pi, f, and theta based on priors
#' Note: Depends on prior_pi, prior_f, and prior_theta being defined in
#' this namespace
#'
#' @export
#' @return list of ititial values for logit_pi, logit_f, and logit_theta
#'
initialize_chains <- function() {
  prior_pi <- get("prior_pi", envir = HickoryEnv)
  prior_f <- get("prior_f", envir = HickoryEnv)
  prior_theta <- get("prior_theta", envir = HickoryEnv)
  N_loci <- get("N_loci", envir = HickoryEnv)

  logit_prior_pi <- logit_prior(prior_pi)
  logit_prior_f <- logit_prior(prior_f)
  logit_prior_theta <- logit_prior(prior_theta)

  logit_pi <- stats::rnorm(N_loci, logit_prior_pi$mu, logit_prior_pi$sd)
  logit_f <- stats::rnorm(1, logit_prior_f$mu, logit_prior_f$sd)
  logit_theta <- stats::rnorm(1, logit_prior_theta$mu, logit_prior_theta$sd)

  return(list(logit_pi = logit_pi,
              logit_f = logit_f,
              logit_theta = logit_theta))
}

#' Initialize locus-specific thetas
#' Note: Depends on prior_theta, prior_alpha, and N_loci being defined in this
#' namespace. Implicitly depends on prior_pi and prior_f being defined in this
#' namespace.
#'
#' @export
#' @return Initial values for all internal parameters
#'
initialize_locus_pop <- function() {
  prior_theta <- get("prior_theta", envir = HickoryEnv)
  prior_alpha <- get("prior_alpha", envir = HickoryEnv)
  N_loci <- get("N_loci", envir = HickoryEnv)
  N_pops <- get("N_pops", envir = HickoryEnv)

  chains <- initialize_chains()
  theta <- (prior_theta$upper + prior_theta$lower)/2
  theta_l <- numeric(N_loci)
  pi <- numeric(N_loci)
  p <- matrix(nrow = N_pops, ncol = N_loci)
  for (i in 1:N_loci) {
    theta_l[i] <- stats::rbeta(1, theta*prior_alpha, (1-theta)*prior_alpha)
    pi[i] <- inv_logit(chains$logit_pi[i])
    for (j in 1:N_pops) {
      p[j, i] <- stats::rbeta(1, ((1-theta_l[i])/theta_l[i])*pi[i],
                       ((1-theta_l[i])/theta_l[i])*(1-pi[i]))
    }
  }
  return(list(logit_pi = chains$logit_pi,
              logit_f = chains$logit_f,
              logit_theta = chains$logit_theta,
              theta_l = theta_l,
              pi = pi,
              p = p))
}
