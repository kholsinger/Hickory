## set default priors
##
## Note:: These will be overwritten by defaults in analyze_*()
##
prior_pi <- list(lower = 0.1, upper = 0.9)
prior_f <- list(lower = 0.01, upper = 0.2)
prior_theta <- list(lower = 0.01, upper = 0.2)
N_loci <- 5

#' Calculate the logit of p
#'
#' @export
#' @param p A frequency in (0, 1)
#' @return log(p/(1-p))
#'
logit <- function(p) {
  return(log(p/(1-p)))
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
  lo <- logit(prior$lo)
  upper <- logit(prior$upper)
  ## set standard deviation as half of distance from mean to upper bound
  ##
  mean <- (upper + lo)/2.0
  sd <- (upper - lo)/4.0
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
#' @return None
#'
set_priors <- function(prior_pi, prior_f, prior_theta, N_loci) {
  prior_pi <<- prior_pi
  prior_f <<- prior_f
  prior_theta <<- prior_theta
  N_loci <<- N_loci
}

#' Initialize pi, f, and theta based on priors
#' Note: Depends on prior_pi, prior_f, and prior_theta being defined in
#' this namespace
#'
#' @export
#' @return None
#'
initialize_chains <- function() {
  logit_prior_pi <- logit_prior(prior_pi)
  logit_prior_f <- logit_prior(prior_f)
  logit_prior_theta <- logit_prior(prior_theta)

  logit_pi <- rnorm(N_loci, logit_prior_pi$mu, logit_prior_pi$sd)
  logit_f <- rnorm(1, logit_prior_f$mu, logit_prior_f$sd)
  logit_theta <- rnorm(1, logit_prior_theta$mu, logit_prior_theta$sd)

  list(logit_pi = logit_pi,
       logit_f = logit_f,
       logit_theta = logit_theta)
}

