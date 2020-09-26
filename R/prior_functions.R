## set default priors
##
## Note:: These will be overwritten by defaults in analyze_*()
##
prior_pi <- list(mean = 0.5, upper = 0.9)
prior_f <- list(mean = 0.1, upper = 0.2)
prior_theta <- list(mean = 0.1, upper = 0.2)
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
#' @export
#' @param prior A list with two components: mean and upper
#' @return A vector with element 1 = mean and element 2 = sd
#'
logit_prior <- function(prior) {
  mean <- logit(prior$mean)
  upper <- logit(prior$upper)
  ## set standard deviation as half of distance from mean to upper bound
  ##
  sd <- (upper - mean)/2.0
  return(c(mean, sd))
}

#' Put priors in global namespace so that they are available to
#' `initialize_chains()`.
#' NOTE: This needs to be changed so that they are restricted to
#' the package namespace
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

  logit_pi <- rnorm(N_loci, logit_prior_pi[1], logit_prior_pi[2])
  logit_f <- rnorm(1, logit_prior_f[1], logit_prior_f[2])
  logit_theta <- rnorm(1, logit_prior_theta[1], logit_prior_theta[2])

  list(logit_pi = logit_pi,
       logit_f = logit_f,
       logit_theta = logit_theta)
}

