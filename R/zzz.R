.onLoad <- function(libname, pkgname) {
  ## set default priors
  ##
  ## Note:: These will be overwritten by defaults in analyze_*()
  ##
  prior_pi <<- list(lower = 0.1, upper = 0.9)
  prior_f <<- list(lower = 0.01, upper = 0.2)
  prior_theta <<- list(lower = 0.01, upper = 0.2)
  prior_alpha <<- 1
  N_loci <<- 5
  N_pops <<- 10
}
