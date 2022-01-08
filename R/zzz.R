HickoryEnv <- new.env()

.onLoad <- function(libname, pkgname) {
  ## set default priors
  ##
  ## Note:: These will be overwritten by defaults in analyze_*()
  ##
  assign("prior_pi", list(lower = 0.1, upper = 0.9), envir = HickoryEnv)
  assign("prior_f", list(lower = 0.01, upper = 0.2), envir = HickoryEnv)
  assign("prior_theta", list(lower = 0.01, upper = 0.2), envir = HickoryEnv)
  assign("prior_alpha", 1, envir = HickoryEnv)
  assign("N_loci", 5, envir = HickoryEnv)
  assign("N_pops", 10, envir = HickoryEnv)
}
