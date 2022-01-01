.onLoad <- function(libname, pkgname) {
  ## set default priors
  ##
  ## Note:: These will be overwritten by defaults in analyze_*()
  ##
  cat("In .onLoad()")
  Hickory <- new.env()
  assign("prior_pi", list(lower = 0.1, upper = 0.9), envir = Hickory)
  assign("prior_f", list(lower = 0.01, upper = 0.2), envir = Hickory)
  assign("prior_theta", list(lower = 0.01, upper = 0.2), envir = Hickory)
  assign("N_loci", 5, envir = Hickory)
  assign("N_pops", 5, envir = Hickory)
}

.onAttach <- function(libname, pkgname) {
  ## set default priors
  ##
  ## Note:: These will be overwritten by defaults in analyze_*()
  ##
  cat("In .onAttach()")
  Hickory <- new.env()
  assign("prior_pi", list(lower = 0.1, upper = 0.9), envir = Hickory)
  assign("prior_f", list(lower = 0.01, upper = 0.2), envir = Hickory)
  assign("prior_theta", list(lower = 0.01, upper = 0.2), envir = Hickory)
  assign("N_loci", 5, envir = Hickory)
  assign("N_pops", 5, envir = Hickory)
}
