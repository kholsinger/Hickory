globalVariables(names = c("prior_pi",
                          "priof_f",
                          "prior_theta",
                          "N_loci",
                          "N_pops",
                          "Hickory"),
                package = "Hickory",
                add = FALSE)
Hickory <- new.env()

.onLoad <- function(libname, pkgname) {
  ## set default priors
  ##
  ## Note:: These will be overwritten by defaults in analyze_*()
  ##
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
  assign("prior_pi", list(lower = 0.1, upper = 0.9), envir = Hickory)
  assign("prior_f", list(lower = 0.01, upper = 0.2), envir = Hickory)
  assign("prior_theta", list(lower = 0.01, upper = 0.2), envir = Hickory)
  assign("N_loci", 5, envir = Hickory)
  assign("N_pops", 5, envir = Hickory)
}
