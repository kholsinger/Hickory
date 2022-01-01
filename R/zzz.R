## set default priors
##
## Note:: These will be overwritten by defaults in analyze_*()
##
##HickoryEnv <- new.env()
prior_pi <- list(lower = 0.1, upper = 0.9)
prior_f <- list(lower = 0.01, upper = 0.2)
prior_theta <- list(lower = 0.01, upper = 0.2)
N_loci <- 5
N_pops <- 5
