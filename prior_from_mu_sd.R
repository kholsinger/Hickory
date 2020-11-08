prior_from_mu_sd <- function(mu, sd) {
  lo <- mu - 2*sd
  hi <- mu + 2*sd
  return(list(lower = max(lo, 0.001), upper = min(hi, 0.999)))
}
