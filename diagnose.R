pairs_stan <- function(chain, stan_model, pars) {
  energy_vec <- as.matrix(sapply(get_sampler_params(stan_model,
                                                    inc_warmup = F),
                                 function(x) x[,"energy__"]))
  pars_vec <- extract(stan_model, pars = pars, permuted = F)
  energy <- energy_vec[, chain]
  pars <- rep(dimnames(pars_vec)$parameters[1], nrow(energy_vec))
  value <- pars_vec[, chain, 1]
  for (i in 2:dim(pars_vec)[3]) {
    energy <- c(energy, energy_vec[, chain])
    pars <- c(pars, rep(dimnames(pars_vec)$parameters[i], nrow(energy_vec)))
    value <- c(value, pars_vec[, chain, i])
  }
  df <- tibble(energy = energy, pars = pars, value = value)
  p <- ggplot(df, aes(x = energy, y = value)) +
    geom_point() +
    facet_wrap(~ pars) +
    ggtitle(paste("Chain ", chain, sep = ""))
  print(p)
  return(p)
}


diagnose_bulk_ess <- function(fit, threshold = 400) {
  fit_df <- as.data.frame(fit)
  for (name in colnames(fit_df)) {
    tmp <- ess_bulk(fit_df[[name]])
    if (tmp < threshold) {
      cat(name, ": ", tmp, "\n", sep = "")
    }
  }
}

diagnose_tail_ess <- function(fit, threshold = 400) {
  fit_df <- as.data.frame(fit)
  for (name in colnames(fit_df)) {
    tmp <- ess_tail(fit_df[[name]])
    if (tmp < threshold) {
      cat(name, ": ", tmp, "\n", sep = "")
    }
  }
}
