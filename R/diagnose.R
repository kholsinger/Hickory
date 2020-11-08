pairs_stan <- function(chain, stan_model, pars) {
  energy_vec <- as.matrix(sapply(rstan::get_sampler_params(stan_model,
                                                           inc_warmup = F),
                                 function(x) x[,"energy__"]))
  pars_vec <- rstan::extract(stan_model, pars = pars, permuted = F)
  energy <- energy_vec[, chain]
  pars <- rep(dimnames(pars_vec)$parameters[1], nrow(energy_vec))
  value <- pars_vec[, chain, 1]
  for (i in 2:dim(pars_vec)[3]) {
    energy <- c(energy, energy_vec[, chain])
    pars <- c(pars, rep(dimnames(pars_vec)$parameters[i], nrow(energy_vec)))
    value <- c(value, pars_vec[, chain, i])
  }
  df <- data.frame(energy = energy, pars = pars, value = value)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = energy, y = value)) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~ pars) +
    ggplot2::ggtitle(paste("Chain ", chain, sep = ""))
  print(p)
  return(p)
}


diagnose_bulk_ess <- function(fit, threshold = 400) {
  fit_df <- as.data.frame(fit)
  for (name in colnames(fit_df)) {
    tmp <- rstan::ess_bulk(fit_df[[name]])
    if (!is.na(tmp) && (tmp < threshold)) {
      cat(name, ": ", tmp, "\n", sep = "")
    }
  }
}

diagnose_tail_ess <- function(fit, threshold = 400) {
  fit_df <- as.data.frame(fit)
  for (name in colnames(fit_df)) {
    tmp <- rstan::ess_tail(fit_df[[name]])
    if (!is.na(tmp) && (tmp < threshold)) {
      cat(name, ": ", tmp, "\n", sep = "")
    }
  }
}
