## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, echo = FALSE------------------------------------------------------
library(Hickory)
library(readr)

## From https://roxygen2.r-lib.org/articles/rd-formatting.html
##
tabular <- function(df, ...) {
  stopifnot(is.data.frame(df))

  align <- function(x) if (is.numeric(x)) "r" else "l"
  col_align <- purrr::map_chr(df, align)

  cols <- lapply(df, format, ...)
  contents <- do.call("paste",
    c(cols, list(sep = " \\tab ", collapse = "\\cr\n#'   ")))

  paste("#' \\tabular{", paste(col_align, collapse = ""), "}{\n#'   ",
    paste0("\\strong{", names(df), "}", sep = "", collapse = " \\tab "), " \\cr\n#'   ",
    contents, "\n#' }\n", sep = "")
}

## ---- echo = FALSE------------------------------------------------------------
dat <- read_csv(system.file("extdata", "protea_repens.csv", package = "Hickory"),
                col_types = cols())
knitr::kable(dat[c(1:5, 25:30), 1:10])

## ---- echo = FALSE------------------------------------------------------------
dat <- read_csv(system.file("extdata", "dominant_example.csv", package = "Hickory"),
                col_types = cols())
knitr::kable(dat[c(1:5, 25:30), 1:10])

## -----------------------------------------------------------------------------
options(mc.cores = parallel::detectCores())

## set the random number seed to ensure reproducible results
##
set.seed(1234)
genos <- read_marker_data(system.file("extdata", "protea_repens.csv", package = "Hickory"))
fit_free <- analyze_codominant(genos)

## ----echo = FALSE-------------------------------------------------------------
f <- rstan::extract(fit_free)$f
theta <- rstan::extract(fit_free)$theta
mu_f <- round(mean(f), 3)
mu_theta <- round(mean(theta), 3)
f_int <- round(quantile(f, c(0.025, 0.975)), 3)
theta_int <- round(quantile(theta, c(0.025, 0.975)), 3)

## -----------------------------------------------------------------------------
options(mc.cores = parallel::detectCores())

## set the random number seed to ensure reproducible results
##
set.seed(1234)
genos <- read_marker_data(system.file("extdata", "dominant_example.csv", package = "Hickory"))
fit_free <- analyze_dominant(genos)

## ----echo = FALSE-------------------------------------------------------------
f <- rstan::extract(fit_free)$f
theta <- rstan::extract(fit_free)$theta
mu_f <- round(mean(f), 3)
mu_theta <- round(mean(theta), 3)
f_int <- round(quantile(f, c(0.025, 0.975)), 3)
theta_int <- round(quantile(theta, c(0.025, 0.975)), 3)

