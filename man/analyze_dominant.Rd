% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analyze_dominant.R
\name{analyze_dominant}
\alias{analyze_dominant}
\title{Estimate f and theta from dominant marker data}
\usage{
analyze_dominant(
  genos,
  prior_pi = list(lower = 0.1, upper = 0.9),
  prior_f = list(lower = 0.01, upper = 0.2),
  prior_theta = list(lower = 0.01, upper = 0.2),
  f_zero = FALSE,
  f_one = FALSE,
  theta_zero = FALSE,
  theta_lp = FALSE,
  f_pop = FALSE,
  alpha_l = 0.1,
  alpha_p = 0.1,
  ...
)
}
\arguments{
\item{genos}{A list in the format returned by read_marker_data()}

\item{prior_pi}{A vector specifying lower and upper limits on pi}

\item{prior_f}{A vector specifying lower and upper limits on f}

\item{prior_theta}{A vector specifying lower and upper limits on theta}

\item{f_zero}{TRUE for f = 0 model}

\item{f_one}{TRUE for f = 1 model}

\item{theta_zero}{TRUE for theta = 0 model}

\item{theta_lp}{TRUE to estimate locus- and population-specific effects on
theta}

\item{f_pop}{TRUE to estimate population-specific f}

\item{alpha_l}{"tightness" of prior on locus-specific differences in theta}

\item{alpha_p}{"tightness" of prior on population-specific differences in
theta}

\item{...}{Optional arguments passed to \code{rstan::sampling()}}
}
\value{
An object of class \code{stanfit} returned by \code{rstan::sampling()}
}
\description{
Estimate f and theta from dominant marker data
}
