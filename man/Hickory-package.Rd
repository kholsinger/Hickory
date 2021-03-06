% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Hickory-package.R
\docType{package}
\name{Hickory-package}
\alias{Hickory-package}
\alias{_PACKAGE}
\alias{Hickory}
\title{The 'Hickory' package.}
\description{
Hickory provides functions to read genotype data from simple
CSV files and to estimate the within population inbreeding coefficient, f,
Wright's Fst, theta.

\code{read_marker_data()} reads data from a CSV file. Individuals are in rows.
The population label is in a column named \code{pop}. All other columns are
assumed to be marker information.

For co-dominant markers genotypes are labeled 0, 1, 2, with 0 and 2 being
the alternative homozygotes and 1 being the heterozygote. I plan (hope)
to extend the model for codominant markers to allow for multiallelic
genotypes.

For dominant markers phenotypes are labeled 0 and 1, with 1 being the
dominant phenotype, i.e., the presence of the marker.
\subsection{Setting priors}{

The priors on the mean allele frequency across loci, f, and theta are
set by specifying their means and standard deviations on the logit scale.
To make it easier to think about them, they are specified by setting a
lower and an upper bound. These roughly the central 95\% interval for the
prior probability of each of the parameters.

For example, given \code{lo} and \code{upper} for one of these paramenters, the
prior mean and prior standard deviation are set as follows\preformatted{mean = (logit(upper) + logit(lo))/2.0
sd <- (logit(upper) - logit(lo))/4.0
}
}

\subsection{Checking convergence}{

It's conceivable that you'll see a warning message that says something
about Bayesian fraction missing information. If you do, let me know. I'd
like to explore those data. You may be able to get rid of the warning by
passing \code{control = list(adapt_delta = 0.9, max_treedepth = 20)} as a
parameter to \code{analyze_codominant()} or \code{analyze_dominant()}.

There's a good chance that you'll get a warning about the Bulk Effective
Sample Size and/or the Tail Effective Sample Size being too low. If you do,
take a look at the \code{n_eff} column in the output. Find the smallest one
and increast the number of iterations from the default 2000 by a factor
big enough to make the smallest \code{n_eff} comfortably bigger than 400. For
example, if the smallest \code{n_eff} is 162, you need at least a factor of 3
increase in the number of iterations. So pass \code{iter = 6000} to
\code{analyze_codominant()}or \code{analyze_dominant()}.
}
}
\examples{

genos <- read_marker_data(system.file("extdata", "protea_repens.csv", package = "Hickory"))
fit_free <- analyze_codominant(genos)
fit_f0 <- analyze_codominant(genos, f_zero = TRUE)
## compare the models using leave-one-out cross validation
##
## NOTE: The Pareto k diagnostics are something to worry about. The model
## comparisons may not be reliable. Don't trust them unless the difference
## in elpd_loo is much greater than twice its standard error
##
loo_free <- loo(fit_free)
loo_f0 <- loo(fit_f0)
## model1 in the output is the first model listed, i.e. the model where
## theta is estimated from the data
##
compare(loo_free, loo_f0)
## save the example file to disk
##
write.csv(genos, file = "protea_repens.csv", row.names = FALSE)

genos <- read_marker_data(system.file("extdata", "dominant_example.csv", package = "Hickory"))
fit_free <- analyze_dominant(genos)
fit_f0 <- analyze_dominant(genos, f_zero = TRUE)
## compare the models using leave-one-out cross validation
##
## NOTE: The Pareto k diagnostics are something to worry about. The model
## comparisons may not be reliable. Don't trust them unless the difference
## in elpd_loo is much greater than twice its standard error
##
loo_free <- loo(fit_free)
loo_f0 <- loo(fit_f0)
## model1 in the output is the first model listed, i.e. the model where
## theta is estimated from the data
##
compare(loo_free, loo_f0)
## save the example file to disk
##
write.csv(genos, file = "dominant_example.csv", row.names = FALSE)

}
\references{
Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.21.2. https://mc-stan.org
}
\author{
Kent Holsinger, \email{kent.holsinger@uconn.edu}, \url{https://orcid.org/0000-0003-4312-3804}
}
