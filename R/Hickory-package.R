#' The 'Hickory' package.
#'
#' @description Hickory provides functions to read genotype data from simple
#' CSV files and to estimate the within population inbreeding coefficient, f,
#' Wright's Fst, theta.
#'
#' `read_marker_data()` reads data from a CSV file. Individuals are in rows.
#' The population label is in a column named `pop`. All other columns are
#' assumed to be marker information.
#'
#' For co-dominant markers genotypes are labeled 0, 1, 2, with 0 and 2 being
#' the alternative homozygotes and 1 being the heterozygote. I plan (hope)
#' to extend the model for codominant markers to allow for multiallelic
#' genotypes.
#'
#' For dominant markers phenotypes are labeled 0 and 1, with 1 being the
#' dominant phenotype, i.e., the presence of the marker.
#'
#' @docType package
#' @name Hickory-package
#' @aliases Hickory
#' @useDynLib Hickory, .registration = TRUE
#' @import methods
#' @import tidyverse
#' @import bayesplot
#' @importFrom rstan sampling
#'
#' @examples
#'
#' genos <- read_marker_data(system.file("extdata", "protea_repens.csv", package = "Hickory"))
#' fit_cod <- analyze_codominant(genos)
#' ## save the example file to disk
#' ##
#' write.csv(genos, file = "protea_repens.csv", row.names = FALSE)
#'
#' genos <- read_marker_data(system.file("extdata", "dominant_example.csv", package = "Hickory"))
#' fit_dom <- analyze_dominant(genos)
#' ## save the example file to disk
#' ##
#' write.csv(genos, file = "dominant_example.csv", row.names = FALSE)
#'
#' @author Kent Holsinger, \email{kent.holsinger@uconn.edu}, <https://orcid.org/0000-0003-4312-3804>
#'
#' @references
#' Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.21.2. https://mc-stan.org
#'
#'
"_PACKAGE"

