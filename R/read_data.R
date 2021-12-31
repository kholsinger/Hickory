#' Count the genotypes at a single locus in a single population
#'
#' @export
#' @param x Vector of genotypes at a single locus in a single population
#' - 0, 1, 2 for codominant markers
#' - 0, 1 for dominant markers
#' @return A table with the number of each genotype
#'
count_genos <- function(x) {
  genos <- numeric(3)
  genos[1] <- sum(x==0, na.rm=TRUE)
  genos[2] <- sum(x==1, na.rm=TRUE)
  genos[3] <- sum(x==2, na.rm=TRUE)
  genos
}

#' Read marker data from a CSV file
#'
#' @export
#' @param filename Full path name to a CSV file containing marker data
#' @return A list with components
#' - N_pops for the number of populations in the sample
#' - N_loci for the number of genotypes in the sample
#' - n A N_pops x N_loci x 3 array of genotype counts
#' - N A N_pops x N_loci matrix of sample sizes
#'
read_marker_data <- function(filename) {
  markers <- readr::read_csv(filename, na = ".", col_types = readr::cols())

  N_pops <-length(unique(markers$pop))
  N_loci <- ncol(markers) - 1

  locus <- colnames(markers)

  ## needed to shut up devtools::check()
  ##
  utils::globalVariables("pop")

  n <- array(dim=c(N_pops, N_loci, 3))
  N <- matrix(nrow=N_pops, ncol=N_loci)
  for (i in 1:N_pops) {
    pop_n <- unique(markers$pop)[i]
    ## x has one row for each individual in the population
    ## each locus is in a different column
    ##
    x <- subset(markers, pop==pop_n)
    for (j in 1:N_loci) {
      ## +1 because first columnt is population
      ##
      n[i, j, ] <- count_genos(x[, j+1])
      N[i,j] <- sum(n[i, j, ])
    }
  }
  rownames(N) <- unique(markers$pop)
  ## -1 because first column is pop
  ##
  colnames(N) <- colnames(markers)[-1]
  list(N_pops = N_pops,
       N_loci = N_loci,
       n = n,
       N = N)
}


