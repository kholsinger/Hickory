library(tidyverse)
library(rstan)

count_genotypes <- function(genos) {
  n_pop <- length(unique(genos$pop))
  n_loci <- ncol(genos) - 1
  counts <- array(c(n_pop, n_loci, 3))
  for (i in 1:n_pop) {
    for (j in 1:n_loci) {
      tmp
      print(as.numeric(table(genos[i,j])))
      counts[i, j, ] <-
    }
  }
  return(counts)
}

genos <- read_csv("protea_repens.csv", na = ".")
genos$pop <- gsub("^([A-Z]+)_.*$", "\\1", genos$indiv)

counts <- count_genotypes(genos[, -1])
