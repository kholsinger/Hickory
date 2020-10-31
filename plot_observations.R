library(ggplot2)

plot_observations <- function(genos) {
  n <- genos$n
  N <- genos$N
  N_loci <- genos$N_loci
  N_pops <- genos$N_pops
  pops <- rownames(genos$N)
  loci <- colnames(genos$N)
  for (i in 1:N_pops) {
    for (j in 1:N_loci) {
      cat(pops[i], " ", loci[j], ": ", round(n[i, j, 1]/N[i, j], 2), ",",
          round(n[i, j, 2]/N[i, j], 2), "\n", sep = "")
    }
  }
}
