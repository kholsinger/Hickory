// estimate allele frequencies, inbreeding coefficient, and theta (Fst)
// from dominant markers
//

data {
  int<lower=0> N_loci;        // number of loci
  int<lower=0> N_pops;        // number of populations
  int<lower=0> n[N_pops, N_loci]; // counts of dominant phenotypes
  int<lower=0> N[N_pops, N_loci]; // counts of sample sizes
}

parameters {
  real<lower=0, upper=1> f;     // within-population inbreeding coefficient
  real<lower=0, upper=1> theta; // Fst - Weir & Cockerham
  vector<lower=0, upper=1>[N_loci] pi;  // mean allele frequencies
  vector<lower=0, upper=1>[N_loci] p[N_pops]; // allele frequencies by loc & pop
}

transformed parameters {
  real<lower=0, upper=1> x[N_loci, N_pops]; // dominant phenotype frequencies

  for (i in 1:N_loci) {
    for (j in 1:N_pops) { 
      x[i,j] = (p[j][i]^2)*(1.0 - f) + f*p[j][i] +
               2.0*p[j][i]*(1.0 - p[j][i])*(1-f);
    }
  }
}

model {
  // likelihood
  //
  for (i in 1:N_loci) {
    for (j in 1:N_pops) {
      n[j,i] ~ binomial(n[j,i], x[i,j]);
    }
  }

  // priors
  //
  pi ~ beta(1.0, 1.0);
  f ~ beta(1.0, 1.0);
  theta ~ beta(1.0, 1.0);
  for (i in 1:N_loci) {
    for (j in 1:N_pops) {
      p[j][i] ~ beta(((1.0 - theta)/theta)*pi[i],
                    ((1.0 - theta)/theta)*(1.0 - pi[i]));
    }
  }
}
