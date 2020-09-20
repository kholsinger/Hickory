// estimate allele frequencies, inbreeding coefficient, and theta (Fst)
// from co-dominant markers (currently only two alleles per locus)
//

data {
  int<lower=0> N_loci;        // number of loci
  int<lower=0> N_pops;        // number of populations
  int<lower=0> n[N_pops, N_loci, 3]; // genotype counts by locus and pop
}

parameters {
  real<lower=0, upper=1> f;     // within-population inbreeding coefficient
  real<lower=0, upper=1> theta; // Fst - Weir & Cockerham
  vector<lower=0, upper=1>[N_loci] pi;  // mean allele frequencies
  vector<lower=0, upper=1>[N_loci] p[N_pops]; // allele frequencies by loc & pop
}

transformed parameters {
  vector<lower=0, upper=1>[3] x[N_loci, N_pops]; // vector of geno frequencies

  for (i in 1:N_loci) {
    for (j in 1:N_pops) { 
      x[i,j][1] = (p[j][i]^2)*(1.0 - f) + f*p[j][i];
      x[i,j][2] = 2.0*p[j][i]*(1.0 - p[j][i])*(1-f);
      x[i,j][3] = ((1.0 - p[j][i])^2)*(1.0 - f) + f*(1.0 - p[j][i]);
    }
  }
}

model {
  // likelihood
  //
  for (i in 1:N_loci) {
    for (j in 1:N_pops) {
      n[j,i] ~ multinomial(x[i,j]);
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
