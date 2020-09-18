// estimate allele frequencies, inbreeding coefficient, and theta (Fst)
//

data {
  int<lower=0> N_loci;        // number of loci
  int<lower=0> N_pop;         // number of populations
  vector<lower=0>[3] n[N_loci, N_pop] // genotype counts by locus and pop
}

parameters {
  real<lower=0, upper=1> w;     // convenience parameter for prior on f
  real<lower=0, upper=1> theta; // Fst - Weir & Cockerham
  vector<lower=0, upper=1> pi;  // mean allele frequencies
  vector<lower=0, upper=1>[N_loci] p[N_pop]; // allele frequencies by loc & pop
}

transformed parameters {
  real<lower=-1, upper=1> f;   // inbreeding coefficient
  vector<lower=0, upper=1>[3] x[N_loci, N_pop]; // vector of geno frequencies

  f <- 2.0*w - 1.0;
  for (i in 1:N_loci) {
    for (j in 1:N_pop) { 
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
    for (j in 1:N_pop) {
      n[i,k] ~ multinomial(x[i,j]);
    }
  }

  // priors
  //
  pi ~ beta(1.0, 1.0);
  w ~ beta(1.0, 1.0);
  theta ~ beta(1.0, 1.0);
  p ~ beta(((1.0 - theta)/theta)*pi, ((1.0 - theta)/theta)*(1.0 - pi));
}
