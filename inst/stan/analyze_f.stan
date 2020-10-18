// estimate allele frequencies and inbreeding coefficients at every locus
// in every population in a sample
//

data {
  int<lower=0> N_loci;        // number of loci
  int<lower=0> N_pops;        // number of populations
  int<lower=0> n[N_pops, N_loci, 3]; // genotype counts by locus and pop
}

parameters {
  real<lower=0, upper=1> p[N_pops, N_loci]; // allele frequencies
  real<lower=0, upper=1> w[N_pops, N_loci]; // used to allow negative f
}

transformed parameters {
  vector<lower=0, upper=1>[3] x[N_pops, N_loci]; // genotype frequencies
  real<lower=-1, upper=1> f[N_pops, N_loci]; // inbreeding coefficients

  for (i in 1:N_pops) {
    for (j in 1:N_loci) {
      real f_min = max({-p[i,j]/(1.0 - p[i,j]), -(1.0 - p[i,j])/p[i,j]});
      f[i,j] = w[i,j]*(1.0 - f_min) + (1 - w[i,j])*f_min;
      x[i,j][1] = (p[i,j]^2)*(1.0 - f[i,j]) + f[i,j]*p[i,j];
      x[i,j][2] = 2.0*p[i,j]*(1.0 - p[i,j])*(1-f[i,j]);
      x[i,j][3] = ((1.0 - p[i,j])^2)*(1.0 - f[i,j]) + f[i,j]*(1.0 - p[i,j]);
    }
  }
}

model {
  // likelihood
  //
  for (i in 1:N_pops) {
    for (j in 1:N_loci) {
      n[i,j] ~ multinomial(x[i,j]);
    }
  }

  // priors
  //
  for (i in 1:N_pops) {
    for (j in 1:N_loci) {
      p[i,j] ~ uniform(0.0, 1.0);
      w[i,j] ~ uniform(0.0, 1.0);
    }
  }
}
