// estimate allele frequencies and common inbreeding coefficient
//

data {
  int<lower=0> N_loci;        // number of loci
  int<lower=0> n[N_loci, 3];       // matrix of genotype counts
}

parameters {
  real<lower=0, upper=1> w;    // auxiliary parameter transformed to f
  vector<lower=0, upper=1>[N_loci] p; // allele frequencies
}

transformed parameters {
  vector<lower=0, upper=1>[3] x[N_loci]; // vector of genotype frequencies
  real<lower=-1, upper=1> f;   // inbreeding coefficient

  f = 2.0*w - 1.0;
  for (i in 1:N_loci) {
    x[i][1] = (p[i]^2)*(1.0 - f) + f*p[i];
    x[i][2] = 2.0*p[i]*(1.0 - p[i])*(1-f);
    x[i][3] = ((1.0 - p[i])^2)*(1.0 - f) + f*(1.0 - p[i]);
  }
}

model {
  // likelihood
  //
  for (i in 1:N_loci) {
    n[i] ~ multinomial(x[i]);
  }

  // priors
  //
  p ~ beta(1.0, 1.0);
  w ~ beta(1.0, 1.0);
}
