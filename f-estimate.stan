// estimate allele frequency and inbreeding coefficient at one locus
//

data {
  int<lower = 0> n[3];          // vector of genotype counts
}

parameters {
  real<lower=0, upper=1> p;    // allele frequency
  real<lower=-1, upper=1> f;   // inbreeding coefficient
}

transformed parameters {
  vector<lower=0, upper=1>[3] x; // vector of genotype frequencies

  x[1] = (p^2)*(1.0 - f) + f*p;
  x[2] = 2.0*p*(1.0 - p)*(1-f);
  x[3] = ((1.0 - p)^2)*(1.0 - f) + f*(1.0 - p);
}

model {
  // likelihood
  //
  n ~ multinomial(x);

  // priors
  //
  p ~ uniform(0.0, 1.0);
  f ~ uniform(-1.0, 1.0);
}
