// estimate allele frequency and inbreeding coefficient at one locus
//

data {
  int<lower = 0> n[3];          // vector of genotype counts
}

transformed data {
  vector<lower=0>[3] alpha;

  for (i in 1:3) {
    alpha[i] = 1.0;
  }
}

parameters {
  simplex[3] x;                // vector of genotype frequencies
}

transformed parameters {
  real<lower=0, upper=1> p;    // allele frequency
  real<lower=-1, upper=1> f;   // inbreeding coefficient

  p = x[1] + x[2]/2.0;
  f = 1.0 - x[2]/(2.0*p*(1.0 - p));
}

model {
  // likelihood
  //
  n ~ multinomial(x);

  // prior
  //
  x ~ dirichlet(alpha);
}
