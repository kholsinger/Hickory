// estimate allele frequencies and inbreeding coefficient
// from dominant markers
//

data {
  int<lower=0> N_loci;        // number of loci
  int<lower=0> N_pops;        // number of populations
  int<lower=0> n[N_pops, N_loci]; // counts of dominant phenotypes
  int<lower=0> N[N_pops, N_loci]; // counts of sample sizes
  // prior parameters
  //
  real mu_pi;                 // logit of mean(pi)
  real<lower=0> sd_pi;        // sd of logit(pi)
  real mu_f;                  // logit of mean(f)
  real<lower=0> sd_f;         // sd of logit(f)
}

parameters {
  real logit_f;               // logit of f
  vector[N_loci] logit_pi;    // logit of pi
}

transformed parameters {
  real<lower=0, upper=1> f;     // within-population inbreeding coefficient
  vector<lower=0, upper=1>[N_loci] pi;  // mean allele frequencies
  // allele frequencies by locus & population
  //
  vector<lower=0, upper=1>[N_loci] p[N_pops];
  // dominant phenotype frequencies
  //
  real<lower=0, upper=1> x[N_loci, N_pops]; 

  f = inv_logit(logit_f);
  pi = inv_logit(logit_pi);

  for (i in 1:N_pops) {
    p[i] = pi;
  }

  for (i in 1:N_loci) {
    for (j in 1:N_pops) {
      // x[i,j] = (p[j][i]^2)*(1.0 - f) + f*p[j][i] +
      //          2.0*p[j][i]*(1.0 - p[j][i])*(1-f);
      x[i, j] = 1.0 - (((1.0 - p[j][i])^2)*(1.0 - f) +
                f*(1.0 - p[j][i]));
    }
  }
}

model {
  // likelihood
  //
  for (i in 1:N_loci) {
    for (j in 1:N_pops) {
      n[j,i] ~ binomial(N[j,i], x[i,j]);
    }
  }

  // priors
  //
  logit_pi ~ normal(mu_pi, sd_pi);
  logit_f ~ normal(mu_f, sd_f);
}

generated quantities {
  vector[N_pops*N_loci] log_lik; 

  for (i in 1:N_loci) {
    for (j in 1:N_pops) {
      log_lik[(i-1)*N_pops + j] = binomial_lpmf(n[j,i] | N[j,i], x[i,j]);
    }
  }
}
