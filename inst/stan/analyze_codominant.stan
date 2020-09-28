// estimate allele frequencies, inbreeding coefficient, and theta (Fst)
// from co-dominant markers (currently only two alleles per locus)
//

data {
  int<lower=0> N_loci;        // number of loci
  int<lower=0> N_pops;        // number of populations
  int<lower=0> n[N_pops, N_loci, 3]; // genotype counts by locus and pop
  // prior parameters
  //
  real mu_pi;                 // logit of mean(pi)
  real<lower=0> sd_pi;        // sd of logit(pi)
  real mu_f;                  // logit of mean(f)
  real<lower=0> sd_f;         // sd of logit(f)
  real mu_theta;              // logit of mean(theta)
  real<lower=0> sd_theta;     // sd of logit(theta)
  int<lower=0, upper=1> f_zero;  // 0 = estimate f, 1 = fix it to 0
}

parameters {
  real logit_f;               // logit of f
  real logit_theta;           // logit of theta
  vector[N_loci] logit_pi;    // logit of pi
  // allele frequencies by locus & population
  //
  vector<lower=0, upper=1>[N_loci] p[N_pops];
}

transformed parameters {
  real<lower=0, upper=1> f;     // within-population inbreeding coefficient
  real<lower=0, upper=1> theta; // Fst - Weir & Cockerham
  vector<lower=0, upper=1>[N_loci] pi;  // mean allele frequencies
  vector<lower=0, upper=1>[3] x[N_loci, N_pops]; // vector of geno frequencies

  if (f_zero == 0) {
    f = inv_logit(logit_f);
  } else {
    f = 0.0;
  }
  theta = inv_logit(logit_theta);
  pi = inv_logit(logit_pi);

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
  logit_pi ~ normal(mu_pi, sd_pi);
  logit_f ~ normal(mu_f, sd_f);
  logit_theta ~ normal(mu_theta, sd_theta);
  for (i in 1:N_loci) {
    for (j in 1:N_pops) {
      p[j][i] ~ beta(((1.0 - theta)/theta)*pi[i],
                    ((1.0 - theta)/theta)*(1.0 - pi[i]));
    }
  }
}

generated quantities {
  real log_lik[N_pops, N_loci]; 

  for (i in 1:N_loci) {
    for (j in 1:N_pops) {
      log_lik[j,i] = multinomial_lpmf(n[j,i] | x[i,j]);
    }
  }
}
