// estimate allele frequencies, inbreeding coefficient, and theta (Fst)
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
  real mu_theta;              // logit of mean(theta)
  real<lower=0> sd_theta;     // sd of logit(theta)
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
  real<lower=0, upper=1> x[N_loci, N_pops]; // dominant phenotype frequencies

  f = inv_logit(logit_f);
  theta = inv_logit(logit_theta);
  pi = inv_logit(logit_pi);

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
      n[j,i] ~ binomial(N[j,i], x[i,j]);
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
  real log_lik;

  log_lik = 0.0;
  for (i in 1:N_loci) {
    for (j in 1:N_pops) {
      log_lik += binomial_lpmf(n[j,i] | N[j,i], x[i,j]);
    }
  }
}
