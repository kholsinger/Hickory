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
  real<lower=0, upper=1> alpha_l;  // "tightness" of among-locus theta
  real<lower=0, upper=1> alpha_p;  // "tightness" of among-pop theta
}

parameters {
  real logit_f;               // logit of f
  real logit_theta;           // logit of theta
  vector[N_loci] logit_pi;    // logit of pi
  vector<lower=0, upper=1>[N_loci] theta_i; // locus-specific theta
  vector<lower=0, upper=1>[N_pops] theta_j; // population-specific theta
  // allele frequencies by locus & population
  //
  vector<lower=0, upper=1>[N_loci] p[N_pops];
}

transformed parameters {
  real<lower=0, upper=1> f;     // within-population inbreeding coefficient
  real<lower=0, upper=1> theta; // Fst - Weir & Cockerham
  vector<lower=0, upper=1>[N_loci] pi;  // mean allele frequencies
  vector<lower=0, upper=1>[3] x[N_loci, N_pops]; // vector of geno frequencies
  // locus- and population-specific theta
  //
  real<lower=0, upper=1> theta_ij[N_loci, N_pops]; 

  f = inv_logit(logit_f);
  theta = inv_logit(logit_theta);
  pi = inv_logit(logit_pi);

  for (i in 1:N_loci) {
    for (j in 1:N_pops) { 
      theta_ij[i,j] = (theta_i[i] + theta_j[j])/2.0;
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
  for (j in 1:N_pops) {
    theta_j[j] ~ beta(((1.0 - alpha_p)/alpha_p)*theta,
                      ((1.0 - alpha_p)/alpha_p)*(1.0 - theta));
  }
  for (i in 1:N_loci) {
    theta_i[i] ~ beta(((1.0 - alpha_l)/alpha_l)*theta,
                      ((1.0 - alpha_l)/alpha_l)*(1.0 - theta));
  }
  for (i in 1:N_loci) {
    for (j in 1:N_pops) {
      p[j][i] ~ beta(((1.0 - theta_ij[i,j])/theta_ij[i,j])*pi[i],
                     ((1.0 - theta_ij[i,j])/theta_ij[i,j])*(1.0 - pi[i]));
    }
  }
}

generated quantities {
  vector[N_pops*N_loci] log_lik; 

  for (i in 1:N_loci) {
    for (j in 1:N_pops) {
      log_lik[(i-1)*N_pops + j] = multinomial_lpmf(n[j,i] | x[i,j]);
    }
  }
}
