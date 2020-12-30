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
  int<lower=0, upper=1> f_zero;  // 1 = fix it to 0
  int<lower=0, upper=1> f_one;   // 1 = fix it to 1
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

  if ((f_zero == 0) && (f_one == 0)) {
    f = inv_logit(logit_f);
  } else if ((f_zero == 1) && (f_one == 0)) {
    f = 0.0;
  } else if ((f_zero == 0) && (f_one == 1)) {
    f = 1.0;
  } else {
    reject("Inconsistent specification of f_zero and f_one: f_zero=", f_zero,
           ", f_one=", f_one);
  }
  theta = inv_logit(logit_theta);
  pi = inv_logit(logit_pi);

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
      if (N[j,i] > 0) {
        n[j,i] ~ binomial(N[j,i], x[i,j]);
      }
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
  vector[N_pops*N_loci] log_lik; 

  for (i in 1:N_loci) {
    for (j in 1:N_pops) {
      log_lik[(i-1)*N_pops + j] = binomial_lpmf(n[j,i] | N[j,i], x[i,j]);
    }
  }
}
