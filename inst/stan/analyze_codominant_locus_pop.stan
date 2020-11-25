// estimate allele frequencies, inbreeding coefficient, and theta (Fst)
// from co-dominant markers (currently only two alleles per locus)
//

functions {
  real beta_scale(real p) {
    real scale;
    
    if (p < 0.5) {
      scale = 1.0/p;
    } else {
      scale = 1.0/(1.0 - p);
    }
    return scale;
  }
}

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
  int<lower=0, upper=1> f_zero;  // 1 = fix it to 0
  int<lower=0, upper=1> f_one;   // 1 = fix it to 1
  real<lower=0, upper=1> alpha_l;  // "tightness" of among-locus theta
  real<lower=0, upper=1> alpha_p;  // "tightness" of among-pop theta
}

parameters {
  real logit_f;               // logit of f
  real logit_theta;           // logit of theta
  real<lower=0, upper=1> alpha_ll;  // "tightness" of among-locus theta
  real<lower=0, upper=1> alpha_pp;  // "tightness" of among-pop theta
  vector[N_loci] logit_pi;    // logit of pi
  vector<lower=0, upper=1>[N_loci] theta_l; // locus-specific theta
  vector<lower=0, upper=1>[N_pops] theta_p; // population-specific theta
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
  real<lower=0, upper=1> theta_lp[N_loci, N_pops]; 

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
      theta_lp[i,j] = (theta_l[i] + theta_p[j])/2.0;
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
  alpha_pp ~ beta(beta_scale(alpha_p)*alpha_p,
                  beta_scale(alpha_p)*(1.0 - alpha_p));
  alpha_ll ~ beta(beta_scale(alpha_l)*alpha_l,
                  beta_scale(alpha_l)*(1.0 - alpha_l));
  for (j in 1:N_pops) {
    theta_p[j] ~ beta(((1.0 - alpha_pp)/alpha_pp)*theta,
                      ((1.0 - alpha_pp)/alpha_pp)*(1.0 - theta));
  }
  for (i in 1:N_loci) {
    theta_l[i] ~ beta(((1.0 - alpha_ll)/alpha_ll)*theta,
                      ((1.0 - alpha_ll)/alpha_ll)*(1.0 - theta));
  }
  for (i in 1:N_loci) {
    for (j in 1:N_pops) {
      p[j][i] ~ beta(((1.0 - theta_lp[i,j])/theta_lp[i,j])*pi[i],
                     ((1.0 - theta_lp[i,j])/theta_lp[i,j])*(1.0 - pi[i]));
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
