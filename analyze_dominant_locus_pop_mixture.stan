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
  real<lower=0, upper=1> alpha;  // "tightness" of among-locus theta
}

parameters {
  real logit_f;               // logit of f
  real logit_theta;           // logit of theta
  vector[N_loci] logit_pi;    // logit of pi
  vector<lower=0, upper=1>[N_loci] theta_i; // locus-specific theta
  // allele frequencies by locus & population
  //
  vector<lower=0, upper=1>[N_loci] p[N_pops];
}

transformed parameters {
  real<lower=0, upper=1> f;     // within-population inbreeding coefficient
  real<lower=0, upper=1> theta; // Fst - Weir & Cockerham
  vector<lower=0, upper=1>[N_loci] pi;  // mean allele frequencies
  real<lower=0, upper=1> x[N_loci, N_pops]; // dominant phenotype frequencies
  real<lower=0, upper=1> hw[N_loci, N_pops]; // dominant frequency under H-W

  f = inv_logit(logit_f);
  theta = inv_logit(logit_theta);
  pi = inv_logit(logit_pi);

  for (i in 1:N_loci) {
    for (j in 1:N_pops) {
      x[i,j] = (p[j][i]^2)*(1.0 - f) + f*p[j][i] +
               2.0*p[j][i]*(1.0 - p[j][i])*(1-f);
      hw[i,j] = p[j][i]^2 + 2.0*p[j][i]*(1.0 - p[j][i]);
    }
  }
}

model {
  // likelihood
  //
  for (i in 1:N_loci) {
    for (j in 1:N_pops) {
      if (N[j,i] > 0) {
        // n[j,i] ~ binomial(N[j,i], x[i,j]);
        target += log_sum_exp(log(1.0 - f) + binomial_lpmf(n[j,i] | N[j,i],
                                                                    hw[i,j]),
                              log(f) + binomial_lpmf(n[j,i] | N[j,i],
                                                              p[j][i]));
      }
    }
  }

  // priors
  //
  logit_pi ~ normal(mu_pi, sd_pi);
  logit_f ~ normal(mu_f, sd_f);
  logit_theta ~ normal(mu_theta, sd_theta);
  for (i in 1:N_loci) {
    theta_i[i] ~ beta(((1.0 - alpha)/alpha)*theta,
                      ((1.0 - alpha)/alpha)*(1.0 - theta));
    for (j in 1:N_pops) {
      p[j][i] ~ beta(((1.0 - theta_i[i])/theta_i[i])*pi[i],
                     ((1.0 - theta_i[i])/theta_i[i])*(1.0 - pi[i]));
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
