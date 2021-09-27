data {
  int<lower = 1> n_obs;         // number of days observed
  int<lower = 1> n_weeks;       // number of weeks
  int<lower = 1> n_pop;         // population
  int<lower=0> S0;              // initially Susceptible
  int<lower=0> I0;              // initially Infected
  int<lower=0> yD[n_obs];       // data, total number of new deaths
  int<lower=0> yC[n_obs];       // data, total number of new cases
  int<lower = 1> sigma_cp1;
  int<lower = 1> sigma_cp2;
  int<lower = 1> ifr_cp1;
  int<lower = 1> ifr_cp2;
  int<lower = 1> ifr_cp3;
  real I_D[n_obs];              // discretized infection to death distribution
  real<lower=0> IFRmu[3];
}

transformed data {
  real I_D_rev[n_obs];    // reversed discretized infection to death distribution
  
  for(t in 1:n_obs){
    I_D_rev[t] = I_D[n_obs-t+1];
  }
}

parameters {
  real eta_w0;                     // initial log transmission rate
  real eta_w[n_weeks];             // weekly log transmission rate     
  real<lower = -1,upper=1> rho;    // AR(1) coefficient
  real<lower = 0> reciprocal_phiD;  // overdispersion of NB on deaths
  real<lower=0> sigma[3];         
  real<lower=0> c_tot[n_obs];       // true cases
  real<lower = 0> ifr[3];            // probability of death given infection
  
}  

transformed parameters {
  real eta[n_obs];                     // log transmission rate 
  real<lower = 0> beta[n_obs];         // transmission rate   
  real mu_eta_w[n_weeks];              // mean AR 
  real<lower = 0> phiD;                // 1/reciprocal_phiD
  real<lower = 0> S[n_obs];            // Susceptible individuals
  real<lower = 0> removals[n_obs];     // Removals
  real<lower = 0> I_tot[n_obs];        // Infected individuals
  real<lower = 0, upper=1> r[n_obs];   // probability of Binomial
  real<lower = 0> muD[n_obs];          // mean of poisson on deaths
  
  for (t in 1:n_obs){
    eta[t] = eta_w[1+(t/8)];          // weekly transmission rate
    beta[t] = exp(eta[t]);
  }
  
  mu_eta_w[1] = rho*eta_w0;
  for (t in 2:n_weeks){
    mu_eta_w[t] = rho*eta_w[t-1];
  }
  
  S[1] = S0 - c_tot[1];
  for(t in 2:n_obs){
    S[t]  = S[t-1] - c_tot[t];
  }
  
  for(t in 1:n_obs){
    if (t<7)
      removals[t] = 0;
    else
      removals[t] = c_tot[t-6];
  }
  
  I_tot[1] = I0 + c_tot[1] - removals[1];
  for(t in 2:n_obs){
    I_tot[t] = I_tot[t-1] + c_tot[t] - removals[t];
  }
  
  r[1] = 1-exp(-(beta[1]*I0/n_pop));
  for(t in 2:n_obs){
    r[t] = 1-exp(-(beta[t]*I_tot[t-1]/n_pop));
  }
  
  muD[1] = 0.00001; 
  for (i in 2:n_obs){
    if (i<ifr_cp1)
      muD[i] =  ifr[1] * dot_product(head(c_tot,i-1),tail(I_D_rev, i-1));
    else if (ifr_cp1<=i<ifr_cp2)
      muD[i] =  ifr[2] * dot_product(head(c_tot,i-1),tail(I_D_rev, i-1));
    else if (ifr_cp2<=i<ifr_cp3)
      muD[i] =  ifr[1] * dot_product(head(c_tot,i-1),tail(I_D_rev, i-1));
    else
      muD[i] =  ifr[3] * dot_product(head(c_tot,i-1),tail(I_D_rev, i-1));
  }
  
  phiD = 1. / reciprocal_phiD;
}

model {
  //priors
  rho ~ uniform(-1,1);
  reciprocal_phiD ~ cauchy(0,5);
  for (i in 1:3){
    ifr[i] ~ beta_proportion(IFRmu[i],10000000000);
    sigma[i] ~ cauchy(0,5);
  }
  
  eta_w0 ~ normal(0,1);
  for (t in 1:n_weeks){
    if (t<sigma_cp1)
      eta_w[t] ~ normal(mu_eta_w[t],sigma[1]);
    else if (sigma_cp1<=t<sigma_cp2)
      eta_w[t] ~ normal(mu_eta_w[t],sigma[2]);
    else
      eta_w[t] ~ normal(mu_eta_w[t],sigma[3]);
  }
  
  c_tot[1] ~ gamma(S0*r[1]*(1/(1-r[1])), 1/(1-r[1]));  
  for(t in 2:n_obs){
    c_tot[t] ~ gamma(S[t-1]*r[t]*(1/(1-r[t])), 1/(1-r[t]));
  }
  
  yD ~ neg_binomial_2(muD,phiD);
}

generated quantities {
  
  vector[n_obs] R_t;               // Reproduction number
  real<lower = 0> c_unrep[n_obs];  // Unreported individuals
  vector[n_obs] log_lik;           // log-likelihood for loo package
  real dev;                        // deviance
  
  for (t in 1:n_obs) {
    R_t[t] = beta[t]*6;
  }
  
  for (t in 1:6){ 
    c_unrep[t] = 0;
  }
  
  for (t in 7:n_obs) {
    c_unrep[t] = c_tot[t-6] - yC[t];
    if (c_unrep[t] < 0)
      c_unrep[t] = 0;
  }
  
  for (t in 1:n_obs) {
    log_lik[t] = neg_binomial_2_lpmf(yD[t]| muD[t],phiD);
  }
  
  dev = (-2) * sum(log_lik);
}