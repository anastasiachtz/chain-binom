data {
  int<lower = 1> n_obs;         // number of days observed
  int<lower = 1> n_groups;      // number of age groups 
  int<lower = 1> ifr_cp1;       // change point ifr
  int<lower = 1> N_g[n_groups]; // population of each age group          
  int yD[n_obs];                // data, total number of new deaths
  int<lower=0> S0[n_groups];    // initially susceptible individuals
  vector<lower=0>[n_groups] I0; // initially infected individuals
  matrix<lower=0>[n_groups,n_groups] C_M;  // mean contact rates
  real I_D[n_obs];              // discretized infection to death distribution
  real IFRmu[2];                // mean  ifr
  real v0;                      // initial time-varying factor
}

transformed data {
  real x_r[0];
  int x_i[0];
  real I_D_rev[n_obs];    // reversed discretized infection to death distribution
  
  for(t in 1:n_obs){
    I_D_rev[t] = I_D[n_obs-t+1];
  }
}

parameters {
  real<lower = 0, upper=1> p;        // probability of infection given contact
  matrix<lower = 0>[n_groups, n_groups] mu;    // Contact Matrix
  real v[n_obs];                     // time-varying factor (OU)  
  real<lower = 0> sigma_mu;          // sd of contacts
  real<lower = 0> reciprocal_phiD;   // overdispersion of NB on deaths
  real<lower = 0> reciprocal_phi;    // 1/speed of reversion
  real<lower = 0> sigma;             // instantaneous diffusion term
  real<lower = 0> c_i[n_groups,n_obs];     // new cases per age group
  real<lower = 0> raw_c_i[n_groups,n_obs]; // raw new cases per age group
  real<lower = 0> ifr[2];            // probability of death given infection
}  

transformed parameters {
  real<lower = 0> w[n_obs];             // time-varying factor 
  real mu_v[n_obs];                     // mean of OU
  real<lower = 0> phiD;                 // 1/reciprocal_phiD
  real<lower = 0> phi;                  // 1/reciprocal_phi
  real<lower = 0> k[n_groups];          // p/N_g
  real<lower = 0> sigmaOU;              // sd of OU
  real<lower = 0> S[n_groups,n_obs];    // susceptible per age group
  real<lower = 0> removals[n_groups,n_obs]; // removals per age group
  matrix<lower = 0>[n_groups,n_obs] I_tot;  // infected per age group
  real<lower = 0> M[n_groups,n_obs];     // sum(C*I)
  real<lower = 0> pr[n_groups,n_obs];    // k*w*M
  real<lower = 0, upper=1> r[n_groups,n_obs]; // probability of CB
  real<lower = 0> c_tot[n_obs];         // total cases
  real<lower = 0> muD[n_obs];           // mean of NB on deaths
  
  phi = 1. / reciprocal_phi;
  phiD = 1. / reciprocal_phiD;
  
  w = exp(v);
  
  mu_v[1] = v0*exp(-phi);
  for (t in 2:n_obs){
    mu_v[t] = v[t-1]*exp(-phi);
  }
  
  for(i in 1:n_groups){
    k[i] = p/N_g[i];
  }
  
  sigmaOU = ((sigma^2)/(2*phi))*(1- exp(-2*phi));
  
  for (i in 1:n_groups){
    S[i,1] = S0[i] - c_i[i,1];
  }
  
  for (i in 1:n_groups){
    for(t in 2:n_obs){
      S[i,t]  = S[i,t-1] - c_i[i,t];
    }
  }
  
  for (i in 1:n_groups){
    for(t in 1:n_obs){
      if (t<7)
        removals[i,t] = 0;
      else
        removals[i,t] = c_i[i,t-6];
    }
  }
  for (i in 1:n_groups){
    I_tot[i,1] = I0[i] + c_i[i,1] - removals[i,1];
  }
  
  for (i in 1:n_groups){
    for(t in 2:n_obs){
      I_tot[i,t] = I_tot[i,t-1] + c_i[i,t] - removals[i,t];
    }
  }
  
  for (i in 1:n_groups){
    M[i,1] = dot_product(mu[,i],I0);
    pr[i,1] = k[i]*w[1]*M[i,1];
    r[i,1] = 1-exp(-pr[i,1]);
  }
  
  for (i in 1:n_groups){
    for(t in 2:n_obs){
      M[i,t] = dot_product(mu[,i],I_tot[,t-1]);
      pr[i,t] = k[i]*w[t]*M[i,t];
      r[i,t] = 1-exp(-pr[i,t]);
    }
  }
  
  for (t in 1:n_obs){
    c_tot[t] = sum(raw_c_i[,t]);
  }
  
  muD[1] = 3; 
  for (i in 2:n_obs){
    if (i<ifr_cp1)
      muD[i] =  ifr[1] * dot_product(head(c_tot,i-1),tail(I_D_rev, i-1));
    else
      muD[i] =  ifr[2] * dot_product(head(c_tot,i-1),tail(I_D_rev, i-1));
  }
  
}

model {
  //priors
  p ~ beta(1,1);
  for (i in 1:2){
    ifr[i] ~ beta_proportion(IFRmu[i],10000000000);
  }
  reciprocal_phiD ~ cauchy(0,5);
  reciprocal_phi ~ cauchy(0,5);
  sigma ~ cauchy(0,5);
  sigma_mu ~ cauchy(0,5);
  
  for (t in 1:n_obs){
    v[t]~normal(mu_v[t],sigmaOU);
  }  
  
  for(i in 1:n_groups){
    for(j in 1:n_groups){
      mu[i,j] ~ normal(C_M[i,j],sigma_mu);
    }}
  
  
  for (i in 1:n_groups){
    c_i[i,1] ~ gamma(S0[i]*r[i,1]*(1/(1-r[i,1])), 1/(1-r[i,1]));
  }
  
  for (i in 1:n_groups){
    for(t in 2:n_obs){
      c_i[i,t] ~ gamma(S[i,t-1]*r[i,t]*(1/(1-r[i,t])), 1/(1-r[i,t]));
    }
  }
  
  for (t in 1:n_obs){
    for(i in 1:n_groups){
      //raw_c_i[i,t] ~ normal(c_i[i,t],sqrt(c_i[i,t]));
      raw_c_i[i,t] ~ gamma(c_i[i,t],1);
    }
  }
  
  yD ~ neg_binomial_2(muD,phiD);
}