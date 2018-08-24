data {
  // declare data
  int N; // number of observations
  int Nsites; // number of sites
  int Nevents; // number of sampling events (date x site)
  vector<lower=0>[Nsites] site; // site [numeric index]
  vector<lower=0>[Nevents] event; // sampling event [numeric index]
  vector<lower=0>[N] do_flux; // observed DO flux [g m^-2]
  vector<lower=0>[N] par; // observed PAR [photon flux density]
  vector<lower=0>[N] temp; // observed tempreature [C]
  real temp_ref; // reference temperature
}
parameters {
  
}
transformed parameters {
  // structure parameters
  beta0 = beta_site[site] + sig_beta*z_beta[event];
  alpha = alpha_site[site] + sig_alpha*z_alpha[event];
  rho = rho_site[site] + sig_rho*z_rho[event]; 
  // ecosystem metabolism
  beta = beta0*gamma_1^(temp - temp_ref);
  gpp = beta*par/
    ((beta/(alpha*opt^2))*par^2 + 
      (1 - 2*(beta/(alpha*opt)))*par + 
        beta/alpha);
  er = rho*gamma_2^(temp - temp_ref);
  nep = gpp - er;
}
model {
  // priors
  beta_site ~ normal(prior_beta_mean, prior_beta_sig); 
  alpha_site ~ normal(prior_alpha_mean, prior_alpha_sig); 
  rho_site ~ normal(prior_rho_mean, prior_rho_sig);
  sig_beta ~ normal(prior_sig_beta_mean, prior_sig_beta_sig);
  sig_alpha ~ normal(prior_sig_alpha_mean, prior_sig_alpha_sig);
  sig_rho ~ normal(prior_sig_rho_mean, prior_sig_rho_sig);
  // z-values for non-centered paraemterization of random effects
  z_beta ~ normal(0, 1); 
  z_alpha ~ normal(0, 1); 
  z_rho ~ normal(0, 1); 
  // observation model
  do_flux ~ normal(nep, sig_obs);
}
