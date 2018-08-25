data {
  // declare data
  // indices
  int N; // number of observations
  int Nsites; // number of sites
  int Nevents; // number of sampling events (date x site)
  // actual data
  int event_sites[Nevents]; // site for each event [numeric index]
  int event[N]; // sampling event [numeric index]
  vector[N] do_flux; // observed DO flux [g m^-2]
  vector[N] par; // observed PAR [photon flux density]
  vector[N] temp; // observed tempreature [C]
  real temp_ref; // reference temperature [C]
  // priors
  vector[2] prior_beta_site;
  vector[2] prior_alpha_site;
  vector[2] prior_opt_site;
  vector[2] prior_rho_site;
  vector[2] prior_sig_beta;
  vector[2] prior_sig_alpha;
  vector[2] prior_sig_opt;
  vector[2] prior_sig_rho;
  vector[2] prior_gamma_1;
  vector[2] prior_gamma_2;
  vector[2] prior_sig_obs;
}
parameters {
  // delcare variables
  vector[Nsites] beta_site;
  vector[Nsites] alpha_site;
  vector[Nsites] opt_site;
  vector[Nsites] rho_site;
  real sig_beta;
  real sig_alpha;
  real sig_opt;
  real sig_rho;
  real gamma_1;
  real gamma_2;
  vector[Nevents] z_beta;
  vector[Nevents] z_alpha;
  vector[Nevents] z_opt;
  vector[Nevents] z_rho;
  real<lower=0> sig_obs;
}
transformed parameters {
  // declar variables
  vector[Nevents] beta0;
  vector[N] beta;
  vector[Nevents] alpha;
  vector[Nevents] opt;
  vector[Nevents] rho;
  vector[N] gpp;
  vector[N] er;
  vector[N] nep;
  // structure parameters
  beta0 = beta_site[event_sites[1:Nevents]] + sig_beta*z_beta[1:Nevents];
  alpha = alpha_site[event_sites[1:Nevents]] + sig_alpha*z_alpha[1:Nevents];
  opt = opt_site[event_sites[1:Nevents]] + sig_opt*z_opt[1:Nevents];
  rho = rho_site[event_sites[1:Nevents]] + sig_rho*z_rho[1:Nevents]; 
  // ecosystem metabolism
  for(n in 1:N){
    beta[n] = beta0[event[n]]*gamma_1^(temp[n] - temp_ref);
    gpp[n] = beta[n]*par[n]/
      ((beta[n]/(alpha[event[n]]*opt[event[n]]^2))*par[n]^2 + 
        (1 - 2*(beta[n]/(alpha[event[n]]*opt[event[n]])))*par[n] + 
          beta[n]/alpha[event[n]]);
    er[n] = rho[event[n]]*gamma_2^(temp[n] - temp_ref);
    nep[n] = gpp[n] - er[n];
  }
}
model {
  // priors
  for(s in 1:Nsites){
    beta_site[s] ~ normal(prior_beta_site[1], prior_beta_site[2]) T[0,]; 
    alpha_site[s] ~ normal(prior_alpha_site[1], prior_alpha_site[2]) T[0,]; 
    opt_site[s] ~ normal(prior_opt_site[1], prior_opt_site[2]) T[0,]; 
    rho_site[s] ~ normal(prior_rho_site[1], prior_rho_site[2]) T[0,];
  }
  sig_beta ~ normal(prior_sig_beta[1], prior_sig_beta[2]) T[0,];
  sig_alpha ~ normal(prior_sig_alpha[1], prior_sig_alpha[2]) T[0,];
  sig_opt ~ normal(prior_sig_opt[1], prior_sig_opt[2]) T[0,];
  sig_rho ~ normal(prior_sig_rho[1], prior_sig_rho[2]) T[0,];
  gamma_1 ~ normal(prior_gamma_1[1], prior_gamma_1[2]) T[1,];
  gamma_2 ~ normal(prior_gamma_2[1], prior_gamma_2[2]) T[1,];
  sig_obs ~ normal(prior_sig_obs[1], prior_sig_obs[2]) T[0,];
  // z-values for non-centered paraemterization of random effects
  z_beta ~ normal(0, 1); 
  z_alpha ~ normal(0, 1); 
  z_opt ~ normal(0, 1); 
  z_rho ~ normal(0, 1); 
  // observation model
  do_flux ~ normal(nep, sig_obs);
}
