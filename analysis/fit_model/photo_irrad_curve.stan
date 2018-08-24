data {
  // declare data
  // indices
  int N; // number of observations
  int Nsites; // number of sites
  int Nevents; // number of sampling events (date x site)
  // actual data
  vector<lower=0>[Nsites] site; // site [numeric index]
  vector<lower=0>[Nevents] event; // sampling event [numeric index]
  vector<lower=0>[N] do_flux; // observed DO flux [g m^-2]
  vector<lower=0>[N] par; // observed PAR [photon flux density]
  vector<lower=0>[N] temp; // observed tempreature [C]
  real temp_ref; // reference temperature [C]
  // priors: means
  real<lower=0> prior_beta_site;
  real<lower=0> prior_alpha_site;
  real<lower=0> prior_opt_site;
  real<lower=0> prior_rho_site;
  real<lower=0> prior_sig_beta;
  real<lower=0> prior_sig_alpha;
  real<lower=0> prior_sig_opt;
  real<lower=0> prior_sig_rho;
  real<lower=1> prior_gamma_1;
  real<lower=1> prior_gamma_2;
  real<lower=0> prior_sig_obs;
  // priors: standard deviations
  real<lower=0> prior_beta_site_sig;
  real<lower=0> prior_alpha_site_sig;
  real<lower=0> prior_opt_site_sig;
  real<lower=0> prior_rho_site_sig;
  real<lower=0> prior_sig_beta_sig;
  real<lower=0> prior_sig_alpha_sig;
  real<lower=0> prior_sig_opt_sig;
  real<lower=0> prior_sig_rho_sig;
  real<lower=1> prior_gamma_1_sig;
  real<lower=1> prior_gamma_2_sig;
  real<lower=0> prior_sig_obs_sig;
}
parameters {
  // delcare variables
  vector<lower=0>[Nsites] beta_site;
  vector<lower=0>[Nsites] alpha_site;
  vector<lower=0>[Nsites] opt_site;
  vector<lower=0>[Nsites] rho_site;
  real<lower=0> sig_beta;
  real<lower=0> sig_alpha;
  real<lower=0> sig_opt;
  real<lower=0> sig_rho;
  real<lower=1> gamma_1;
  real<lower=1> gamma_2;
  vector[Nevents] z_beta;
  vector[Nevents] z_alpha;
  vector[Nevents] z_opt;
  vector[Nevents] z_rho;
  real<lower=0> sig_obs;
}
transformed parameters {
  // declar variables
  vector[N] beta0;
  vector[N] beta;
  vector[N] alpha;
  vector[N] opt;
  vector[N] rho;
  vector[N] gpp;
  vector[N] er;
  vector[N] nep;
  // structure parameters
  beta0 = beta_site[site] + sig_beta*z_beta[event];
  alpha = alpha_site[site] + sig_alpha*z_alpha[event];
  opt = opt_site[site] + sig_opt*z_opt[event];
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
  beta_site ~ normal(prior_beta, prior_beta_sig) T[0,]; 
  alpha_site ~ normal(prior_alpha, prior_alpha_sig) T[0,]; 
  opt_site ~ normal(prior_opt, prior_opt_sig) T[0,]; 
  rho_site ~ normal(prior_rho, prior_rho_sig) T[0,];
  sig_beta ~ normal(prior_sig_beta, prior_sig_beta_sig) T[0,];
  sig_alpha ~ normal(prior_sig_alpha, prior_sig_alpha_sig) T[0,];
  sig_opt ~ normal(prior_sig_opt, prior_sig_opt_sig) T[0,];
  sig_rho ~ normal(prior_sig_rho, prior_sig_rho_sig) T[0,];
  gamma_1 ~ normal(prior_gamma_1, prior_gamma_1_sig) T[1,];
  gamma_2 ~ normal(prior_gamma_2, prior_gamma_2_sig) T[1,];
  sig_obs ~ normal(prior_sig_obs, prior_sig_obs_sig) T[0,];
  // z-values for non-centered paraemterization of random effects
  z_beta ~ normal(0, 1); 
  z_alpha ~ normal(0, 1); 
  z_opt ~ normal(0, 1); 
  z_rho ~ normal(0, 1); 
  // observation model
  do_flux ~ normal(nep, sig_obs);
}
