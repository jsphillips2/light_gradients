functions {
  real pi_curve(real b_f, real a_f, real d_f, real l_f) {
    real theta1;
    real theta2;
    real theta3;
    real theta4;
    real pred;
    theta1 = b_f*l_f;
    theta2 = (b_f*l_f*l_f)/(a_f*d_f*d_f);
    theta3 = (1 - 2*b_f/(a_f*d_f))*l_f;
    theta4 = b_f/a_f;
    pred = theta1/(theta2 + theta3 + theta4);
    return pred;
  }
}
data {
  int n_obs;
  int n_events;
  int n_sites;
  int event_site[n_events];
  int event[n_obs];
  real do_flux[n_obs]; 
  real<lower=0> temp[n_obs];
  real<lower=0> par[n_obs];
  real<lower=0> temp_ref;
}
transformed data {
  real mu;
  real<lower=0> tau;
  real<lower=0> eta;
  real x[n_obs];
  real<lower=0> l[n_obs];
  // scale data
  mu = mean(do_flux);
  tau = sd(do_flux);
  eta = mean(par);
  for (n in 1:n_obs){
    x[n] = (do_flux[n] - mu)/tau;
    l[n] = par[n]/eta;
  }
}
parameters {
  // priors for fixed values
  real<lower=1> g_b; 
  real<lower=1> g_r;
  real<lower=0> sig_obs; 
  // priors for means
  real log_b0_mean[n_sites];
  real log_a_mean[n_sites];
  real log_d_mean[n_sites];
  real log_r_mean[n_sites];
  // priors for sd's
  real<lower=0> sig_b0;
  real<lower=0> sig_a;
  real<lower=0> sig_d;
  real<lower=0> sig_r;
  // z's
  real z_b0[n_events];
  real z_a[n_events];
  real z_d[n_events];
  real z_r[n_events]; 
}
transformed parameters {
  real log_b0[n_events];
  real log_a[n_events];
  real log_d[n_events];
  real log_r[n_events];
  real b0[n_events];
  real a[n_events];
  real d[n_events];
  real r[n_events];
  real b[n_obs];
  real chi[n_obs];
  real kappa[n_obs];
  real phi[n_obs];
  // variable parameters
  for (i in 1:n_events){
    log_b0[i] = sig_b0*z_b0[i];
    log_a[i] = sig_a*z_a[i];
    log_d[i] = sig_d*z_d[i];
    log_r[i] = sig_r*z_r[i];
  }
  // exponentiate parameters
  b0 = exp(log_b0); 
  a = exp(log_a);
  d = exp(log_d);
  r = exp(log_r);
  // predicted o2 flux
    for (n in 1:n_obs){
    b[n] = b0[event[n]]*g_b^(temp[n] - temp_ref);
    chi[n] = pi_curve(b[n], a[event[n]], d[event[n]], l[n]);
    kappa[n] = r[event[n]]*g_r^(temp[n] - temp_ref);
    phi[n] = chi[n] - kappa[n];
  }
}
model {
  // priors for fixed values
  g_b ~ normal(1, 1) T[1, ]; 
  g_r ~ normal(1, 1) T[1, ];
  sig_obs ~ normal(0, 1) T[0, ]; 
  // priors for sd's
  sig_b0 ~ normal(0, 1) T[0, ]; 
  sig_a ~ normal(0, 1) T[0, ];
  sig_d ~ normal(0, 1) T[0, ];
  sig_r ~ normal(0, 1) T[0, ];
  // z's
  z_b0 ~ normal(0, 1); 
  z_a ~ normal(0, 1); 
  z_d ~ normal(0, 1); 
  z_r ~ normal(0, 1);
  // observation process
  x ~ normal(phi, sig_obs);
}
generated quantities {
  real beta0[n_events];
  real alpha[n_events];
  real delta[n_events];
  real rho[n_events];
  real beta[n_obs];
  real gpp[n_obs];
  real er[n_obs];
  real nep[n_obs];
  real log_lik[n_obs];
  for(i in 1:n_events){
    beta0[i] = tau*b0[i];
    alpha[i] = tau*a[i]/eta;
    delta[i] = eta*d[i];
    rho[i] = tau*r[i];
  }
  for(n in 1:n_obs){
    beta[n] = tau*b[n];
    gpp[n] = tau*chi[n];
    er[n] = tau*kappa[n];
    nep[n] = tau*phi[n];
    log_lik[n] = normal_lpdf(x[n]|phi[n], sig_obs);
  }
}