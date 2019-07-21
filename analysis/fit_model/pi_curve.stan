functions {
  real pi_curve(real b_f, real a, real d_f, real l_f) {
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
  
}
transformed data {
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
  
}
transformed parameters {
  // variable parameters
  for (i in 1:n_events){
    log_b0 = log_b0_mean[event_site[i]] + sig_b0*z_b0[i];
    log_a = log_a_mean[event_site[i]] + sig_a*z_a[i];
    log_d = log_d_mean[event_site[i]] + sig_d*z_d[i];
    log_r = log_r_mean[event_site[i]] + sig_r*z_r[i];
  }
  // exponentiate parameters
  b0 = exp(log_b0); 
  a = exp(log_a);
  d = exp(log_d);
  r = exp(log_r);
  // predicted o2 flux
    for (n in 1:n_obs){
    b[n] = b0[events[n]]*g_b^(temp[n] - temp_ref[n]);
    chi[n] = pi_curve(b[n], a[events[n]], d[events[n]], l[n]);
    kappa[n] = r[events[n]]*g_r^(temp[n] - temp_ref[n]);
    phi[n] = chi[n] - kappa[n];
  }
}
model {
  // priors for fixed values
  g_b ~ normal(1, 1) T[1, ]; 
  g_r ~ normal(1, 1) T[1, ];
  sig_obs ~ normal(0, 1) T[0, ]; 
  // priors for means
  log_b0_mean ~ normal(0, 1); 
  log_a_mean ~ normal(0, 1); 
  log_d_mean ~ normal(0, 1); 
  log_r_mean ~ normal(0, 1);
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
