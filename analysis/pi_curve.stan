//=======================================================================================

// Variable PI curves

//=======================================================================================




//=======================================================================================

data {
  // Declare variables
  int n_obs; // number of observations
  int x[n_obs, 3]; // parameter groupings
  vector[n_obs] par; // PAR observations
  vector[n_obs] tscale;
  vector[n_obs] do_flux; // DO flux observations
  int n_sum; // number of levels over which to generate fitted curves
  int x_sum[n_sum, 3]; // parameter groupings (fitted curves)
  vector[n_sum] par_sum; // PAR for fitted curves
  
}


//=======================================================================================


transformed data {
  
  // Declare variables
  int xs[3, 3]; // number of groups per parameter
  real mu; // mean observed DO flux
  real tau; // sd observed DO flux
  real eta; // maximum observed PAR
  real y[n_obs]; // z-scored DO flux as response variable
  vector[n_obs] l; // scale observed PAR
  
  // Number of groups per parameter
  for (i in 1 : 3) {
    xs[1, i] = max(x[, i]);
    if(xs[1, i] > 1) xs[2, i] = 1; else xs[2, i] = 0;
    if(xs[1, i] > 1) xs[3, i] = xs[1, i]; else xs[3, i] = 0;
  }
  
  // Calculate values for scaling
  mu = mean(do_flux);
  tau = sd(do_flux);
  eta = mean(par);
  
  // Scale variables
  for (n in 1 : n_obs){
    y[n] = (do_flux[n] - mu)/tau;
    l[n] = par[n]/eta;
  }
  
}


//=======================================================================================


parameters {
  
  // Declare variables
  real lb_m;
  real la_m;
  real lr_m;
  real<lower=0> lb_s[xs[2, 1]];
  real<lower=0> la_s[xs[2, 2]];
  real<lower=0> lr_s[xs[2, 3]];
  vector[xs[3, 1]] lb_d;
  vector[xs[3, 2]] la_d;
  vector[xs[3, 3]] lr_d;
  real<lower=0> e_b;
  real<lower=0> e_r;
  real<lower=0> s; // residual standard deviation
  
}


//=======================================================================================


transformed parameters {
  
  // Declare variables 
  vector[xs[1, 1]] b_z; // beta (maximum GPP when PAR is saturating)
  vector[xs[1, 2]] a_z; // alpha (increase in GPP with PAR when PAR is limiting)
  vector[xs[1, 3]] r_z; // rho (ER)
  vector[n_obs] nep_z; // NEP (on z-scored response scale)
  vector[n_obs] y_p; // predicted values (on z-scored response scale)
  
  
  // Calculate parameters
  if (xs[1, 1] > 1) {
    b_z = exp(lb_m + lb_s[1] * lb_d);
  } 
  else {
    b_z = exp(rep_vector(lb_m, xs[1, 1]));
  }
  if (xs[1, 2] > 1) {
    a_z = exp(la_m + la_s[1] * la_d);
  } 
  else {
    a_z = exp(rep_vector(la_m, xs[1, 2]));
  }
  if (xs[1, 3] > 1) {
    r_z = exp(lr_m + lr_s[1] * lr_d);
  } 
  else {
    r_z = exp(rep_vector(lr_m, xs[1, 3]));
  }
  
  // Calculate metablism
  {
    real b_t;
    real a_t;
    real r_t;
    for(n in 1:n_obs){
      b_t = b_z[x[n, 1]] * exp(-e_b * tscale[n]);
      a_t = a_z[x[n, 2]];
      r_t = r_z[x[n, 3]] * exp(-e_r * tscale[n]);
      nep_z[n] = b_t * tanh((a_t / b_t) * l[n]) - r_t;
      y_p[n] = nep_z[n] - mu / tau;
  }
  }
  
  
}


//=======================================================================================


model {
  
  // Priors
  lb_m ~ normal(0, 1); 
  la_m ~ normal(0, 1); 
  lr_m ~ normal(0, 1); 
  lb_s ~ gamma(1.5, 1.5 / 0.5);
  la_s ~ gamma(1.5, 1.5 / 0.5);
  lr_s ~ gamma(1.5, 1.5 / 0.5);
  s ~ gamma(1.5, 1.5 / 0.5);
  e_b ~ exponential(1);
  e_r ~ exponential(1);
  
  // Random deviates
  lb_d ~ normal(0, 1); 
  la_d ~ normal(0, 1); 
  lr_d ~ normal(0, 1); 
  
  
  // Likelihood
  y ~ normal(y_p, s);
  
}


//=======================================================================================


generated quantities {
  
  // Declare variables 
  vector[xs[1, 1]] b; // beta (maximum GPP when PAR is saturating)
  vector[xs[1, 2]] a; // alpha (increase in GPP with PAR when PAR is limiting)
  vector[xs[1, 3]] r; // rho (ER)
  vector[n_obs] nep; // NEP
  vector[n_sum] gpp_sum; // GPP (fitted curves)
  vector[n_sum] nep_sum; // NEP (fitted curves)
  vector[n_sum] sat; // half-saturation constant
  real log_lik [n_obs]; // pointwise log-likelihood
  real log_lik_sum; // total log-likelihood
  
  
  // Backscale parameters
  b = tau * b_z;
  a = tau * a_z / eta;
  r = tau * r_z;
  
  
  // Backscale metabolism
  nep = tau * nep_z;
  
  // Calculate fitted curves
  for(n in 1:n_sum){
    gpp_sum[n] = b[x_sum[n, 1]] * tanh((a[x_sum[n, 2]] / b[x_sum[n, 1]]) * par_sum[n]);
    nep_sum[n] = gpp_sum[n] - r[x_sum[n, 3]];
    sat[n] = (b[x_sum[n, 1]] / a[x_sum[n, 2]]) * atanh(0.5);
  }
  
  // Pointwise log-likelihood
    for (n in 1 : n_obs) {
      log_lik[n] = normal_lpdf(y[n] | y_p[n], s);
    }
    
  // Total log-likelihood
  log_lik_sum = sum(log_lik);

}

//=======================================================================================
