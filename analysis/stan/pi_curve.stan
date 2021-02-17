//=======================================================================================

// Variable photosynthesis-irradiance curves

//=======================================================================================




//=======================================================================================

functions {
  real curve(real b_f, real o_f, real l_f) {
    real p;
    p = b_f * (l_f / o_f) * exp(1 - (l_f / o_f)); // photoinhibition curve
    return p;
  }
}

//=======================================================================================

data {
  // Declare variables
  int n_obs;                      // number of observations
  int x[n_obs, 3];                // parameter groupings
  vector[n_obs] par;              // PAR observations
  vector[n_obs] do_flux;          // DO flux observations
  int n_sum;                      // number of levels over which to generate fitted curves
  int x_sum[n_sum, 3];            // parameter groupings (fitted curves)
  vector[n_sum] par_sum;          // PAR for fitted curves
  
}


//=======================================================================================


transformed data {
  
  // Declare variables
  int xs[3, 3];                   // number of groups per parameter
  real mu;                        // mean observed DO flux
  real tau;                       // sd observed DO flux
  real eta;                       // maximum observed PAR
  real y[n_obs];                  // z-scored DO flux as response variable
  vector[n_obs] l;                // scale observed PAR
  
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
  real lb_m;                      // log beta mean (max GPP)
  real lo_m;                      // log omega mean (optimum light)
  real lr_m;                      // log rho mean (respiration)
  real<lower=0> lb_s[xs[2, 1]];   // log beta sd (max GPP)
  real<lower=0> lo_s[xs[2, 2]];   // log omega sd (optimum light)
  real<lower=0> lr_s[xs[2, 3]];   // log rho sd (respiration)
  vector[xs[3, 1]] lb_d;          // log beta deviate (max GPP)
  vector[xs[3, 2]] lo_d;          // log omega deviate (optimum light)
  vector[xs[3, 3]] lr_d;          // log rho deviate (respiration)
  real<lower=0> s;                // residual standard deviation
  
}


//=======================================================================================


transformed parameters {
  
  // Declare variables 
  vector[xs[1, 1]] b_z;          // beta (maximum GPP when PAR is saturating)
  vector[xs[1, 2]] o_z;          // omega (light level where GPP is optimized)
  vector[xs[1, 3]] r_z;          // rho (ER)
  vector[n_obs] nep_z;           // NEP (on z-scored response scale)
  vector[n_obs] y_p;             // predicted values (on z-scored response scale)
  
  
  // Calculate parameters
  if (xs[1, 1] > 1) {
    b_z = exp(lb_m + lb_s[1] * lb_d);
  } 
  else {
    b_z = exp(rep_vector(lb_m, xs[1, 1]));
  }
  if (xs[1, 2] > 1) {
    o_z = exp(lo_m + lo_s[1] * lo_d);
  } 
  else {
    o_z = exp(rep_vector(lo_m, xs[1, 2]));
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
    real o_t;
    real r_t;
    for(n in 1:n_obs){
      b_t = b_z[x[n, 1]];
      o_t = o_z[x[n, 2]];
      r_t = r_z[x[n, 3]];
      nep_z[n] = curve(b_t, o_t, l[n]) - r_t;
      y_p[n] = nep_z[n] - mu / tau;
  }
  }
  
  
}


//=======================================================================================


model {
  
  // Priors
  lb_m ~ normal(0, 1); 
  lo_m ~ normal(0, 1); 
  lr_m ~ normal(0, 1); 
  lb_s ~ gamma(1.5, 1.5 / 0.5);
  lo_s ~ gamma(1.5, 1.5 / 0.5);
  lr_s ~ gamma(1.5, 1.5 / 0.5);
  s ~ gamma(1.5, 1.5 / 0.5);

  // Random deviates
  lb_d ~ normal(0, 1); 
  lo_d ~ normal(0, 1); 
  lr_d ~ normal(0, 1); 
  
  // Likelihood
  y ~ normal(y_p, s);
  
}


//=======================================================================================


generated quantities {
  
  // Declare variables 
  vector[xs[1, 1]] b;             // beta (maximum GPP when PAR is saturating)
  vector[xs[1, 2]] o;             // omega (light level where GPP is optimized)
  vector[xs[1, 3]] r;             // rho (ER)
  vector[n_obs] nep;              // NEP
  vector[n_sum] gpp_sum;          // GPP (fitted curves)
  vector[n_sum] nep_sum;          // NEP (fitted curves)
  real log_lik [n_obs];           // pointwise log-likelihood
  real log_lik_sum;               // total log-likelihood
  
  
  // Backscale parameters
  b = tau * b_z;
  o = eta * o_z;
  r = tau * r_z;
  
  
  // Backscale metabolism
  nep = tau * nep_z;
  
  // Calculate fitted curves
  {
    real b_t;
    real o_t;
    real r_t;
    for(n in 1:n_sum){
      b_t = b[x_sum[n, 1]];
      o_t = o[x_sum[n, 2]];
      r_t = r[x_sum[n, 3]];
      gpp_sum[n] = curve(b_t, o_t, par_sum[n]);
      nep_sum[n] = gpp_sum[n] - r_t;
    }
  }
  
  // Pointwise log-likelihood
    for (n in 1 : n_obs) {
      log_lik[n] = normal_lpdf(y[n] | y_p[n], s);
    }
    
  // Total log-likelihood
  log_lik_sum = sum(log_lik);

}

//=======================================================================================
