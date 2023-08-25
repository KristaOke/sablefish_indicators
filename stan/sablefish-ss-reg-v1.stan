// State-space Regression for Sablefish Indicator

data {
  int<lower=0> n_year;
  vector[n_year] rec_ln;
  vector[n_year] sd_ln;
  
  int<lower=0> n_trends;
  matrix[n_year, n_trends] trends;
  matrix[n_year, n_trends] trends_se;
}


parameters {
  // Regression Parameters
  real incpt;  // Intercept
  real slp[n_trends]; // Slope coefficients
  
  // Estimated trend
  // matrix[n_year, n_trends] pred_trends;
  vector[n_trends] pred_trends[n_year];
  
}

transformed parameters {
  vector[n_year] pred_rec_ln;
  
  real temp_eff;
  temp_eff=0.0;
  
  for(y in 1:n_year) {
    for(t in 1:n_trends) {
      temp_eff += slp[t]*pred_trends[y,t];
    } // next t
    pred_rec_ln[y] = incpt + temp_eff;
  } // next y
  
}


model {
  // PRIORS
  incpt ~ normal(0,10);
  slp ~ normal(0,10);
  
  // LIKELIHOODS
  // Observation Model
  rec_ln ~ normal(pred_rec_ln, sd_ln);
  // Process Model
  for(t in 1:n_trends) {
    for(y in 1:n_year) {
      trends[y,t] ~ normal(pred_trends[y,t], trends_se[y,t]);
    } // next y
  } // next t
  
}

generated quantities {
  
}

