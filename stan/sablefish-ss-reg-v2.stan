// State-space Regression for Sablefish Indicator - v2
//   Adds capacity to generate predictions for new DFA trend value(s) and their SE's

data {
  int<lower=0> n_year;
  vector[n_year] rec_ln;
  vector[n_year] sd_ln;
  
  int<lower=0> n_trends;
  matrix[n_year, n_trends] trends;
  matrix[n_year, n_trends] trends_se;
  
  // // Prediction quantities
  int<lower=0> n_year_fcst;
  matrix[n_year_fcst, n_trends] trends_fcst;
  matrix[n_year_fcst, n_trends] trends_se_fcst;

}


parameters {
  // Regression Parameters
  real incpt;  // Intercept
  real slp[n_trends]; // Slope coefficients
  // vector[n_trends] slp;
  
  // Estimated trend
  matrix[n_year, n_trends] pred_trends;
  // vector[n_trends] pred_trends[n_year];
  
}

transformed parameters {
  vector[n_year] pred_rec_ln;
  
  vector[n_year] temp_eff; // Temporary additive effect of DFA trends
  
  
  for(y in 1:n_year) {
    temp_eff[y]=0.0;
    for(t in 1:n_trends) {
      temp_eff[y] += slp[t]*pred_trends[y,t];
    } // next t
    pred_rec_ln[y] = incpt + temp_eff[y];
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
  // Forecasted Quantities   
  vector[n_year_fcst] pred_rec_ln_fcst;
  vector[n_year_fcst] temp_eff_fcst; // Temporary additive effect of DFA trends
  // Posterior predicted quantities 
  // vector[n_trends] pred_trends_fcst[n_year_fcst];
  matrix[n_year_fcst, n_trends] pred_trends_fcst;
  vector[n_year_fcst] pred_rec_ln_fcst_postPred;
  vector[n_year_fcst] temp_eff_fcst_postPred; // Temporary additive effect of DFA trends

  
  // Randomly sample value for forecasted trends   
  // Process Model
  for(t in 1:n_trends) {
    for(y in 1:n_year_fcst) {
      // trends[y,t] ~ normal(pred_trends[y,t], trends_se[y,t]);
      pred_trends_fcst[y,t] = normal_rng(trends_fcst[y,t], trends_se_fcst[y,t]);
    } // next y
  } // next t
  
  for(y in 1:n_year_fcst) {
    temp_eff_fcst[y]=0.0;
    temp_eff_fcst_postPred[y]=0.0;
    for(t in 1:n_trends) {
      temp_eff_fcst_postPred[y] += slp[t]*pred_trends_fcst[y,t];
      temp_eff_fcst[y] += slp[t]*trends_fcst[y,t];
    } // next t
    pred_rec_ln_fcst_postPred[y] = incpt + temp_eff_fcst_postPred[y];
    pred_rec_ln_fcst[y] = incpt + temp_eff_fcst[y];
  } // next y
}

