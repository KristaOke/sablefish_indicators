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
  real incpt;
  real slp[n_trends];
  
  // Estimated trend
  vector[n_trends] real[n_year] pred_trends;
  
  
}


model {
  
  
  
}

