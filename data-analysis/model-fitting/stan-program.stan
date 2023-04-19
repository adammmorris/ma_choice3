functions {
  real get_loglik(real inv_temp, row_vector weights, matrix option_diffs, row_vector choices) {
    row_vector[cols(option_diffs)] utilities = inv_temp * (weights * option_diffs);
    return sum(utilities .* choices - log1p_exp(utilities));
  }
}

data {
  int<lower=1> numChoices;
  int<lower=1> numAtts;       // Number of observations
  matrix[numAtts, numChoices] option_diffs; // Options matrix
  row_vector[numChoices] choices;  // Choices vector (should be 0 or 1)
}

transformed data {
 real uniform_wts_prior = -0.6931472 * numAtts;
}

parameters {
  real<lower=0> inv_temp;
  row_vector<lower=-1, upper=1>[numAtts] weights;
}


model {
  target += gamma_lpdf(inv_temp | 4, 1) + uniform_wts_prior + get_loglik(inv_temp, weights, option_diffs, choices);
}

generated quantities {
  vector[numChoices] log_lik;
  row_vector[numChoices] utilities = inv_temp * (weights * option_diffs);
  for (n in 1:numChoices) {
    log_lik[n] = utilities[n] * choices[n] - log1p_exp(utilities[n]);
  }
}
