functions {
  real get_loglik(real inv_temp, vector utilities_untemped, vector utilities_num_untemped) {
     return sum(inv_temp * utilities_num_untemped - log1p_exp(inv_temp * utilities_untemped));
  }
}

data {
  int<lower=1> numChoices;
  vector[numChoices] utilities_untemped;
  vector[numChoices] utilities_num_untemped;
}

parameters {
  real<lower=0> inv_temp;
}

model {
  target += gamma_lpdf(inv_temp | 4, 1) + get_loglik(inv_temp, utilities_untemped, utilities_num_untemped);
}

generated quantities {
  vector[numChoices] log_lik;
  for (n in 1:numChoices) {
    log_lik[n] = inv_temp * utilities_num_untemped[n] - log1p_exp(inv_temp * utilities_untemped[n]);
  }
}
