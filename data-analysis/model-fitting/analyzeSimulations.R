# analyze results ---------------------------------------------------------
numGeneratingModels = 6

### import data
fits_parsed = vector('list', numGeneratingModels)
lmes_parsed = vector('list', numGeneratingModels)
true_temps = vector('list', numGeneratingModels)
true_weights = vector('list', numGeneratingModels)
best_fit_temps = vector('list', numGeneratingModels)
best_fit_weights = vector('list', numGeneratingModels)
for (i in 1:numGeneratingModels) {
  load(paste0('simulations/', i, '/analysis_output_stan.rdata'))
  
  fits_parsed[[i]] = list(stan_results_full_parsed[[1]],
                          stan_results_binatts_parsed[[1]],
                          stan_results_binwts_parsed[[1]],
                          stan_results_binattswts_parsed[[1]],
                          stan_results_oneatt_parsed[[1]],
                          stan_results_oneattbinatts_parsed[[1]])
  lmes_parsed[[i]] = list(stan_results_full_parsed[[2]],
                          stan_results_binatts_parsed[[2]],
                          stan_results_binwts_parsed[[2]],
                          stan_results_binattswts_parsed[[2]],
                          stan_results_oneatt_parsed[[2]],
                          stan_results_oneattbinatts_parsed[[2]])
  
  true_temps[[i]] = agent_inv_temps
  true_weights[[i]] = agent_weights
  
  # get best-fitting params from appropriate model
  best_fit_temps[[i]] = numeric(numSubj)
  best_fit_weights[[i]] = matrix(NA, nrow = numSubj, ncol = numAtts)
  for (subj in 1:numSubj) {
    fits_parsed_subj = fits_parsed[[i]][[i]][[subj]]
    if (!is.null(fits_parsed_subj)) {
      best_fit_params = fits_parsed_subj
      best_fit_temps[[i]][subj] = best_fit_params[1]
      best_fit_weights[[i]][subj,] = best_fit_params[2:(numAtts+1)] 
    } else {
      best_fit_temps[[i]][subj] = NA
    }
  }
}

### check model recoverability
numFittingModels = 6

bayes_matrix <- matrix(0, nrow = numGeneratingModels, ncol = numFittingModels)

numQuestions <- 3
q_indiv <- array(0, dim = c(numQuestions, numGeneratingModels, numSubj))
q_families <- list(
  list(1:4, 5:6),
  list(1:2, 3:4),
  list(c(1,3,5), c(2,4,6))
)

for (data_generating_model in 1:numGeneratingModels) {
  best_model <- numeric(numSubj)
  for (agent in 1:numSubj) {
    cur_lmes <- numeric(numFittingModels)
    for (fitted_model in 1:numFittingModels) {
      cur_lmes[fitted_model] <- lmes_parsed[[data_generating_model]][[fitted_model]][agent]
    }
    if (any(!is.na(cur_lmes))) {
      best_model[agent] <- which.max(cur_lmes)
    } else {
      best_model[agent] = NA
    }
    
    cur_mes <- exp(cur_lmes)
    
    for (question in 1:numQuestions) {
      q_indiv[question, data_generating_model, agent] <- mean(cur_mes[q_families[[question]][[1]]]) / (mean(cur_mes[q_families[[question]][[1]]]) + mean(cur_mes[q_families[[question]][[2]]]))
    }
  }
  
  for (fitted_model in 1:numFittingModels) {
    bayes_matrix[data_generating_model, fitted_model] <- sum(best_model == fitted_model, na.rm=T)
  }
}

q_avg <- apply(q_indiv, c(1, 2), mean, na.rm = T)
q_correct <- apply(q_indiv > 0.5, c(1, 2), mean, na.rm = T)
q_correct2 <- numeric(numQuestions)
for (question in 1:numQuestions) {
  q_correct2[question] <- sum(q_correct[question, q_families[[question]][[1]]]) / (sum(q_correct[question, q_families[[question]][[1]]]) + sum(q_correct[question, q_families[[question]][[2]]]))
}
q_correct2

bayes_matrix_normed1 <- matrix(0, nrow = numGeneratingModels, ncol = numFittingModels) # prob(fitted models | true models)
bayes_matrix_normed2 <- matrix(0, nrow = numGeneratingModels, ncol = numFittingModels) # prob(true models | fitted models)

for (i in 1:numGeneratingModels) {
  bayes_matrix_normed1[i,] <- bayes_matrix[i,] / sum(bayes_matrix[i,],na.rm=T)
}

for (i in 1:numFittingModels) {
  bayes_matrix_normed2[,i] <- bayes_matrix[,i] / sum(bayes_matrix[,i],na.rm=T)
}

mean(diag(bayes_matrix_normed1))
mean(diag(bayes_matrix_normed2))

### check parameter recoverability
att_cors = matrix(NA, nrow = numGeneratingModels, ncol = numAtts)
subj_cors = matrix(NA, nrow = numGeneratingModels, ncol = numSubj)
for (i in 1:numGeneratingModels) {
  for (att in 1:numAtts) {
    att_cors[i,att] = cor(true_weights[[i]][,att], best_fit_weights[[i]][,att], use = 'complete.obs')
  }
  
  for (subj in 1:numSubj) {
    if (all(!is.na(best_fit_weights[[i]][subj,]))) {
      subj_cors[i,subj] = cor(true_weights[[i]][subj,], best_fit_weights[[i]][subj,], use = 'complete.obs')
    }
  }
}

save.image('simulations/simulations_analysis_output.rdata')
