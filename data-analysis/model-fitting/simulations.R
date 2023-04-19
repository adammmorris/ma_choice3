
# setup -------------------------------------------------------------------
require(this.path)

numAtts = 9
numTrials = 100

source(paste0(here(),'/model-fitting.R')) # get right filepath for this
load(paste0(here(),'/option_diffs.rdata')) # get right filepath for this

numSubj = 3

models_to_simulate = 1
numGeneratingModels = length(models_to_simulate)

stan_model_full <- stan_model(paste0(here(),"/stan-program.stan"))
stan_model_binwts <- stan_model(paste0(here(),"/stan-program-binwts.stan"))


# run simulations ---------------------------------------------------------

for (simulating_model in 1:numGeneratingModels) {
  stan_output_path = paste0('simulations/', simulating_model, '/')
  option_diffs_touse = sample(option_diffs_allsubj, numSubj, replace = T)
  agent_inv_temps = rgamma(numSubj, 4)
  
  if (simulating_model %in% c(1,2)) {
    agent_weights = runifm(numSubj, numAtts, min = -1, max = 1) 
  } else if (simulating_model %in% c(3,4)) {
    agent_weights = matrix(unlist(all_combinations_binwts[sample(num_comb_binwts, numSubj, replace = T)]), nrow = numSubj, byrow = T)
  } else {
    agent_weights = matrix(unlist(all_combinations_singleatt[sample(num_comb_singleatt, numSubj, replace = T)]), nrow = numSubj, byrow = T)
  }
  
  choices_touse = vector(mode = 'list', numSubj)
  for (subj in 1:numSubj) {
    option_diffs_touse_subj = option_diffs_touse[[subj]]
    if (simulating_model %in% c(2,4,6)) {option_diffs_touse_subj = sign(option_diffs_touse_subj)}
    choices_touse[[subj]] = generateChoices(agent_inv_temps[subj], agent_weights[subj,], option_diffs_touse_subj)
  }
  
  ## FIT FULL
  stan_results_full = doFitting(stan_model_full, option_diffs_touse, choices_touse, 0, 0)
  #stan_results_full_parsed = parseFittingResults(stan_results_full)
  
  ## FIT BINATT
  stan_results_binatts <- doFitting(stan_model_full, option_diffs_touse, choices_touse, 1, 0)
  #stan_results_binatts_parsed = parseFittingResults(stan_results_binatts)
  
  ## FIT BINWTS
  stan_results_binwts <- doFitting(stan_model_binwts, option_diffs_touse, choices_touse, 0, 1, all_combinations_binwts)
  #stan_results_binwts_parsed = parseFittingResults(stan_results_binwts, all_combinations_binwts)

  ## FIT BINATTSWTS
  stan_results_binattswts <- doFitting(stan_model_binwts, option_diffs_touse, choices_touse, 1, 1, all_combinations_binwts)
  #stan_results_binattswts_parsed = parseFittingResults(stan_results_binattswts, all_combinations_binwts)

  ## FIT 1ATT
  stan_results_oneatt <- doFitting(stan_model_binwts, option_diffs_touse, choices_touse, 0, 1, all_combinations_singleatt)
  #stan_results_oneatt_parsed = parseFittingResults(stan_results_oneatt, all_combinations_singleatt)
  
  ## FIT 1ATTBINATTS
  stan_results_oneattbinatts <- doFitting(stan_model_binwts, option_diffs_touse, choices_touse, 1, 1, all_combinations_singleatt)
  #stan_results_oneattbinatts_parsed = parseFittingResults(stan_results_oneattbinatts, all_combinations_singleatt)
  
  # save(simulating_model, numSubj, numAtts,
  #      option_diffs_touse, choices_touse, agent_inv_temps, agent_weights,
  #      stan.seed, numChains, numIter,
  #      generateChoices, doFitting, parseFittingResults,
  #      all_combinations_binwts, all_combinations_singleatt,
  #      stan_results_full, stan_results_full_parsed,
  #      stan_results_binatts, stan_results_binatts_parsed,
  #      stan_results_binwts, stan_results_binwts_parsed,
  #      stan_results_binattswts, stan_results_binattswts_parsed,
  #      stan_results_oneatt, stan_results_oneatt_parsed,
  #      stan_results_oneattbinatts, stan_results_oneattbinatts_parsed, stan_output_path,
  #      file = paste0(stan_output_path, 'analysis_output_stan.rdata'))

  rm(stan_results_full, stan_results_full_parsed,
       stan_results_binatts, stan_results_binatts_parsed,
       stan_results_binwts, stan_results_binwts_parsed,
       stan_results_binattswts, stan_results_binattswts_parsed,
       stan_results_oneatt, stan_results_oneatt_parsed,
       stan_results_oneattbinatts, stan_results_oneattbinatts_parsed)
  gc()
}

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
