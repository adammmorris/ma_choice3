# setup -------------------------------------------------------------------
if (!require('pacman')) {
  install.packages('pacman')
  require('pacman')
}

p_load(this.path)

numAtts = 9
numTrials = 100

source(paste0(here(),'/model-fitting.R')) # get right filepath for this
load(paste0(here(),'/option_diffs.rdata')) # get right filepath for this

numSubj = 32

models_to_simulate = 1
numGeneratingModels = length(models_to_simulate)

stan_model_full <- stan_model(paste0(here(),"/stan-program.stan"))
stan_model_binwts <- stan_model(paste0(here(),"/stan-program-binwts.stan"))


# run simulations ---------------------------------------------------------

for (simulating_model in 1:numGeneratingModels) {
  stan_output_path = paste0(here(), '/simulations/', simulating_model, '/')
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
  
  stan_results_full = doFitting(stan_model_full, option_diffs_touse, choices_touse, 0, 0)
  stan_results_binatts <- doFitting(stan_model_full, option_diffs_touse, choices_touse, 1, 0)
  stan_results_binwts <- doFitting(stan_model_binwts, option_diffs_touse, choices_touse, 0, 1, all_combinations_binwts)
  stan_results_binattswts <- doFitting(stan_model_binwts, option_diffs_touse, choices_touse, 1, 1, all_combinations_binwts)
  stan_results_oneatt <- doFitting(stan_model_binwts, option_diffs_touse, choices_touse, 0, 1, all_combinations_singleatt)
  stan_results_oneattbinatts <- doFitting(stan_model_binwts, option_diffs_touse, choices_touse, 1, 1, all_combinations_singleatt)

  save(simulating_model, numSubj, numAtts,
       option_diffs_touse, choices_touse, agent_inv_temps, agent_weights,
       stan.seed, numChains, numIter,
       generateChoices, doFitting, parseFittingResults,
       all_combinations_binwts, all_combinations_singleatt,
       stan_results_full, stan_results_binatts, stan_results_binwts,
       stan_results_binattswts, stan_results_oneatt, stan_results_oneattbinatts,
       stan_output_path,
       file = paste0(stan_output_path, 'analysis_output_stan.rdata'))

  rm(stan_results_full, stan_results_binatts, stan_results_binwts,
     stan_results_binattswts, stan_results_oneatt, stan_results_oneattbinatts)
  gc()
}