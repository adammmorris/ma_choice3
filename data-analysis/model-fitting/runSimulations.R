# setup -------------------------------------------------------------------
if (!require('pacman')) {
  install.packages('pacman')
  require('pacman')
}

p_load(this.path, parallel)

numAtts = 9
numTrials = 100

source(paste0(here(),'/model-fitting.R')) # get right filepath for this
load(paste0(here(),'/option_diffs.rdata')) # get right filepath for this

numSubj = 1000
numCores = 32

models_to_simulate = 1:6

stan_model_full <- stan_model(paste0(here(),"/stan-program.stan"))
stan_model_binwts <- stan_model(paste0(here(),"/stan-program-binwts.stan"))

# run simulations ---------------------------------------------------------

for (simulating_model in models_to_simulate) {
  stan_output_path = paste0(here(), '/simulations/', simulating_model, '/')
  fitting_results = vector(mode = 'list', numSubj)
    
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
  
  fitting_results = mclapply(
    1:numSubj,
    function(subj) {return(fitAllModels(option_diffs_touse[[subj]], choices_touse[[subj]]))},
    mc.cores = numCores
  )

  save(simulating_model, numSubj, numAtts,
       option_diffs_touse, choices_touse, agent_inv_temps, agent_weights,
       stan.seed, numChains, numIter,
       generateChoices, fitSubj,
       all_combinations_binwts, all_combinations_singleatt,
       fitting_results,
       file = paste0(stan_output_path, 'analysis_output_stan.rdata'))
}