# setup -------------------------------------------------------------------
if (!require('pacman')) {
  install.packages('pacman')
  require('pacman')
}

p_load(this.path, parallel, matrixStats)

numAtts = 9
numTrials = 100

source(paste0(here(),'/model-fitting.R')) # get right filepath for this
load(paste0(here(),'/option_diffs.rdata')) # get right filepath for this

numSubj = 100
numCores = 32

models_to_simulate = 1

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
  
  #debugging
  # test <- fitSubj(stan_model_binwts, option_diffs_touse[[subj]], choices_touse[[subj]], 0, 1, all_combinations_binwts)
  # weights = agent_weights[subj,]
  # utilities_untemped = as.vector(weights %*% option_diffs_touse[[subj]])
  # utilities_num_untemped = utilities_untemped
  # utilities_num_untemped[choices_touse[[subj]] == 0] = 0
  # 
  # begin = Sys.time()
  # stan_fit <- sampling(stan_model_binwts,
  #                                data = list(numChoices = 100,
  #                                            utilities_untemped = utilities_untemped,
  #                                            utilities_num_untemped = utilities_num_untemped),
  #                                iter = numIter,
  #                                chains = numChains,
  #                                seed = stan.seed,
  #                                refresh = 0)
  # stan_lme = bridge_sampler(stan_fit, silent = T)$logml
  # Sys.time() - begin
  # 
  # 
  # f_creator = function(u1, u2) {return(function(x) {dgamma(x, 4) * rowProds(exp(x %*% t(u1)) / (1+exp(x %*% t(u2))))})}
  # 
  # begin = Sys.time()
  # f = f_creator(utilities_num_untemped, utilities_untemped)
  # log(integrate(f, 0, 20)$value)
  # test = optimize(f, c(0,20), maximum = T)
  # Sys.time() - begin
  
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