# needs: numAtts, numTrials

# setup -------------------------------------------------------------------
#require(groundhog)

#pkg.names.modelfitting = c('rstan', 'bridgesampling', 'iterpc', 'foreach', 'doParallel', 'matricks','loo')
#groundhog.library(pkg.names, '2023-03-11')
#lapply(pkg.names.modelfitting, require, character.only = TRUE)

if (!require('pacman')) {
  install.packages('pacman')
  require('pacman')
}

p_load(rstan, bridgesampling, iterpc, foreach, doParallel, matricks, loo)

options(mc.cores = 1)
rstan_options(auto_write = TRUE)
registerDoParallel(cores = 32)

stan.seed = 12345
numChains = 2
numIter = 2000

numReturnsPerFit = 8
numTopInds = 100

# get all combinations, for binwts & single att
elements <- c(-1L, 0L, 1L)
n <- 9
I <- iterpc(length(elements), n, ordered = TRUE, replace = TRUE)
all_combinations_matrix <- getall(I)
all_combinations_binwts <- split(all_combinations_matrix, row(all_combinations_matrix))
all_combinations_binwts <- lapply(all_combinations_binwts, function(x) elements[x])
num_comb_binwts = length(all_combinations_binwts)

all_combinations_singleatt = NULL
ind = 1
for (i in 1:numAtts) {
  all_combinations_singleatt[[ind]] = rep(0,numAtts)
  all_combinations_singleatt[[ind]][i] = -1
  all_combinations_singleatt[[ind+1]] = rep(0,numAtts)
  all_combinations_singleatt[[ind+1]][i] = 1
  ind = ind + 2
}
num_comb_singleatt = length(all_combinations_singleatt)

compute_lpdf = function(inv_temp, weights, option_diffs, choices) {
  utilities = inv_temp * as.vector(weights %*% option_diffs)
  utilities_num = utilities
  utilities_num[choices == 0] = 0
  return(sum(utilities_num - log(exp(utilities)+1)))
}

# critical functions ------------------------------------------------------

generateChoices = function(inv_temp, weights, option_diffs) {
  num_choices <- ncol(option_diffs)
  utilities <- inv_temp * (weights %*% option_diffs)
  choice_probs <- exp(utilities) / (exp(utilities) + 1)
  
  choices = as.numeric(runif(num_choices) < choice_probs)
  return(choices)
}

# adapted from Michael Betancort, https://betanalpha.github.io/assets/case_studies/rstan_workflow.html
get_all_diagnostics = function(fit) {
  output = vector('character', 5)
  sampler_params <- get_sampler_params(fit, inc_warmup=FALSE)
  
  # divergence
  divergent <- do.call(rbind, sampler_params)[,'divergent__']
  n = sum(divergent)
  N = length(divergent)
  
  output[1] = as.character(n)#sprintf('%s of %s iterations ended with a divergence (%s%%)',
                      #n, N, 100 * n / N)
  
  # tree depth
  max_depth = 10
  treedepths <- do.call(rbind, sampler_params)[,'treedepth__']
  n = length(treedepths[sapply(treedepths, function(x) x == max_depth)])
  N = length(treedepths)
  
  output[2] = sprintf('%s of %s iterations saturated the maximum tree depth of %s (%s%%)',
                      n, N, max_depth, 100 * n / N)
  
  # E-BFMI
  no_warning <- TRUE
  for (n in 1:length(sampler_params)) {
    energies = sampler_params[n][[1]][,'energy__']
    numer = sum(diff(energies)**2) / length(energies)
    denom = var(energies)
    if (numer / denom < 0.2) {
      print(sprintf('Chain %s: E-BFMI = %s', n, numer / denom))
      no_warning = FALSE
    }
  }
  
  if (no_warning)
    output[3] = 'E-BFMI indicated no pathological behavior'
  else
    output[3] = 'E-BFMI below 0.2 indicates you may need to reparameterize your model'
  
  # neff and rhat
  fit_summary <- summary(fit, probs = c(0.5))$summary
  N <- dim(fit_summary)[[1]]
  
  iter <- dim(extract(fit)[[1]])[[1]]
  
  for (n in 1:N) {
    ratio <- fit_summary[,5][n] / iter
    if (ratio < 0.001) {
      output[4] = paste0(output[4], sprintf('n_eff / iter for parameter %s is %s!',
                                            rownames(fit_summary)[n], ratio))
    }
    
    rhat <- fit_summary[,6][n]
    if (rhat > 1.1 || is.infinite(rhat) || is.nan(rhat)) {
      output[5] = paste0(output[5], sprintf('Rhat for parameter %s is %s!',
                                            rownames(fit_summary)[n], rhat))
    }
  }
  
  if (output[4] == '') output[4] = 'n_eff is good'
  if (output[5] == '') output[5] = 'Rhat is good'
  
  return(output)
}

## fitting function
doFitting = function(stan_model_obj, option_diffs, choices, binary_atts, binary_wts, wts_combinations_list = NULL) {
  stan_results <- foreach(subj = 1:length(option_diffs), .combine = 'c', .multicombine = T, .packages = c("rstan", "bridgesampling", "dplyr")) %dopar% {
    option_diffs_subj = option_diffs[[subj]]
    if (binary_atts) {option_diffs_subj = sign(option_diffs_subj)}
    choices_subj = choices[[subj]]
    
    # Fit the model
    if (length(choices_subj) == numTrials) {
      if (!binary_wts) {
        stan_fit <- sampling(stan_model_obj,
                             data = list(numChoices = numTrials,
                                         numAtts = nrow(option_diffs_subj),
                                         option_diffs = option_diffs_subj,
                                         choices = choices_subj),
                             iter = numIter,
                             chains = numChains,
                             seed = stan.seed,
                             refresh = 0)
        stan_params = summary(stan_fit, probs = 0.5, pars = c('inv_temp', 'weights', 'lp__'))$summary[,1]
        stan_lme = bridge_sampler(stan_fit, silent = T)$logml
        stan_diagnostics = get_all_diagnostics(stan_fit)
        stan_loo = loo(stan_fit)
        best_wts_ind = NA
        good_inds_to_test = NA
        stan_lmes_all = NA
        stan_lp_bad = NA
        stan_lme_bad = NA
      } else {
        if (length(wts_combinations_list) > 18) {
          untemped_lpdfs = sapply(wts_combinations_list, function(x) compute_lpdf(1, x, option_diffs_subj, choices_subj))
          good_inds_to_test = which(untemped_lpdfs %in% tail(sort(untemped_lpdfs), numTopInds))
          bad_inds_to_test = sample(which(!(untemped_lpdfs %in% tail(sort(untemped_lpdfs), numTopInds))), 10)
        } else {
          good_inds_to_test = 1:length(wts_combinations_list)
          bad_inds_to_test = NA
        }
        
        stan_fits_all = vector('list', length(wts_combinations_list))
        stan_params = vector('list', length(wts_combinations_list))
        stan_lmes_all = rep(NA, length = length(wts_combinations_list))
        
        best_lp = -Inf
        best_wts_ind = NA
        for (i in good_inds_to_test) {
          weights = wts_combinations_list[[i]]
          utilities_untemped = as.vector(weights %*% option_diffs_subj)
          utilities_num_untemped = utilities_untemped
          utilities_num_untemped[choices_subj == 0] = 0
          stan_fits_all[[i]] <- sampling(stan_model_obj,
                                    data = list(numChoices = numTrials,
                                                utilities_untemped = utilities_untemped,
                                                utilities_num_untemped = utilities_num_untemped),
                                    iter = numIter,
                                    chains = numChains,
                                    seed = stan.seed,
                                    refresh = 0)
          stan_lmes_all[i] = bridge_sampler(stan_fits_all[[i]], silent = T)$logml
          
          cur_lp = summary(stan_fits_all[[i]], probs = 0.5, pars=c('lp__'))$summary[1]
          if (cur_lp > best_lp) {
            best_wts_ind = i
            best_lp = cur_lp
          }
        }
        
        # if we're in binwts, get random sample of bad ones to approximate default lme
        if (any(!is.na(bad_inds_to_test))) {
          stan_lps_bad_all = rep(NA, length = length(bad_inds_to_test))
          stan_lmes_bad_all = rep(NA, length = length(bad_inds_to_test))
          for (i in 1:length(bad_inds_to_test)) {
            weights = wts_combinations_list[[bad_inds_to_test[i]]]
            utilities_untemped = as.vector(weights %*% option_diffs_subj)
            utilities_num_untemped = utilities_untemped
            utilities_num_untemped[choices_subj == 0] = 0
            stan_fit_bad_cur <- sampling(stan_model_obj,
                                           data = list(numChoices = numTrials,
                                                       utilities_untemped = utilities_untemped,
                                                       utilities_num_untemped = utilities_num_untemped),
                                           iter = numIter,
                                           chains = numChains,
                                           seed = stan.seed,
                                           refresh = 0)
            stan_lps_bad_all[i] = summary(stan_fit_bad_cur, probs = 0.5, pars=c('lp__'))$summary[1]
            stan_lmes_bad_all[i] = bridge_sampler(stan_fit_bad_cur, silent = T)$logml
          }
          stan_lp_bad = mean(stan_lps_bad_all)
          stan_lme_bad = mean(stan_lmes_bad_all)
          stan_lmes_all[is.na(stan_lmes_all)] = stan_lme_bad
        } else {
          stan_lp_bad = NA
          stan_lme_bad = NA
        }
        
        # for best one, save more stuff
        stan_fit = stan_fits_all[[best_wts_ind]]
        stan_lme = log(mean(exp(stan_lmes_all)))
        stan_params_temp = summary(stan_fit, probs = 0.5, pars=c('inv_temp', 'lp__'))$summary[,1]
        stan_params = c(stan_params_temp['inv_temp'], wts_combinations_list[[best_wts_ind]], stan_params_temp['lp__'])
        stan_diagnostics = get_all_diagnostics(stan_fit)
        stan_loo = loo(stan_fit)
        rm(stan_fits_all)
      }
      
      list(stan_fit, stan_params, stan_lme, stan_diagnostics, stan_loo, best_wts_ind, good_inds_to_test, stan_lmes_all[good_inds_to_test], stan_lp_bad, stan_lme_bad)
    } else {
      list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
    }
  }
  
  return(stan_results)
}

# returns list w/: best-fitting params, LMEs
parseFittingResults = function(stan_results, wts_combinations_list = NULL) {
  numSubj = length(stan_results) / numReturnsPerFit
  stan_fits = vector('list', numSubj)
  stan_lmes = rep(NA, length = numSubj)
  stan_diagnostics = vector('list', numSubj)
  
  ind = 1
  for (i in 1:numSubj) {
    if (!is.null(stan_results[[ind]])) {
      if (is.null(wts_combinations_list)) {
        stan_fits[[i]] = get_posterior_mean(stan_results[[ind]], pars=c('inv_temp', 'weights', 'lp__'))[, ifelse(numChains == 1, 1, numChains+1)]
        stan_lmes[i] = stan_results[[ind+1]]$logml
        stan_diagnostics[[i]] = stan_results[[ind+2]]
      } else {
        stan_fitresults_all = stan_results[[ind]]
        stan_lmeresults_all = stan_results[[ind+1]]
        stan_diagnostics_all = stan_results[[ind+2]]
        
        best_lp = -Inf
        best_wts_ind = NA
        stan_lmes_all = rep(log(.5^100), length = length(wts_combinations_list))
        for (j in 1:length(wts_combinations_list)) {
          if (!is.null(stan_fitresults_all[[j]])) {
            cur_lp = stan_fitresults_all[[j]]['lp__']
            if (cur_lp > best_lp) {
              best_wts_ind = j
              best_lp = cur_lp
            }
            
            stan_lmes_all[j] = stan_lmeresults_all[[j]]$logml
          }
        }
        temp_var = stan_fitresults_all[[best_wts_ind]]
        stan_fits[[i]] = c(temp_var['inv_temp'], wts_combinations_list[[best_wts_ind]], temp_var['lp__'])
        stan_lmes[i] = log(mean(exp(stan_lmes_all)))
        stan_diagnostics[[i]] = stan_diagnostics_all[[best_wts_ind]]
      }
    }
    ind = ind + numReturnsPerFit
  } 
  
  return(list(stan_fits, stan_lmes, stan_diagnostics))
}
