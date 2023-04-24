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

p_load(rstan, bridgesampling, iterpc, matricks, loo)

options(mc.cores = 1)
rstan_options(auto_write = TRUE)

stan.seed = 12345
numChains = 4
numIter = 2000

numTopInds = 1000
numBadInds = 100

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

post_creator = function(u_num, u) {
  return(function(x) {
    dgamma(x, 4) * rowProds(exp(x %*% t(u_num)) / (1+exp(x %*% t(u))))
  })
}

fitSubj = function(stan_model_obj, option_diffs, choices, binary_atts, binary_wts, wts_combinations_list = NULL) {
  if (binary_atts) {option_diffs = sign(option_diffs)}
  
  # Fit the model
  if (length(choices) == numTrials) {
    if (!binary_wts) {
      stan_fit <- sampling(stan_model_obj,
                           data = list(numChoices = numTrials,
                                       numAtts = nrow(option_diffs),
                                       option_diffs = option_diffs,
                                       choices = choices),
                           iter = numIter,
                           chains = numChains,
                           seed = stan.seed,
                           refresh = 0)
      best_params = summary(stan_fit, probs = 0.5, pars = c('inv_temp', 'weights', 'lp__'))$summary[,1]
      lme = bridge_sampler(stan_fit, silent = T)$logml
      stan_samples = NULL#extract(stan_fit, pars = c('inv_temp', 'weights', 'lp__'))
      stan_diagnostics = get_all_diagnostics(stan_fit)
      stan_loo = loo(stan_fit)
      best_wts_ind = NA
      good_inds_to_test = NA
      lmes_goodinds = NA
      lps_badinds = NA
      lmes_badinds = NA
    } else {
      if (length(wts_combinations_list) > 18) {
        untemped_lpdfs = sapply(wts_combinations_list, function(x) compute_lpdf(1, x, option_diffs, choices))
        good_inds_to_test = which(untemped_lpdfs %in% tail(sort(untemped_lpdfs), numTopInds))
        bad_inds_to_test = sample(which(!(untemped_lpdfs %in% tail(sort(untemped_lpdfs), numTopInds))), numBadInds)
        do_sampling = T
      } else {
        good_inds_to_test = 1:length(wts_combinations_list)
        bad_inds_to_test = NA
        do_sampling = F
      }
      
      fits_goodinds = vector('list', length(good_inds_to_test))
      lmes_goodinds = rep(NA, length = length(good_inds_to_test))
      
      best_lp = -Inf
      best_wts_ind = NA
      for (i in 1:length(good_inds_to_test)) {
        weights = wts_combinations_list[[good_inds_to_test[i]]]
        utilities_untemped = as.vector(weights %*% option_diffs)
        utilities_num_untemped = utilities_untemped
        utilities_num_untemped[choices == 0] = 0
        f_post = post_creator(utilities_num_untemped, utilities_untemped)
        
        lmes_goodinds[i] = log(integrate(f_post, 0, 20)$value)
        fits_goodinds[[i]] = optimize(f_post, c(0, 20), maximum = T)
        
        cur_lp = log(fits_goodinds[[i]]$objective)
        if (cur_lp > best_lp) {
          best_wts_ind = i
          best_lp = cur_lp
        }
      }
      
      # if we're in binwts, get random sample of bad ones to approximate default lme
      lps_badinds = rep(NA, length = length(bad_inds_to_test))
      lmes_badinds = rep(NA, length = length(bad_inds_to_test))
      
      if (do_sampling) {
        for (i in 1:length(bad_inds_to_test)) {
          weights = wts_combinations_list[[bad_inds_to_test[i]]]
          utilities_untemped = as.vector(weights %*% option_diffs)
          utilities_num_untemped = utilities_untemped
          utilities_num_untemped[choices == 0] = 0
          f_post = post_creator(utilities_num_untemped, utilities_untemped)

          lmes_badinds[i] = log(integrate(f_post, 0, 20)$value)
          lps_badinds[i] = log(optimize(f_post, c(0, 20), maximum = T)$objective)
        }
        
        lme = log(mean(exp(
          c(lmes_goodinds,
            sample(lmes_badinds, length(wts_combinations_list) - length(lmes_goodinds), replace = T)
          )
        )))
      } else {
        lme = log(mean(exp(lmes_goodinds)))
      }
      
      # for best one, save more stuff
      best_fit = fits_goodinds[[best_wts_ind]]
      best_wts = wts_combinations_list[[good_inds_to_test[best_wts_ind]]]
      best_params = c(best_fit$maximum, best_wts, log(best_fit$objective))
      
      utilities_untemped = as.vector(best_wts %*% option_diffs)
      utilities_num_untemped = utilities_untemped
      utilities_num_untemped[choices == 0] = 0
      stan_fit <- sampling(stan_model_obj,
                                     data = list(numChoices = numTrials,
                                                 utilities_untemped = utilities_untemped,
                                                 utilities_num_untemped = utilities_num_untemped),
                                     iter = numIter,
                                     chains = numChains,
                                     seed = stan.seed,
                                     refresh = 0)
      stan_samples = extract(stan_fit, pars = c('inv_temp', 'lp__'))
      stan_diagnostics = get_all_diagnostics(stan_fit)
      stan_loo = loo(stan_fit)
    }
    
    # for space... too big to store
    #stan_fit@sim$samples = NULL
    
    stan_results = list(best_params, lme)#, stan_fit, stan_samples, stan_diagnostics, stan_loo, best_wts_ind, good_inds_to_test, lmes_goodinds, lmes_badinds, lps_badinds)
  } else {
    stan_results = list(NULL, NULL)#, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
  }
  
  return(stan_results)
}

fitAllModels = function(option_diffs, choices) {
  fitting_results = vector('list', 6)
  
  fitting_results[[1]] = fitSubj(stan_model_full, option_diffs, choices, 0, 0)
  fitting_results[[2]] <- fitSubj(stan_model_full, option_diffs, choices, 1, 0)
  fitting_results[[3]] <- fitSubj(stan_model_binwts, option_diffs, choices, 0, 1, all_combinations_binwts)
  fitting_results[[4]] <- fitSubj(stan_model_binwts, option_diffs, choices, 1, 1, all_combinations_binwts)
  fitting_results[[5]] <- fitSubj(stan_model_binwts, option_diffs, choices, 0, 1, all_combinations_singleatt)
  fitting_results[[6]] <- fitSubj(stan_model_binwts, option_diffs, choices, 1, 1, all_combinations_singleatt)
  gc()
  return(fitting_results)
}