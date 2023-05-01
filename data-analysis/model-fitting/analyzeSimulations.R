# analyze results ---------------------------------------------------------
require(this.path)
setwd(here())

numGeneratingModels = 6
numFittingModels = 6

### import data
lmes = vector('list', numGeneratingModels)
true_temps = vector('list', numGeneratingModels)
true_weights = vector('list', numGeneratingModels)
best_fit_temps = vector('list', numGeneratingModels)
best_fit_weights = vector('list', numGeneratingModels)

for (generating_model in 1:numGeneratingModels) {
  load(paste0('simulations/', generating_model, '/analysis_output_stan.rdata'))
  
  # get true parameters
  true_temps[[generating_model]] = agent_inv_temps
  true_weights[[generating_model]] = agent_weights
  
  # get best-fitting params
  best_fit_temps[[generating_model]] = rep(NA, numSubj)
  best_fit_weights[[generating_model]] = matrix(NA, nrow = numSubj, ncol = numAtts)
  lmes[[generating_model]] = matrix(NA, nrow = numSubj, ncol = numGeneratingModels)
  for (subj in 1:numSubj) {
    fit_subj = fitting_results[[subj]][[generating_model]] # from corresponding fitted model
    if (!is.null(fit_subj[[1]])) {
      best_fit_params = fit_subj[[1]]
      best_fit_temps[[generating_model]][subj] = best_fit_params[1]
      best_fit_weights[[generating_model]][subj,] = best_fit_params[2:(numAtts+1)] 
      
      for (fitted_model in 1:numFittingModels) {
        lmes[[generating_model]][subj,fitted_model] = fitting_results[[subj]][[fitted_model]][[2]]
      }
    }
  }
}

### check model recoverability
bayes_matrix <- matrix(0, nrow = numGeneratingModels, ncol = numFittingModels)

numQuestions <- 3
q_indiv <- array(0, dim = c(numQuestions, numGeneratingModels, numSubj))
q_families <- list(
  list(1:4, 5:6),
  list(1:2, 3:4),
  list(c(1,3,5), c(2,4,6))
)

for (generating_model in 1:numGeneratingModels) {
  best_model <- numeric(numSubj)
  for (agent in 1:numSubj) {
    cur_lmes <- lmes[[generating_model]][agent,]
    if (any(!is.na(cur_lmes))) {
      best_model[agent] <- which.max(cur_lmes)
    } else {
      best_model[agent] = NA
    }
    
    cur_mes <- exp(cur_lmes)
    
    for (question in 1:numQuestions) {
      q_indiv[question, generating_model, agent] <- mean(cur_mes[q_families[[question]][[1]]]) / (mean(cur_mes[q_families[[question]][[1]]]) + mean(cur_mes[q_families[[question]][[2]]]))
    }
  }
  
  for (fitted_model in 1:numFittingModels) {
    bayes_matrix[generating_model, fitted_model] <- sum(best_model == fitted_model, na.rm=T)
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
