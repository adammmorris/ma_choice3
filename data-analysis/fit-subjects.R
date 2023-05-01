# Setup -------------------------------------------------------------

if (!require('pacman')) {
  install.packages('pacman')
  require('pacman')
}

p_load(this.path, parallel, matrixStats)


# Set version -------------------------------------------------------------

## which version do we want?
# note that we don't fit observers
versions = c('home-fixed', 'home-random', 'movie')
version = versions[3]

filepath = paste0(here(), '/', version, '/')

load(paste0(filepath,'analysis_output.rdata'))

numSubj = length(subjlist)
numAtts = length(atts)
source(paste0(here(),'/model-fitting/model-fitting.R'))
stan_model_full <- stan_model(paste0(here(),"/model-fitting/stan-program.stan"))
stan_model_binwts <- stan_model(paste0(here(),"/model-fitting/stan-program-binwts.stan"))

# Do model-fitting --------------------------------------------------------
numCores = 32

option_diffs_allsubj = vector(mode = 'list', numSubj)
for (subj in 1:numSubj) {
  option_diffs_allsubj[[subj]] = t(as.matrix(df.s1 %>%
                                               filter(subject == subjlist[subj]) %>%
                                               select(all_of(atts.opt.scaled.diff))))
}
option_diffs_touse = option_diffs_allsubj
save(option_diffs_allsubj, file = paste0(here(), '/model-fitting/option_diffs.rdata')) # save for simulations

choices_touse = vector(mode = 'list', numSubj)
for (subj in 1:numSubj) {
  choices_touse[[subj]] = df.s1[df.s1$subject == subjlist[subj],'choice']
}

## FIT FULL
fitting_results = mclapply(
  1:numSubj,
  function(subj) {return(fitAllModels(option_diffs_touse[[subj]], choices_touse[[subj]], full_output = TRUE))},
  mc.cores = numCores
)

# Split so that each saved data file is <100mb
fitting_results_split = split(fitting_results, ceiling(seq_along(fitting_results) / 20))

for (i in 1:length(fitting_results_split)) {
  fitting_results_split_cur = fitting_results_split[[i]]
  save(numAtts, option_diffs_touse, choices_touse, filepath,
       stan.seed, numChains, numIter,
       generateChoices, fitSubj,
       fitting_results_split_cur,
       file = paste0(filepath, 'modeling-output/modeling_output_', i, '.rdata'))
}
