### This code has been tested in R version 4.2.2.

# Setup -------------------------------------------------------------------
# if (!require('pacman')) {
#   install.packages('pacman')
#   require('pacman')
# }

require(groundhog)

pkg.names = c('ggplot2', 'lme4', 'lmerTest', 'tidyverse', 'jsonlite', 'combinat', 'effectsize', 'RColorBrewer', 'scales')
#groundhog.library(pkg.names, '2023-03-11')
lapply(pkg.names, require, character.only = TRUE)

theme_update(strip.background = element_blank(),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             panel.background = element_blank(),
             plot.background = element_blank(),
             axis.text=element_text(size=25, colour = "black"),
             axis.title=element_text(size=30, face = "bold"),
             axis.title.x = element_text(vjust = 0),
             legend.title = element_text(size = 24, face = "bold"),
             legend.text = element_text(size = 20),
             plot.title = element_text(size = 26, face = "bold", vjust = 1),
             panel.margin = unit(1.0, "lines"), 
             plot.margin = unit(c(0.5,  0.5, 0.5, 0.5), "lines"),
             axis.line = element_line(colour = "black", size = 2),
             axis.ticks = element_line(color = 'black', size = 3),
             axis.ticks.length = unit(.25, 'cm')
)

theme_black = function(base_size = 12, base_family = "") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # Specify axis options
      axis.line = element_blank(),  
      axis.text.x = element_text(size = 12, color = "white", lineheight = 0.9),  
      axis.text.y = element_text(size = 12, color = "white", lineheight = 0.9),  
      axis.ticks = element_line(color = "white", size  =  0.2),  
      axis.title.x = element_text(size = 18, color = "white", margin = margin(0, 10, 0, 0)),  
      axis.title.y = element_text(size = 18, color = "white", angle = 90, margin = margin(0, 10, 0, 0)),  
      axis.ticks.length = unit(0.3, "lines"),   
      # Specify legend options
      legend.background = element_rect(color = NA, fill = "black"),  
      legend.key = element_rect(color = "white",  fill = "black"),  
      legend.key.size = unit(1.2, "lines"),  
      legend.key.height = NULL,  
      legend.key.width = NULL,      
      legend.text = element_text(size = base_size*0.8, color = "white"),  
      legend.title = element_text(size = base_size*0.8, face = "bold", hjust = 0, color = "white"),  
      legend.position = "right",  
      legend.text.align = NULL,  
      legend.title.align = NULL,  
      legend.direction = "vertical",  
      legend.box = NULL, 
      # Specify panel options
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA, color = "white"),  
      # Specify facetting options
      strip.background = element_rect(fill = "grey30", color = "grey10"),  
      strip.text.x = element_text(size = base_size*0.8, color = "white"),  
      strip.text.y = element_text(size = base_size*0.8, color = "white",angle = -90),  
      # Specify plot options
      plot.background = element_rect(color = "black", fill = "black"),  
      plot.title = element_text(size = base_size*1.2, color = "white"),  
      plot.margin = unit(rep(1, 4), "lines")
    )
}
se = function(x) {return(sd(x, na.rm = T) / sqrt(sum(!is.na(x))))}
se.prop = function(x) {return(sqrt(mean(x, na.rm = T) * (1-mean(x, na.rm = T)) / sum(!is.na(x))))}
get.ci = function(x) {return(c(mean(x,na.rm = T) - 1.96*se(x), mean(x, na.rm = T), mean(x, na.rm = T) + 1.96*se(x)))}
get.ci.prop = function(x) {return(c(mean(x,na.rm = T) - 1.96*se.prop(x), mean(x, na.rm = T), mean(x, na.rm = T) + 1.96*se.prop(x)))}

as.string.vector = function(x) {
  return(strsplit(x,',')[[1]])
}
as.numeric.vector = function(x) {
  return(as.numeric(strsplit(gsub('\\[|\\]','',x),',')[[1]]))
}
as.string = function(x) {
  return(paste(x, collapse = ','))
}
dodge <- position_dodge(width=0.9)

# for loading model fitting results
combine_lists <- function(directory) {
  # Get a list of .RData files in the directory
  files <- list.files(directory, pattern="\\.rdata$", full.names=TRUE)
  
  # Create an empty list to store the data
  all_data <- list()
  
  # Load the data from each file
  for (i in 1:length(files)) {
    load(paste0(directory, 'modeling_output_', i, '.rdata'))
    all_data <- c(all_data, fitting_results_split_cur)
  }
  
  return(all_data)
}


# Only works in RStudio -- otherwise you have to set the path manually
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Set version -------------------------------------------------------------

## which version do we want?
# note: the observer versions are just in here for completeness, but you (as the reader)
# should never have to run that code -- the output of the observer analysis gets fed into
# the decider analysis (in the "Compare deciders to observers" section below), where we
# do the critical comparisons between versions.
versions = c('home-fixed', 'home-random', 'movie', 'home-fixed-observers', 'movie-observers')
version = versions[3]
task_variant = ifelse(version %in% c('home-fixed', 'home-random', 'home-fixed-observers'),
                      'home',
                      'movie')
participant_type = ifelse(version %in% c('home-fixed', 'home-random', 'movie'),
                          'decider',
                          'observer')
has_observers = version %in% c('home-fixed', 'movie')
filepath = paste0(version, '/')
filepath_data = paste0(filepath, 'data/')
if (participant_type == 'observer') {
  filepath_decider = paste0(sub("-observers", "", version), '/');
} else {
  filepath_observer = paste0(version, '-observers/')
}

# Load data ---------------------------------------------------------------
# h1.levels =c('One', 'Multiple')
# h2.levels =c('Binary', 'Graded')
# h3.levels =c('Binary', 'Graded')
numTrials = 100

# df.demo (short for "demographics") has all the subject-level information, i.e. one row per subject
df.demo = read.csv(paste0(filepath_data, 'demo.csv'), stringsAsFactors = F) %>% arrange(subject) %>%
  rowwise() %>%
  mutate(total_time_real = total_time / 60000,
         instructions_times_list = list(as.numeric.vector(instructions_times) / 1000),
         stamp = parse_date_time(stamp, c('ymd HMS', 'mdy HMS', 'mdy HM'))) %>%
  ungroup() 
for (i in 1:nrow(df.demo)) {
  df.demo$instruction_times_median[i] = median(df.demo$instructions_times_list[i][[1]])
  df.demo$instruction_times_sd[i] = sd(df.demo$instructions_times_list[i][[1]])
}

df.demo.followup = read.csv(paste0(filepath_data, 'demo_followup.csv'), stringsAsFactors = F) %>% arrange(subject) %>%
  rowwise() %>%
  mutate(total_time_real = total_time / 60000,
         instructions_times_list = list(as.numeric.vector(instructions_times) / 1000))

# df.s1 has the data from stage 1 (where people make choices)
df.s1 = read.csv(paste0(filepath_data,'s1.csv'), stringsAsFactors = F)
# df.s2 has the data from stage 2 (where people report how they made their choices)
df.s2 = read.csv(paste0(filepath_data,'s2.csv'), stringsAsFactors = F)
# df.attributes has the info about the set of choice attributes presented to that subject
df.attributes = read.csv(paste0(filepath_data,'attributes.csv'), stringsAsFactors = F)

# filter out anyone who didn't finish
subjlist.demo = unique(df.demo$subject)
subjlist.s1 = unique(df.s1$subject)
subjlist.s2 = unique(df.s2$subject)
subjlist.attributes = unique(df.attributes$subject)
subjlist = Reduce(intersect, list(subjlist.demo, subjlist.s1, subjlist.s2, subjlist.attributes))

df.demo = df.demo %>% filter(subject %in% subjlist) %>%
  mutate(subject = factor(subject), subject.num = as.numeric(subject))
df.s1 = df.s1 %>% filter(subject %in% subjlist, practice == 0) %>%
  arrange(subject) %>% mutate(subject = factor(subject), subject.num = as.numeric(subject))
df.s1.practice = df.s1 %>% filter(subject %in% subjlist, practice == 1) %>%
  arrange(subject) %>% mutate(subject = factor(subject), subject.num = as.numeric(subject))
df.s2 = df.s2 %>% filter(subject %in% subjlist) %>%
  arrange(subject)
# there was a bug in this version where the columns got reversed
if (version %in% c('home-fixed', 'home-fixed-observers')) {
  df.s2 = df.s2 %>% mutate(subject = factor(subject), subject.num = as.numeric(subject),
         least_preferred_temp = least_preferred, least_preferred = most_preferred,
         most_preferred = least_preferred_temp) %>%
  dplyr::select(-least_preferred_temp)
}
df.attributes = df.attributes %>% filter(subject %in% subjlist) %>%
  arrange(subject) %>% mutate(subject = factor(subject), subject.num = as.numeric(subject))

# Get & clean ANT data ---------------------------------------------------------------------

# some versions don't have the ANT
if (file.exists(paste0(filepath_data,'ant.csv'))) {
  df.ant.raw = read.csv(paste0(filepath_data,'ant.csv'), stringsAsFactors = F) %>% arrange(subject) %>%
    mutate(flanker_type = factor(flanker_type, c('neutral', 'congruent', 'incongruent')),
           cue = factor(cue, c('nocue', 'center', 'double', 'spatial')))
  df.ant = df.ant.raw %>% filter(subject %in% subjlist, practice == 0)
  
  mean(df.ant$correct)
  mean(df.ant$timed_out)
  hist(df.ant$rt)
  
  df.ant.subj = df.ant %>% group_by(subject) %>%
    summarize(correct.m = mean(correct), correct.se = se.prop(correct),
              timed_out.m = mean(timed_out), timed_out.se = se.prop(timed_out),
              rt.m = mean(rt, na.rm = T), rt.se = se(rt),
              num.trials = n())
  hist(df.ant.subj$correct.m)
  hist(df.ant.subj$rt.m)
  hist(df.ant.subj$timed_out.m)
  
  # filter out bad subj
  exclude.subj.ant = (df.ant.subj %>% filter(correct.m < .7 | rt.m < 350 | timed_out.m > .1 | df.ant.subj$num.trials != 288))$subject
  
  df.ant.filt = df.ant %>% filter(!(subject %in% exclude.subj.ant), rt < 1200, rt > 200, !timed_out)
  
  hist(df.ant.filt$rt)
  hist(log(df.ant.filt$rt))
  
  # graph
  df.ant.filt.graph = df.ant.filt %>% group_by(flanker_type, cue, subject) %>%
    summarize(correct = mean(correct),
              rt = mean(rt)) %>%
    group_by(flanker_type, cue) %>%
    summarize(correct.m = mean(correct), correct.se = se(correct),
              rt.m = mean(rt), rt.se = se(rt))
  
  df.ant.filt.graph.fake = df.ant.filt %>%
    group_by(flanker_type, cue) %>%
    summarize(correct.m = mean(correct), correct.se = se.prop(correct),
              rt.m = mean(rt), rt.se = se(rt))
  
  ggplot(df.ant.filt.graph, aes(x = flanker_type, group = cue, color = cue, y = rt.m)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = rt.m - rt.se, ymax = rt.m + rt.se))
  
  ggplot(df.ant.filt.graph, aes(x = flanker_type, group = cue, color = cue, y = correct.m)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = correct.m - correct.se, ymax = correct.m + correct.se))
  
  # alerting effect
  df.ant.filt.alert = df.ant.filt %>% filter(cue %in% c('nocue', 'double')) %>%
    group_by(cue, subject) %>%
    summarize(correct = mean(correct),
              rt = mean(rt)) %>%
    arrange(subject) %>%
    group_by(subject) %>%
    mutate(alerting = rt - lead(rt)) %>%
    filter(cue == 'nocue') %>%
    dplyr::select(subject, alerting)
  
  # orienting effect
  df.ant.filt.orient = df.ant.filt %>% filter(cue %in% c('center', 'spatial')) %>%
    group_by(cue, subject) %>%
    summarize(correct = mean(correct),
              rt = mean(rt)) %>%
    arrange(subject) %>%
    group_by(subject) %>%
    mutate(orienting = rt - lead(rt)) %>%
    filter(cue == 'center') %>%
    dplyr::select(subject, orienting)
  
  # executive effect
  df.ant.filt.exec = df.ant.filt %>% filter(flanker_type %in% c('congruent', 'incongruent')) %>%
    group_by(flanker_type, subject) %>%
    summarize(correct = mean(correct),
              rt = mean(rt)) %>%
    arrange(subject) %>%
    group_by(subject) %>%
    mutate(exec = rt - lag(rt)) %>%
    filter(flanker_type == 'incongruent') %>%
    dplyr::select(subject, exec)
  
  # combine them
  df.ant.networks = df.ant.filt.alert
  df.ant.networks$orienting = df.ant.filt.orient$orienting
  df.ant.networks$exec = df.ant.filt.exec$exec
  rm(df.ant.filt.alert, df.ant.filt.orient, df.ant.filt.exec)
}

# Clean up demographics df -----------------------------------------

### MODELS
models = c('Full', 'BinAtts', 'BinWts', 'BinWtsAtts', '1Att', '1AttBinAtts')
model.prob.names = paste0('model.', models, '.prob')
models.one.att = c('Multiple', 'Multiple', 'Multiple', 'Multiple', 'One', 'One')
models.bin.wts = c('Graded', 'Graded', 'Binary', 'Binary', NA, NA)
models.bin.atts = c('Graded', 'Binary', 'Graded', 'Binary', 'Graded', 'Binary')
models.order = c('Full', 'BinAtts', 'BinWts', 'BinWtsAtts', '1Att', '1AttBinAtts')
df.demo = df.demo %>%
  mutate(reported.h1 = 1 - lex_real / 100, # question was reversed in actual task
         reported.h2 = 1 - binwts_real / 100,
         reported.h3 = 1 - binatts_real / 100,
         reported.h1.dichotomized = reported.h1 > 0.5, # did they report using each heuristic
         reported.h2.dichotomized = reported.h2 > 0.5,
         reported.h3.dichotomized = reported.h3 > 0.5,
         reported.model.num = ifelse(reported.h1.dichotomized,
                                     ifelse(reported.h3.dichotomized, 6, 5),
                                     ifelse(reported.h2.dichotomized,
                                            ifelse(reported.h3.dichotomized, 4, 3),
                                            ifelse(reported.h3.dichotomized, 2, 1))),
         reported.model = models[reported.model.num],
         reported.model.fac = factor(reported.model, models.order, models.order),
         norm.h1 = 1 - lex_norm / 100,
         norm.h2 = 1 - binwts_norm / 100,
         norm.h3 = 1 - binatts_norm / 100,
         norm.h1.dichotomized = norm.h1 > 0.5,
         norm.h2.dichotomized = norm.h2 > 0.5,
         norm.h3.dichotomized = norm.h3 > 0.5,
         norm.model.num = ifelse(norm.h1.dichotomized,
                                 ifelse(norm.h3.dichotomized, 6, 5),
                                 ifelse(norm.h2.dichotomized,
                                        ifelse(norm.h3.dichotomized, 4, 3),
                                        ifelse(norm.h3.dichotomized, 2, 1))),
         norm.model = models[norm.model.num],
         norm.model.fac = factor(norm.model, models.order, models.order))

## ADD ANT
df.demo$alerting = NA
df.demo$orienting = NA
df.demo$exec = NA
for (i in 1:nrow(df.demo)) {
  if (exists('df.ant.networks')) {
    ant.row = df.ant.networks$subject == df.demo$subject[i]
    
    if (any(ant.row)) {
      alerting = df.ant.networks$alerting[ant.row]
      if (!is.null(alerting)) df.demo$alerting[i] = alerting;
      
      orienting = df.ant.networks$orienting[ant.row]
      if (!is.null(orienting)) df.demo$orienting[i] = orienting;
      
      exec = df.ant.networks$exec[ant.row]
      if (!is.null(exec)) df.demo$exec[i] = exec;
    }
  }
}

### ADD ICAR
for (i in 1:nrow(df.demo.followup)) {
  df.demo$icar_num_correct[df.demo$subject == df.demo.followup$subject[i]] = df.demo.followup$icar_num_correct[i]
}

### MEDITATION EXP
howoften_options = c('Less than once per week', 'About once per week', '2-4 times per week', 'Daily or almost daily')
amount_options = c('5-15 minutes per day', '15-30 minutes per day', '>30 minutes per day')
years_options = c('0-1 years', '1-3 years', '3-5 years', '5+ years')
df.demo = df.demo %>%
  mutate(meditation = factor(meditation_exp2, c('Yes', 'No', ''), c('Currently', 'Used to', 'Never')),
         meditation.over5 = ifelse(meditation_exp1 == 'No', F, meditation_exp2 == 'More than 5 years'),
         meditation_exp3 = as.numeric(meditation_exp3),
         meditation_exp4 = as.numeric(meditation_exp4))

### EDUCATION
df.demo = df.demo %>%
  mutate(edu.num = as.numeric(factor(edu,
                                     levels=c('Some high school', 'High school', 'Some college',
                                              '2 year degree', '4 year degree', 'Postgraduate/Professional degree/other'))))
df.demo = df.demo %>%
  mutate(choice_exp_num = as.numeric(factor(choice_exp, c("0", '5-Jan', 
                                                          '10-May', '15-Oct', '15+'))))
### SCALES
battery = c('decisionstyle', 'acs', 'mindfulness', 'sris', 'maia')
max.vals = c(5,4,4,6,6)

# which ones are reversed?
ds.reversed = 6:10
acs.reversed = c(1, 2, 3, 6, 7, 8, 11, 12, 15, 16, 20)
mindfulness.reversed = c(2,6,7)
sris.reversed = c(1, 2, 4, 7, 14, 16, 17, 18, 19)
maia.reversed = c()

reversed = list(ds.reversed, acs.reversed, mindfulness.reversed, sris.reversed, maia.reversed)

# which factors?
ds.factors = c(rep('deliberative',5), rep('intuitive',5))
acs.factors = c(rep('focusing',9), rep('shifting',11))
mindfulness.factors = c('attention', 'present focus', 'acceptance', 'acceptance', 'awareness',
                        'attention', 'present focus', 'awareness', 'awareness', 'acceptance', 'present focus', 'attention')
sris.factors = c(rep('tendency',12), rep('insight',8))
maia.factors = c(rep('noticing',4),rep('attention regulation',7),rep('emotional awareness',5),rep('body listening',3),rep('trusting',3))

factors = list(ds.factors, acs.factors, mindfulness.factors, sris.factors, maia.factors)

for (battery.ind in 1:length(battery)) {
  name = battery[battery.ind]
  reversed.cur = reversed[[battery.ind]]
  factors.cur = factors[[battery.ind]]
  factors.unique = unique(factors.cur)
  max.val.cur = max.vals[battery.ind]
  df.demo.followup[,name] = NA
  df.demo[,name] = NA
  
  for (i in 1:nrow(df.demo.followup)) {
    resp.str = df.demo.followup[i, paste0(name, '_responses')]
    if (!is.null(resp.str)) {
      demo.row = df.demo$subject == df.demo.followup$subject[i]
      
      resp = as.numeric.vector(resp.str)
      resp[reversed.cur] = max.val.cur - resp[reversed.cur] + 1
      df.demo.followup[i,name] = mean(resp)
      df.demo[demo.row,name] = mean(resp)
      
      for (fac in factors.unique) {
        which.items = which(factors.cur == fac)
        df.demo.followup[i,paste0(name, '.', fac)] = mean(resp[which.items])
        df.demo[demo.row,paste0(name, '.', fac)] = mean(resp[which.items])
      }
    }
  }
}

### COMPREHENSION CHECKS

demo.colnames = colnames(df.demo)
for (i in 1:nrow(df.demo)) {
  if (df.demo$lex_comp1_number[i] == 2) {
    demo.colnames.temp1 = demo.colnames[grepl("lex_comp1.*", demo.colnames)]
    demo.colnames.temp2 = demo.colnames[grepl("lex_comp2.*", demo.colnames)]
    values.temp = df.demo[i,demo.colnames.temp1]
    df.demo[i,demo.colnames.temp1] = df.demo[i,demo.colnames.temp2]
    df.demo[i,demo.colnames.temp2] = values.temp
  }
  if (df.demo$binwts_comp1_number[i] == 2) {
    demo.colnames.temp1 = demo.colnames[grepl("binwts_comp1.*", demo.colnames)]
    demo.colnames.temp2 = demo.colnames[grepl("binwts_comp2.*", demo.colnames)]
    values.temp = df.demo[i,demo.colnames.temp1]
    df.demo[i,demo.colnames.temp1] = df.demo[i,demo.colnames.temp2]
    df.demo[i,demo.colnames.temp2] = values.temp
  }
  if (df.demo$binatts_comp1_number[i] == 2) {
    demo.colnames.temp1 = demo.colnames[grepl("binatts_comp1.*", demo.colnames)]
    demo.colnames.temp2 = demo.colnames[grepl("binatts_comp2.*", demo.colnames)]
    values.temp = df.demo[i,demo.colnames.temp1]
    df.demo[i,demo.colnames.temp1] = df.demo[i,demo.colnames.temp2]
    df.demo[i,demo.colnames.temp2] = values.temp
  }
}

colnames.temp = c('lex_comp1', 'lex_comp2', 'binwts_comp1', 'binwts_comp2', 'binatts_comp1', 'binatts_comp2')
df.cc = df.demo %>% select(subject, strat_q_order, colnames.temp) %>% pivot_longer(colnames.temp)
df.cc = df.cc %>% mutate(answer = ifelse(grepl(".*1.*", name), 0, 100),
                         correct = (answer == 0 & value < 50) | (answer == 100 & value > 50))
ggplot(df.cc, aes(x = value)) +
  geom_histogram() +
  geom_vline(aes(xintercept = answer), color = 'red', linetype = 'dashed') +
  geom_vline(aes(xintercept = 50), color = 'black', linetype = 'dashed') +
  facet_wrap(~name)

df.cc.question = df.cc %>% group_by(name) %>%
  summarize(correct = mean(correct), correct.se = se.prop(correct))
df.cc.question

df.cc.order = df.cc %>% group_by(strat_q_order, name) %>%
  summarize(correct = mean(correct), correct.se = se.prop(correct))
df.cc.order %>% filter(name == 'binwts_comp1')

df.cc.subj = df.cc %>% group_by(subject) %>%
  summarize(correct.total = sum(correct),
            correct.without.binatts1 = sum(ifelse(name == 'binatts_comp1', 0, correct)),
            correct.lex = sum(ifelse(grepl('lex_comp.*', name), correct, 0)),
            correct.binwts = sum(ifelse(grepl('binwts_comp.*', name), correct, 0)),
            correct.binatts = sum(ifelse(grepl('binatts_comp.*', name), correct, 0)))

for (i in 1:nrow(df.demo)) {
  cc.row = df.cc.subj$subject == df.demo$subject[i]
  df.demo$cc.correct.total[i] = df.cc.subj$correct.total[cc.row]
  df.demo$cc.correct.lex[i] = df.cc.subj$correct.lex[cc.row]
  df.demo$cc.correct.binwts[i] = df.cc.subj$correct.binwts[cc.row]
  df.demo$cc.correct.binatts[i] = df.cc.subj$correct.binatts[cc.row]
}

# Clean up S1 --------------------------------------------------------------
movie_scale = c('Very Bad', 'Bad', 'Moderate', 'Good', 'Very Good')

## Stage 1 choices
# make wide df, put each attribute in its own column
if (participant_type == 'decider') {
  atts = unique(df.attributes$attribute)
} else {
  atts = read.csv(paste0(filepath, 'att_order.csv'), header = F)
  atts = atts$V1
}

numAtts = length(atts)
att.nums = 1:length(atts)
att.nums.str = as.character(att.nums)

atts.opt1 = paste0(atts,'.opt1')
atts.opt2 = paste0(atts,'.opt2')
atts.opt1.enclosed = paste0('`',atts.opt1,'`')
atts.opt2.enclosed = paste0('`',atts.opt2,'`')
atts.opt1.scaled = paste0(atts.opt1,'.scaled')
atts.opt2.scaled = paste0(atts.opt2,'.scaled')
atts.opt.diff = paste0(atts,'.diff')
atts.opt.scaled.diff = paste0(atts,'.scaled.diff')
atts.opt.diff.enclosed = paste0('`',atts.opt.diff,'`')

df.s1[,atts.opt1] = NA
df.s1[,atts.opt2] = NA

df.avail.atts = data.frame(matrix(1,nrow = nrow(df.s1), ncol = length(atts)+1))
colnames(df.avail.atts) = c('subject.num', atts.opt1)
df.avail.atts$subject.num = df.s1$subject.num

for (i in 1:nrow(df.s1)) {
  cur.atts = as.string.vector(df.s1$attributes[i])
  cur.opt1.vals = as.string.vector(df.s1$opt1_values[i])
  cur.opt2.vals = as.string.vector(df.s1$opt2_values[i])
  
  cur.att.nums = numeric(length(cur.atts))
  cur.atts.order = numeric(length(cur.atts))
  
  for (j in 1:length(cur.atts)) {
    cur.att.nums[j] = which(cur.atts[j] == atts)
    
    if (task_variant == 'movie') {
      first = floor(as.numeric(cur.opt1.vals[j])/20)
      if (first < 0) {
        first = 0;
      } else if (first > 4) {
        first = 4;
      }
      second = floor(as.numeric(cur.opt2.vals[j])/20)
      if (second < 0) {
        second = 0;
      } else if (second > 4) {
        second = 4;
      }
      
      df.s1[i,atts.opt1[cur.att.nums[j]]] = movie_scale[first+1]
      df.s1[i,atts.opt2[cur.att.nums[j]]] = movie_scale[second+1]
    } else {
      df.s1[i,atts.opt1[cur.att.nums[j]]] = cur.opt1.vals[j]
      df.s1[i,atts.opt2[cur.att.nums[j]]] = cur.opt2.vals[j]
    }
    
    #df.avail.atts[i,atts.opt1[cur.att.nums[j]]] = 1
    cur.atts.order[j] = which(atts[j] == cur.atts)
  }
  
  df.s1$atts.order[i] = as.string(cur.atts.order)
  df.s1$att.nums[i] = as.string(cur.att.nums)
}

# convert attribute values to numerics
for (i in 1:length(atts)) {
  cur.att.opt1 = atts.opt1[i]
  cur.att.opt2 = atts.opt2[i]
  
  cur.scale = if(task_variant == 'movie') movie_scale else unique(df.attributes$scale[df.attributes$attribute == atts[i]])
  
  if (any(is.na(cur.scale)) || any(cur.scale == "")) {
      if (task_variant == 'movie') {
      df.s1[,cur.att.opt1] = as.numeric(sub('\\%.*', '', df.s1[,cur.att.opt1]))
      df.s1[,cur.att.opt2] = as.numeric(sub('\\%.*', '', df.s1[,cur.att.opt2]))
    } else {
      df.s1[,cur.att.opt1] = as.numeric(sub('\\ .*', '', df.s1[,cur.att.opt1]))
      df.s1[,cur.att.opt2] = as.numeric(sub('\\ .*', '', df.s1[,cur.att.opt2]))
    }
  } else {
    cur.scale = unlist(strsplit(cur.scale, split = ","))
    df.s1[,cur.att.opt1] = as.numeric(factor(df.s1[,cur.att.opt1], cur.scale))
    df.s1[,cur.att.opt2] = as.numeric(factor(df.s1[,cur.att.opt2], cur.scale))
  }
}

# compute diffs
for (i in 1:length(atts)) {
  df.s1[,atts.opt1.scaled[i]] = rescale(df.s1[,atts.opt1[i]])
  df.s1[,atts.opt2.scaled[i]] = rescale(df.s1[,atts.opt2[i]])
  df.s1[,atts.opt.diff[i]] = df.s1[,atts.opt2[i]] - df.s1[,atts.opt1[i]]
  df.s1[,atts.opt.scaled.diff[i]] = df.s1[,atts.opt2.scaled[i]] - df.s1[,atts.opt1.scaled[i]]
}

df.s1.subj = df.s1 %>% 
  mutate(correct.choice = choice == orig_choice) %>% # for observer studies
  group_by(subject) %>%
  summarize(total.time = sum(rt) / 60000,
            pct_left = mean(choice == 0),
            median_rt = median(rt),
            sd_rt = sd(rt),
            num_trials = n(),
            pct_correct = mean(correct.choice))

# Clean up S2 / attributes -------------------------------------------------------------

# This is because of a bug in the movie version
if (task_variant == 'movie') {
  good.names = c('artistic', 'funny', 'romantic', 'visual', 'great acting', 'plot', 'dialogue', 'awesome soundtrack', 'good action')
  bad.names = c('Creativity', 'Humor', 'Romantic Scenes', 'Visuals', 'Acting', 'Plot', 'Dialogue', 'Soundtrack', 'Action Scenes')
  for (i in 1:nrow(df.s2)) {
    for (j in 1:length(good.names)) {
      if (df.s2$attribute[i] == bad.names[j]) {
        df.s2$attribute[i] = good.names[j]
      }
    }
  }
}

for (i in 1:nrow(df.attributes)) {
  cur.att = df.attributes$attribute[i]
  subj = df.attributes$subject[i]
  cur.scale = ifelse(task_variant == 'movie', movie_scale,df.attributes$scale[i])
  
  rows.demo = as.character(df.demo$subject) == as.character(subj)
  rows.s2.wad = as.character(df.s2$subject) == as.character(subj) & df.s2$attribute == cur.att & df.s2$type == 'wad_att_rating'
  rows.s2.ew = as.character(df.s2$subject) == as.character(subj) & df.s2$attribute == cur.att & df.s2$type == 'ew_att_rating'
  rows.s2.lex = as.character(df.s2$subject) == as.character(subj) & df.s2$attribute == cur.att & df.s2$type == 'lex_att_rating'
  rows.s2.direction = as.character(df.s2$subject) == as.character(subj) & df.s2$attribute == cur.att & df.s2$type == 'direction'
  
  # if this is the one...
  if (any(rows.s2.lex)) {
    df.attributes$reported.weight.single[i] = 1 # should be 1
  } else {
    df.attributes$reported.weight.single[i] = 0
  }
  
  if (any(rows.s2.ew)) {
    df.attributes$reported.weight.binary[i] = 1 # should be 1
  } else {
    df.attributes$reported.weight.binary[i] = 0
  }
  
  if (any(rows.s2.wad)) {
    df.attributes$reported.weight.graded[i] = df.s2$rating[rows.s2.wad] / 100
  } else {
    df.attributes$reported.weight.graded[i] = NA
  }
  
  if (any(rows.s2.direction)) {
    if (version %in% c('home-fixed', 'home-fixed-observers')) {
      if (is.na(cur.scale) || cur.scale == "") {
        least = as.numeric(sub('\\ .*', '', df.s2$least_preferred[rows.s2.direction]))
        most = as.numeric(sub('\\ .*', '', df.s2$most_preferred[rows.s2.direction]))
        bounds = c(df.attributes$lb[i], df.attributes$ub[i])
      } else {
        cur.scale = unlist(strsplit(cur.scale, ","))
        least = which(df.s2$least_preferred[rows.s2.direction] == cur.scale)
        most = which(df.s2$most_preferred[rows.s2.direction] == cur.scale)
        bounds = c(1, length(cur.scale))
      }
      
      df.attributes$most[i] = most
      df.attributes$least[i] = least
      df.attributes$direction[i] = sign(most - least)
      df.attributes$same[i] = most == least
      df.attributes$linear[i] = least %in% bounds & most %in% bounds
    } else {
      df.attributes$most[i] = NA
      df.attributes$least[i] = NA
      df.attributes$direction[i] = ifelse(df.s2$rating[rows.s2.direction][1] == 1, 1, -1)
      df.attributes$same[i] = F
      df.attributes$linear[i] = all(df.s2$rating[rows.s2.direction] == 1) | all(df.s2$rating[rows.s2.direction] == 0)
    }
  } else {
    df.attributes$direction[i] = NA
    df.attributes$same[i] = NA
    df.attributes$linear[i] = NA
  }
  
  df.attributes$reported.h1[i] = df.demo$reported.h1[rows.demo] 
  df.attributes$reported.h2[i] = df.demo$reported.h2[rows.demo]  
  df.attributes$reported.h3[i] = df.demo$reported.h3[rows.demo] 
  df.attributes$reported.h1.dichotomized[i] = df.demo$reported.h1.dichotomized[rows.demo] 
  df.attributes$reported.h2.dichotomized[i] = df.demo$reported.h2.dichotomized[rows.demo]  
  df.attributes$reported.h3.dichotomized[i] = df.demo$reported.h3.dichotomized[rows.demo] 
  df.attributes$reported.model[i] = df.demo$reported.model[rows.demo]
  df.attributes$reported.model.num[i] = df.demo$reported.model.num[rows.demo]
  df.attributes$reported.model.fac[i] = df.demo$reported.model.fac[rows.demo]
  df.attributes$reported.weight.reported[i] = ifelse(df.demo$reported.h1.dichotomized[rows.demo] == 'One',
                                                     df.attributes$reported.weight.single[i],
                                                     ifelse(df.demo$reported.h2.dichotomized[rows.demo] == 'Binary',
                                                            df.attributes$reported.weight.binary[i],
                                                            df.attributes$reported.weight.graded[i]))
}

df.attributes = df.attributes %>%
  mutate(reported.weight.graded.signed = reported.weight.graded * direction,
         reported.weight.binary.signed = reported.weight.binary * direction,
         reported.weight.single.signed = reported.weight.single * direction,
         reported.weight.reported.signed = reported.weight.reported * direction)

## Save some data that we need to analyze the observer data
# (this gets saved to the observers folder)
if (has_observers) {
  write.table(df.demo %>% select(subject, subject.num), paste0(filepath_observer, 'observer_mapping.csv'), row.names = F, col.names = F, sep = ",")
  write.table(atts, paste0(filepath_observer, 'att_order.csv'), row.names = F, col.names = F, sep = ",")
}

# Import modeling results ------------------------------------------
# Model-fitting is done in "fit-subjects.R"
if (participant_type == 'decider') {
  fitting_results = combine_lists(paste0(filepath, 'modeling-output/'))

  decider.names = df.demo$subject
  decider.nums = df.demo$subject.num
} else { # observer stuff
  fitting_results = combine_lists(paste0(filepath_decider, 'modeling-output/'))
  
  ## Load the subject mapping data (which should have been saved from the decider analysis)
  df.map = read.csv(paste0(filepath, 'observer_mapping.csv'), header = F)
  colnames(df.map) = c('subject', 'subject.num')
  
  decider.names = df.map$subject
  decider.nums = df.map$subject.num
}

## model comparison

# values to extract:
# best model for each subject, posterior probability of that model, lme of that model
# posterior & lme for reported model
# posterior for each model
# prob for each of the three properties
# loo results
# diagnostic results

numHeuristics <- 3
heuristic_families <- list(
  list(1:4, 5:6),
  list(1:2, 3:4),
  list(c(1,3,5), c(2,4,6))
)

df.demo$actual.model.num = NA
df.demo$actual.model.prob = NA
df.demo$actual.model.lme = NA
df.demo$actual.inv.temp = NA
df.demo$actual.model.loo = NA
df.demo$actual.model.diagnostics = NA
df.demo$process.accuracy.model.prob = NA
df.demo$reported.model.lme = NA
df.demo$norm.model.prob = NA
numSubj = length(subjlist)
model_weights = vector(mode='list', length = numSubj)
for (s in 1:numSubj) {
  decider.num = ifelse(participant_type == 'decider',
                              s,
                              df.demo$target_id_num[s])
  
  if (length(decider.num) > 0) {
    df.demo$decider.subj.num[s] = decider.num
    
    model_temps = rep(NA,length(models))
    model_weights[[s]] = matrix(NA, nrow = length(models), ncol = numAtts)
    colnames(model_weights[[s]]) = atts
    model_probs = rep(NA,length(models))
    model_lmes = rep(NA,length(models))
    for (m in 1:length(models)) {
      if (!is.null(fitting_results[[decider.num]][[1]])) {
        params_cur = fitting_results[[decider.num]][[m]][[1]]
        model_temps[m] = params_cur[1]
        model_weights[[s]][m,] = params_cur[2:(numAtts+1)]
        model_lmes[m] = fitting_results[[decider.num]][[m]][[2]]
      }
    }
    
    model_mes = exp(model_lmes)
    model_probs = model_mes / sum(model_mes)
    
    for (m in 1:length(models)) {
      df.demo[s,paste0('model.', models[m], '.prob')] = model_probs[m]
    }
    
    if (any(!is.na(model_probs))) {
      # stats of best-fitting model
      df.demo$actual.model.num[s] = which.max(model_probs)
      df.demo$actual.model.prob[s] = max(model_probs)
      df.demo$actual.model.lme[s] = model_lmes[df.demo$actual.model.num[s]]
      df.demo$actual.inv.temp[s] = model_temps[df.demo$actual.model.num[s]]
      df.demo$actual.model.loo[s] = fitting_results[[decider.num]][[df.demo$actual.model.num[s]]][[5]]$estimates[1,1]
      df.demo$actual.model.diagnostics[s] = fitting_results[[decider.num]][[df.demo$actual.model.num[s]]][4]
      
      # stats of reported model
      df.demo$process.accuracy.model.prob[s] = model_probs[df.demo$reported.model.num[s]]
      df.demo$reported.model.lme[s] = model_lmes[df.demo$reported.model.num[s]]
      
      # stats of norm model
      df.demo$norm.model.prob[s] = model_probs[df.demo$reported.model.num[s]]
    }
    
    # property questions
    for (heuristic in 1:numHeuristics) {
      df.demo[s, paste0('actual.h', heuristic, '.prob')] = mean(model_mes[heuristic_families[[heuristic]][[2]]]) / (mean(model_mes[heuristic_families[[heuristic]][[1]]]) + mean(model_mes[heuristic_families[[heuristic]][[2]]]))
    }
  } else {
    df.demo$decider.subj.num[s] = NA
  }
}

ll.chance = log(.5 ^ numTrials)

df.demo = df.demo %>%
  mutate(actual.model = models[actual.model.num],
         actual.model.fac = factor(actual.model, models.order, models.order),
         process.accuracy.model.dichotomized = reported.model.num == actual.model.num,
         norm.model.dichotomized = norm.model.num == actual.model.num,
         process.accuracy.div.h1 = abs(actual.h1.prob - reported.h1),
         process.accuracy.div.h2 = abs(actual.h2.prob - reported.h2),
         process.accuracy.div.h3 = abs(actual.h3.prob - reported.h3),
         actual.h1.prob.dichotomized = actual.h1.prob > 0.5,
         actual.h2.prob.dichotomized = actual.h2.prob > 0.5,
         actual.h3.prob.dichotomized = actual.h3.prob > 0.5,
         process.accuracy.h1.dichotomized = reported.h1.dichotomized == actual.h1.prob.dichotomized,
         process.accuracy.h2.dichotomized = reported.h2.dichotomized == actual.h2.prob.dichotomized,
         process.accuracy.h3.dichotomized = reported.h3.dichotomized == actual.h3.prob.dichotomized,
         process.accuracy.model.relprob = process.accuracy.model.prob / actual.model.prob,
         norm.model.relprob = norm.model.prob / actual.model.prob,
         actual.model.loo.readable = exp(actual.model.loo/numTrials)) %>%
  rowwise() %>% mutate(
    process.accuracy.div = mean(c(process.accuracy.div.h1, process.accuracy.div.h2, process.accuracy.div.h3)),
    process.accuracy.qs.dichotomized = mean(c(process.accuracy.h1.dichotomized, process.accuracy.h2.dichotomized, process.accuracy.h3.dichotomized)),
  )

## weights

for (i in 1:nrow(df.attributes)) {
  rows.demo = as.character(df.demo$subject) == as.character(df.attributes$subject[i])
  
  df.attributes$actual.model.num[i] = df.demo$actual.model.num[rows.demo]
  df.attributes$reported.weight.actual[i] = ifelse(df.demo$actual.h1.prob.dichotomized[rows.demo],
                                                   df.attributes$reported.weight.single[i],
                                                   ifelse(df.demo$actual.h2.prob.dichotomized[rows.demo],
                                                          df.attributes$reported.weight.binary[i],
                                                          df.attributes$reported.weight.graded[i]))

  df.attributes$fitted.weight.reported[i] = model_weights[[df.attributes$subject.num[i]]][df.attributes$reported.model.num[i],df.attributes$attribute[i]]
  df.attributes$fitted.weight.actual[i] = model_weights[[df.attributes$subject.num[i]]][df.attributes$actual.model.num[i],df.attributes$attribute[i]]
  
  model.prob = numeric(length(models))
  model.att = numeric(length(models))
  for (m in 1:length(models)) {
    df.attributes[i, paste0('fitted.weight.', models[m])] = model_weights[[df.attributes$subject.num[i]]][m,df.attributes$attribute[i]]
    model.att[m] = df.attributes[i, paste0('fitted.weight.', models[m])]
    model.prob[m] = df.demo[rows.demo, paste0('model.', models[m], '.prob')][[1]]
  }
  
  df.attributes$fitted.weight.averaged[i] = sum(model.prob * model.att)
  df.attributes$reported.weight.averaged[i] = sum(model.prob * c(df.attributes$reported.weight.graded[i], df.attributes$reported.weight.graded[i],
                                                                  df.attributes$reported.weight.binary[i], df.attributes$reported.weight.binary[i],
                                                                  df.attributes$reported.weight.single[i], df.attributes$reported.weight.single[i]))
}

df.attributes = df.attributes %>% mutate(reported.weight.actual.signed = reported.weight.actual * direction,
                                         reported.weight.averaged.signed = reported.weight.averaged * direction,
                                         reported.weight.reported.signed = reported.weight.reported * direction)

# get subject-level accuracies
df.attributes.subj = df.attributes %>% group_by(subject) %>%
  summarize(num_same = sum(same & (reported.weight.graded > .1), na.rm = T),
            weight.accuracy.reported = cor(fitted.weight.reported, reported.weight.reported.signed),
            weight.accuracy.reported.rank = cor(fitted.weight.reported, reported.weight.reported.signed, method = 'kendall'),
            weight.accuracy.best = cor(fitted.weight.actual, reported.weight.actual.signed),
            weight.accuracy.best.rank = cor(fitted.weight.actual, reported.weight.actual.signed, method = 'kendall'),
            weight.accuracy.averaged = cor(fitted.weight.averaged, reported.weight.averaged.signed),
            weight.accuracy.averaged.rank = cor(fitted.weight.averaged, reported.weight.averaged.signed, method = 'kendall')
  )

# transfer these to df.demo

cols.to.xfer = colnames(df.attributes.subj %>% select(!subject))
for (i in 1:nrow(df.demo)) {
  att.row = df.attributes.subj$subject == df.demo$subject[i]
  
  for (cur.col in cols.to.xfer) {
    df.demo[i,cur.col] = df.attributes.subj[att.row, cur.col]
  }
}

# Normative alignment stuff
df.demo = df.demo %>% mutate(
  norm.h1.actual.div = abs(actual.h1.prob - norm.h1),
  norm.h2.actual.div = abs(actual.h2.prob - norm.h2),
  norm.h3.actual.div = abs(actual.h3.prob - norm.h3),
  norm.h1.actual.dichotomized = norm.h1.dichotomized == actual.h1.prob.dichotomized,
  norm.h2.actual.dichotomized = norm.h2.dichotomized == actual.h2.prob.dichotomized,
  norm.h3.actual.dichotomized = norm.h3.dichotomized == actual.h3.prob.dichotomized) %>%
  rowwise() %>% mutate(
    norm.actual.div = mean(c(norm.h1.actual.div, norm.h2.actual.div, norm.h3.actual.div)),
    norm.actual.dichotomized = mean(c(norm.h1.actual.dichotomized, norm.h2.actual.dichotomized, norm.h3.actual.dichotomized))
  )

rm(fitting_results)

# Browser events ----------------------------------------------------------

df.browser = read.csv(paste0(filepath_data,'browser_events.csv'), stringsAsFactors = F)
df.browser.subj = df.browser %>%
  group_by(subject) %>%
  summarize(num.blur = sum(browser_event == 'blur'))

# Get people who did two versions -------------------------------------------
df.demo$did_two = NA
df.demo$which_version_first = NA
demo_movie = read.csv('movie/data/demo.csv') %>%
  mutate(stamp = parse_date_time(stamp, c('ymd HMS', 'mdy HMS', 'mdy HM')))
demo_rand = read.csv('home-random/data/demo.csv') %>%
  mutate(stamp = parse_date_time(stamp, c('ymd HMS', 'mdy HMS', 'mdy HM')))

for (i in 1:nrow(df.demo)) {
  subj = df.demo$subject[i]

  if (subj %in% demo_movie$subject && subj %in% demo_rand$subject) {
    df.demo$did_two[i] = T
    movie_stamp = demo_movie$stamp[demo_movie$subject == subj]
    rand_stamp = demo_rand$stamp[demo_rand$subject == subj]
    if (movie_stamp < rand_stamp) {
      df.demo$which_version_first[i] = 'movie'
    } else {
      df.demo$which_version_first[i] = 'home-random'
    }
  } else {
    df.demo$did_two[i] = F
  }
}

# Do filtering ---------------------------------------------------------

exclude.subj = c()

for (subj in subjlist) {
  demo.row = df.demo$subject == subj
  s1.subj = df.s1.subj$subject == subj
  browser.row = df.browser.subj$subject == subj
  
  if (df.demo$instruction_times_median[demo.row] < 2 || # if they took, on average, <2 sec per instruction screen
      df.s1.subj$num_trials[s1.subj] != numTrials | # if they didn't complete all the trials
      df.demo$attention[demo.row] < 50 | # if they reported themselves as having paid less than 50% attention
      df.demo$cc.correct.total[demo.row] < 4 | # if they got fewer than 4/6 comprehension checks correct
      (df.demo$did_two[demo.row] & df.demo$which_version_first[demo.row] != version) |
      df.attributes.subj$num_same[df.attributes.subj$subject == subj] > 0 |
      #df.demo$cv_result[demo.row] < .5 |
      (any(browser.row) & df.browser.subj$num.blur[df.browser.subj$subject == subj] > 20)) { # if they tabbed away from the experiment more than 20 times
    exclude.subj = c(exclude.subj, subj)
  }
  
  if (participant_type == 'decider') {
    if (df.s1.subj$pct_left[s1.subj] > .8 | # if they chose the left or right option >80% of the time
        df.s1.subj$pct_left[s1.subj] < .2) {
      exclude.subj = c(exclude.subj, subj)
    }
  } else {
    if (df.s1.subj$pct_correct[s1.subj] < 0.95 |
        is.na(df.demo$decider.subj.num[demo.row])) {
      exclude.subj = c(exclude.subj, subj)
    }
  }
  
}

df.demo.filt = df.demo %>% filter(!(subject %in% exclude.subj))
df.s1.filt = df.s1 %>% filter(!(subject %in% exclude.subj))
df.s1.subj.filt = df.s1.subj %>% filter(!(subject %in% exclude.subj))
df.s1.practice.filt = df.s1.practice %>% filter(!(subject %in% exclude.subj))
df.s2.filt = df.s2 %>% filter(!(subject %in% exclude.subj))
df.attributes.filt = df.attributes %>% filter(!(subject %in% exclude.subj))
df.attributes.subj.filt = df.attributes.subj %>% filter(!(subject %in% exclude.subj))
df.cc.filt = df.cc %>% filter(!(subject %in% exclude.subj))
df.cc.subj.filt = df.cc.subj %>% filter(!(subject %in% exclude.subj))

# Data checks ----------------------------------------------------------

## linearity
mean(df.attributes.filt$linear[df.attributes.filt$reported.weight.graded > .1])
df.attributes.byatt = df.attributes.filt %>% filter(reported.weight.graded > .1) %>%
  group_by(attribute) %>%
  summarize(linear = mean(linear))
df.attributes.byatt

## appropriateness
ggplot(df.demo.filt, aes(x = appropriateness)) +
  geom_histogram(color = 'white', bins = 10) +
  theme_black() +
  #geom_vline(xintercept = 100, color = 'red', linetype = 'dashed') +
  #geom_vline(xintercept = 50, color = 'red', linetype = 'dashed') +
  geom_vline(xintercept = mean(df.demo.filt$appropriateness), color = 'red') +
  geom_vline(xintercept = mean(df.demo.filt$appropriateness)+se(df.demo.filt$appropriateness), color = 'red', linetype = 'dashed') +
  geom_vline(xintercept = mean(df.demo.filt$appropriateness)-se(df.demo.filt$appropriateness), color = 'red', linetype = 'dashed') +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100), limits = c(-10,110)) +
  labs(x = '\nRating', y = '# of subjects')
get.ci(df.demo.filt$appropriateness)

## comprehension checks
ggplot(df.cc.filt, aes(x = value)) +
  geom_histogram() +
  geom_vline(aes(xintercept = answer), color = 'red', linetype = 'dashed') +
  geom_vline(aes(xintercept = 50), color = 'black', linetype = 'dashed') +
  facet_wrap(~name) +
  xlab('Response (0-100)') + ylab('# of Subjects') +
  scale_x_continuous(limits = c(-10,110), breaks = c(0,50,100)) +
  scale_y_continuous()

hist(df.cc.subj.filt$correct.total)

df.cc.question.filt = df.cc.filt %>% group_by(name) %>%
  summarize(correct = mean(correct), correct.se = se.prop(correct))
df.cc.question.filt

df.cc.question = df.cc %>% group_by(name) %>%
  summarize(correct = mean(correct), correct.se = se.prop(correct))
df.cc.question

## model fits
# which models were the best fits?
ggplot(df.demo.filt, aes(x = actual.model.fac)) +
  geom_bar(position = 'dodge', color = 'white') +
  theme_black() +
  labs(x = '\nBest-fitting model', y = '# of subjects')

ggplot(df.demo.filt, aes(x = actual.h1.prob)) +
  geom_histogram(color='white') +
  theme_black() +
  labs(x = '\nProbability that subject\nused heuristic', y = '# of subjects')
ggplot(df.demo.filt, aes(x = actual.h2.prob)) +
  geom_histogram(color='white') +
  theme_black() +
  labs(x = '\nProbability that subject\nused heuristic', y = '# of subjects')
ggplot(df.demo.filt, aes(x = actual.h3.prob)) +
  geom_histogram(color='white') +
  theme_black() +
  labs(x = '\nProbability that subject\nused heuristic', y = '# of subjects')

ggplot(df.attributes.filt, aes(x = fitted.weight.Full)) +
  geom_histogram(color='white') +
  theme_black() +
  labs(x = '\nAttribute weight', y = '# of attributes') +
  scale_y_continuous(breaks = NULL) +
  scale_x_continuous(limits = c(-1,1),
                     breaks = c(-1, 0, 1))

# how good were the model fits?
df.probs = df.demo.filt %>%
  select(subject.num, all_of(model.prob.names), actual.model.loo.readable, reported.model.fac)
#df.posts.filt$reported.model.fac = df.demo.filt$reported.model.fac
#df.posts.filt$r.sq = df.demo.filt$actual.model.r2
df.probs.long = df.probs %>% 
  mutate(subject.num = factor(subject.num),
         subject.num = fct_reorder(subject.num, actual.model.loo.readable)) %>%
  pivot_longer(!c(subject.num,actual.model.loo.readable,reported.model.fac), names_to = 'model', values_to = 'prob') %>%
  arrange(subject.num, desc(prob)) %>%
  mutate(reported.model = reported.model.fac == model)
ggplot(df.probs.long,
       aes(x = prob, y = subject.num, group = subject.num, color = subject.num)) +
  geom_point(size = 2) +
  guides(group = 'none', color = 'none') +#,
  #alpha = guide_legend(title = 'Parameter accuracy\n(higher is better)')) +
  scale_color_brewer(palette = 'Set3') +
  theme_black()+
  scale_y_discrete(expand = c(-.5,0)) +
  scale_colour_manual(values=rep(brewer.pal(9,"Set1"),times=100))+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  labs(x = 'Prob',
       y = 'Subject')

ggplot(df.probs.long,
       aes(x = subject.num, y = prob, group = subject.num, color = subject.num)) +
  geom_point(size = 2) +
  guides(group = 'none', color = 'none') +#,
  #alpha = guide_legend(title = 'Parameter accuracy\n(higher is better)')) +
  scale_color_brewer(palette = 'Set3') +
  theme_black()+
  scale_x_discrete(expand = c(-.5,0)) +
  scale_colour_manual(values=rep(brewer.pal(9,"Set1"),times=100))+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(x = 'Subject',
       y = 'Probability of model')

ggplot(df.demo.filt, aes(x = actual.model.loo.readable)) +
  geom_histogram(color = 'white', bins = 20) +
  labs(x = "\nR-squared of best-fitting model",
       y = "# of subjects") +
  geom_vline(xintercept = mean(df.demo.filt$actual.model.loo.readable), color = 'red') +
  geom_vline(xintercept = mean(df.demo.filt$actual.model.loo.readable) + se(df.demo.filt$actual.model.loo.readable), color = 'red', linetype = 'dashed') +
  geom_vline(xintercept = mean(df.demo.filt$actual.model.loo.readable) - se(df.demo.filt$actual.model.loo.readable), color = 'red', linetype = 'dashed') +
  scale_y_continuous(breaks = NULL) +
  theme_black()
get.ci(df.demo.filt$actual.model.loo.readable)

# how certain were the model fits?
df.certainty = df.probs.long %>% group_by(subject.num) %>% filter(row_number() <= 2) %>%
  mutate(diff = prob - lead(prob, default = last(prob))) %>% filter(row_number() < 2) %>%
  arrange(subject.num)
df.demo.filt$actual.model.difftoptwo = df.certainty$diff
ggplot(df.demo.filt, aes(x = actual.model.difftoptwo)) +
  geom_histogram(color = 'white', bins = 20) +
  labs(x = "\nProbability difference b/w\ntop 2 models",
       y = "# of subjects") +
  theme_black() +
  geom_vline(xintercept = mean(df.demo.filt$actual.model.difftoptwo), color = 'red') +
  geom_vline(xintercept = mean(df.demo.filt$actual.model.difftoptwo) + se(df.demo.filt$actual.model.difftoptwo), color = 'red', linetype = 'dashed') +
  geom_vline(xintercept = mean(df.demo.filt$actual.model.difftoptwo) - se(df.demo.filt$actual.model.difftoptwo), color = 'red', linetype = 'dashed')

ggplot(df.demo.filt, aes(x = actual.model.prob)) +
  geom_histogram(color = 'white', bins = 20) +
  labs(x = "\nProbability of best-fitting model",
       y = "# of subjects") +
  theme_black()+
  geom_vline(xintercept = 1/6, color = 'gray') +
  geom_vline(xintercept = mean(df.demo.filt$actual.model.prob), color = 'red') +
  geom_vline(xintercept = mean(df.demo.filt$actual.model.prob) + se(df.demo.filt$actual.model.prob), color = 'red', linetype = 'dashed') +
  geom_vline(xintercept = mean(df.demo.filt$actual.model.prob) - se(df.demo.filt$actual.model.prob), color = 'red', linetype = 'dashed')
get.ci(df.demo.filt$actual.model.prob)

# relationship b/w model fits and other things
ggplot(df.demo.filt, aes(x = actual.model.loo.readable, y = appropriateness)) +
  geom_point(color = 'gray') +
  geom_smooth(method='lm', color = 'white') +
  theme_black() +
  labs(x = 'Best-fitting model R^2', y = 'Rating')
summary(lm(scale(actual.model.loo.readable) ~ scale(appropriateness), data = df.demo.filt))

ggplot(df.demo.filt, aes(x = actual.model.loo.readable, y = consistency1)) +
  geom_point(color = 'gray') +
  geom_smooth(method='lm', color = 'white') +
  theme_black() +
  labs(x = 'Best-fitting model R^2',
       y = 'Self-reported consistency\nin heuristic use')
df.demo.filt %$% cor.test(actual.model.r2, consistency1)

ggplot(df.demo.filt, aes(x = actual.model.loo.readable, y = consistency2)) +
  geom_point(color = 'gray') +
  geom_smooth(method='lm', color = 'white') +
  theme_black() +
  labs(x = 'Best-fitting model R^2',
       y = 'Self-reported consistency\nin weights')
df.demo.filt %$% cor.test(actual.model.r2, consistency2)

ggplot(df.demo.filt, aes(x = actual.inv.temp, y = actual.model.loo.readable)) +
  geom_point(color = 'gray') +
  geom_smooth(method='lm', color = 'white') +
  theme_black() +
  labs(x = 'Inverse temperature\n(i.e., model-estimated noisiness)',
       y = 'Best-fitting model R^2')
summary(lm(scale(actual.model.r2) ~ scale(appropriateness), data = df.demo.filt))

ggplot(df.demo.filt, aes(x = actual.model.loo.readable, y = process.accuracy.div)) +
  geom_point(color = 'gray') +
  geom_smooth(method='lm', color = 'white') +
  theme_black() +
  labs(x = '\nRating', y = 'Out-of-sample likelihood\nof best model')
ggplot(df.demo.filt, aes(x = actual.model.loo.readable, y = weight.accuracy.averaged)) +
  geom_point(color = 'gray') +
  geom_smooth(method='lm', color = 'white') +
  theme_black() +
  labs(x = '\nRating', y = 'Out-of-sample likelihood\nof best model')

summary(lm(process.accuracy.div ~ scale(actual.model.loo.readable), data = df.demo.filt))
summary(lm(weight.accuracy.averaged ~ scale(actual.model.loo.readable), data = df.demo.filt))

# Heuristic awareness -------------------------------------------------------
## what models did people report using?
ggplot(df.demo.filt, aes(x = reported.h1)) +
  geom_histogram(color='white') +
  theme_black() +
  labs(x = '\nReported probability\nof using heuristic', y = '# of subjects')
ggplot(df.demo.filt, aes(x = reported.h2)) +
  geom_histogram(color='white') +
  theme_black() +
  labs(x = '\nReported probability\nof using heuristic', y = '# of subjects')
ggplot(df.demo.filt, aes(x = reported.h3)) +
  geom_histogram(color='white') +
  theme_black() +
  labs(x = '\nReported probability\nof using heuristic', y = '# of subjects')

# get version w/o single-att people
df.demo.filt.multiatt = df.demo.filt %>% filter(actual.h1.prob.dichotomized == 'Multiple')
df.attributes.filt.multiatt = df.attributes.filt %>% filter(subject %in% df.demo.filt.multiatt$subject)

## heat map comparison
df.demo.heat = df.demo.filt %>% group_by(reported.model.fac, actual.model.fac) %>%
  summarize(num.subj = n())
ggplot(df.demo.heat, aes(x = actual.model.fac, y = reported.model.fac,
                         fill = num.subj)) +
  geom_tile() +
  labs(y = '\nSelf-reported model', x = 'Best-fitting model') +
  #scale_fill_brewer(palette = 'YlOrRd') +
  guides(fill = guide_colorbar(title = '# of subjects')) +
  theme_black()
df.demo.heat.normed = df.demo.filt %>% group_by(reported.model.fac, actual.model.fac) %>%
  summarize(num.subj = n()) %>%
  group_by(reported.model.fac) %>%
  mutate(num.subj.norm = num.subj / sum(num.subj))
ggplot(df.demo.heat.normed, aes(x = actual.model.fac, y = reported.model.fac,
                                fill = num.subj.norm)) +
  geom_tile() +
  geom_text(aes(label = round(num.subj.norm, 2))) +
  labs(y = '\nSelf-reported model', x = 'Best-fitting model') +
  #scale_fill_brewer(palette = 'YlOrRd') +
  scale_fill_continuous(limits = c(0,1), low = 'black', high = 'white') +
  guides(fill = 'none') +
  theme_black()

## How good or bad were the models that subjects reported using?

# rel prob
set.seed(12345)
process.accuracy.model.prob.rand = numeric(1000)
process.accuracy.qs.rand = numeric(1000)
for (i in 1:length(process.accuracy.model.prob.rand)) {
  nrow.probs = nrow(df.probs)
  rnd.choices = sample(length(models),nrow.probs,replace=T)
  test = numeric(nrow.probs)
  for (j in 1:nrow.probs) {
    test[j] = as.numeric(df.probs[j,model.prob.names[rnd.choices[j]]] / max(df.probs[j,model.prob.names[1:6]]))
  }
  process.accuracy.model.prob.rand[i] = mean(test)
}

ggplot(df.demo.filt, aes(x = process.accuracy.model.relprob)) +
  geom_histogram(color = 'white', bins = 25) +
  labs(x = "\nReported model BF",
       y = "# of subjects") +
  scale_y_continuous(breaks = NULL) +
  theme_black() +
  geom_vline(xintercept = mean(df.demo.filt$process.accuracy.model.relprob), linetype = 1, color = 'red') +
  geom_vline(xintercept = mean(df.demo.filt$process.accuracy.model.relprob)+se(df.demo.filt$process.accuracy.model.relprob), linetype = 'dashed', color = 'red') +
  geom_vline(xintercept = mean(df.demo.filt$process.accuracy.model.relprob)-se(df.demo.filt$process.accuracy.model.relprob), linetype = 'dashed', color = 'red') +
  geom_vline(xintercept = mean(process.accuracy.model.prob.rand), color = 'gray')
t.test(df.demo.filt$process.accuracy.model.relprob - mean(process.accuracy.model.prob.rand)) 

# heuristic error
process.accuracy.div.rand = numeric(nrow(df.demo.filt))
process.accuracy.div.h1.rand = numeric(nrow(df.demo.filt))
process.accuracy.div.h2.rand = numeric(nrow(df.demo.filt))
process.accuracy.div.h3.rand = numeric(nrow(df.demo.filt))
for (subj in 1:nrow(df.demo.filt)) {
  q.probs = c(df.demo.filt$actual.h1.prob[subj],
              df.demo.filt$actual.h2.prob[subj],
              df.demo.filt$actual.h3.prob[subj])
    
  avg.div = numeric(3)
  for (q in 1:3) {
    avg.div[q] = mean(abs(runif(10000) - q.probs[q]))
  }
  process.accuracy.div.h1.rand[subj] = avg.div[1]
  process.accuracy.div.h2.rand[subj] = avg.div[2]
  process.accuracy.div.h3.rand[subj] = avg.div[3]
  process.accuracy.div.rand[subj] = mean(avg.div)
}
ggplot(df.demo.filt, aes(x = process.accuracy.div)) +
  geom_histogram(color = 'white') +
  geom_vline(xintercept = mean(df.demo.filt$process.accuracy.div), color = 'red') +
  geom_vline(xintercept = mean(df.demo.filt$process.accuracy.div) + se(df.demo.filt$process.accuracy.div), color = 'red', linetype = 'dashed') +
  geom_vline(xintercept = mean(df.demo.filt$process.accuracy.div) - se(df.demo.filt$process.accuracy.div), color = 'red', linetype = 'dashed') +
  geom_vline(xintercept = mean(process.accuracy.div.rand), color = 'gray') +
  labs(x = '\nHeuristic error', y = '# of subjects') +
  theme_black()
get.ci(df.demo.filt$process.accuracy.div)
t.test(df.demo.filt$process.accuracy.div - mean(process.accuracy.div.rand))

ggplot(df.demo.filt, aes(x = process.accuracy.div.h1)) +
  geom_histogram(color = 'white') +
  geom_vline(xintercept = mean(df.demo.filt$process.accuracy.div.h1), color = 'red') +
  geom_vline(xintercept = mean(df.demo.filt$process.accuracy.div.h1) + se(df.demo.filt$process.accuracy.div.h1), color = 'red', linetype = 'dashed') +
  geom_vline(xintercept = mean(df.demo.filt$process.accuracy.div.h1) - se(df.demo.filt$process.accuracy.div.h1), color = 'red', linetype = 'dashed') +
  geom_vline(xintercept = mean(process.accuracy.div.h1.rand), color = 'gray') +
  labs(x = '\nHeuristic error', y = '# of subjects') +
  theme_black()
t.test(df.demo.filt$process.accuracy.div.h1 - mean(process.accuracy.div.h1.rand))
ggplot(df.demo.filt, aes(x = process.accuracy.div.h2)) +
  geom_histogram(color = 'white') +
  geom_vline(xintercept = mean(df.demo.filt$process.accuracy.div.h2), color = 'red') +
  geom_vline(xintercept = mean(df.demo.filt$process.accuracy.div.h2) + se(df.demo.filt$process.accuracy.div.h2), color = 'red', linetype = 'dashed') +
  geom_vline(xintercept = mean(df.demo.filt$process.accuracy.div.h2) - se(df.demo.filt$process.accuracy.div.h2), color = 'red', linetype = 'dashed') +
  geom_vline(xintercept = mean(process.accuracy.div.h2.rand), color = 'gray') +
  labs(x = '\nHeuristic error', y = '# of subjects') +
  theme_black()
t.test(df.demo.filt$process.accuracy.div.h2 - mean(process.accuracy.div.h2.rand))
ggplot(df.demo.filt, aes(x = process.accuracy.div.h3)) +
  geom_histogram(color = 'white') +
  geom_vline(xintercept = mean(df.demo.filt$process.accuracy.div.h3), color = 'red') +
  geom_vline(xintercept = mean(df.demo.filt$process.accuracy.div.h3) + se(df.demo.filt$process.accuracy.div.h3), color = 'red', linetype = 'dashed') +
  geom_vline(xintercept = mean(df.demo.filt$process.accuracy.div.h3) - se(df.demo.filt$process.accuracy.div.h3), color = 'red', linetype = 'dashed') +
  geom_vline(xintercept = mean(process.accuracy.div.h3.rand), color = 'gray') +
  labs(x = '\nHeuristic error', y = '# of subjects') +
  theme_black()
t.test(df.demo.filt$process.accuracy.div.h3 - mean(process.accuracy.div.h3.rand))

ggplot(df.demo.filt, aes(x = process.accuracy.div, y = process.accuracy.model.relprob)) +
  geom_point(color = 'white') +
  geom_smooth(method='lm', color = 'white') +
  theme_black()

## comparing actual to reported
# continuous version
ggplot(df.demo.filt, aes(x = actual.h1.prob, y = reported.h1)) +
  geom_point(color='gray') +
  geom_smooth(method='lm', color = 'white') +
  theme_black() +
  #labs(x='',y='')
  labs(x = '\nProbability that participant used\nsingle attribute heuristic',
       y = 'Degree to which participant\nreported using heuristic')
summary(lm(scale(reported.h1) ~ scale(actual.h1.prob), df.demo.filt))
df.demo.filt %$% cor.test(reported.h1, actual.h1.prob)
ggplot(df.demo.filt, aes(x = actual.h2.prob, y = reported.h2)) +
  geom_point(color='gray') +
  geom_smooth(method='lm', color = 'white') +
  theme_black() +
  #labs(x='',y='')
  labs(x = '\nProbability that participant\nused binary weights heuristic',
       y = 'Degree to which participant\nreported using heuristic')
summary(lm(scale(reported.h2) ~ scale(actual.h2.prob), df.demo.filt))
df.demo.filt %$% cor.test(reported.h2, actual.h2.prob, method = 'kendall')
ggplot(df.demo.filt, aes(x = actual.h3.prob, y = reported.h3)) +
  geom_point(color='gray') +
  geom_smooth(method='lm', color = 'white') +
  theme_black() +
  #labs(x='',y='')
  labs(x = '\nProbability that participant\nused binary atts heuristic',
       y = 'Degree to which participant\nreported using heuristic')
summary(lm(reported.h3 ~ actual.h3.prob, df.demo.filt))
cor.test(df.demo.filt$actual.h3.prob, df.demo.filt$reported.h3)
df.demo.filt %$% cor.test(reported.h3, actual.h3.prob, method = 'kendall')

# dichotomized version
h1.grouped = df.demo.filt %>% group_by(actual.h1.prob.dichotomized) %>%
  summarize(reported.h1.dichotomized.m = mean(reported.h1.dichotomized),
            reported.h1.dichotomized.se = se.prop(reported.h1.dichotomized))
ggplot(h1.grouped, aes(x = actual.h1.prob.dichotomized, y = reported.h1.dichotomized.m)) +
  geom_col(color = 'white') +
  geom_errorbar(aes(ymin = reported.h1.dichotomized.m - reported.h1.dichotomized.se,
                    ymax = reported.h1.dichotomized.m + reported.h1.dichotomized.se),
                width = .2, color = 'white') +
  theme_black() +
  labs(x = 'Used single att heuristic\n(prob > 0.5)',
       y = '% reporting\nusing heuristic')

h2.grouped = df.demo.filt %>% group_by(actual.h2.prob.dichotomized) %>%
  summarize(reported.h2.dichotomized.m = mean(reported.h2.dichotomized),
            reported.h2.dichotomized.se = se.prop(reported.h2.dichotomized))
ggplot(h2.grouped, aes(x = actual.h2.prob.dichotomized, y = reported.h2.dichotomized.m)) +
  geom_col(color = 'white') +
  geom_errorbar(aes(ymin = reported.h2.dichotomized.m - reported.h2.dichotomized.se,
                    ymax = reported.h2.dichotomized.m + reported.h2.dichotomized.se),
                width = .2, color = 'white') +
  theme_black() +
  labs(x = 'Used binary wts heuristic\n(prob > 0.5)',
       y = '% reporting\nusing heuristic')

h3.grouped = df.demo.filt %>% group_by(actual.h3.prob.dichotomized) %>%
  summarize(reported.h3.dichotomized.m = mean(reported.h3.dichotomized),
            reported.h3.dichotomized.se = se.prop(reported.h3.dichotomized))
ggplot(h3.grouped, aes(x = actual.h3.prob.dichotomized, y = reported.h3.dichotomized.m)) +
  geom_col(color = 'white') +
  geom_errorbar(aes(ymin = reported.h3.dichotomized.m - reported.h3.dichotomized.se,
                    ymax = reported.h3.dichotomized.m + reported.h3.dichotomized.se),
                width = .2, color = 'white') +
  theme_black() +
  labs(x = 'Used binary atts heuristic\n(prob > 0.5)',
       y = '% reporting\nusing heuristic')
summary(glm(reported.h1.dichotomized ~ actual.h1.prob.dichotomized, df.demo.filt, family = 'binomial'))
summary(lm(scale(reported.h1) ~ actual.h1.prob.dichotomized, df.demo.filt))
summary(glm(reported.h2.dichotomized ~ actual.h2.prob.dichotomized, df.demo.filt, family = 'binomial'))
summary(lm(scale(reported.h2) ~ actual.h2.prob.dichotomized, df.demo.filt))
summary(glm(reported.h3.dichotomized ~ actual.h3.prob.dichotomized, df.demo.filt, family = 'binomial'))
summary(lm(scale(reported.h3) ~ actual.h3.prob.dichotomized, df.demo.filt))

# combine all 3 together
df.questions1 = df.demo.filt %>%
  select(subject, reported.h1, reported.h2, reported.h3) %>%
  rename(q1 = reported.h1, q2 = reported.h2, q3 = reported.h3) %>%
  pivot_longer(!subject, names_to = 'question', values_to = 'reported')
df.questions2 = df.demo.filt %>%
  select(subject, actual.h1.prob, actual.h2.prob, actual.h3.prob) %>%
  pivot_longer(!subject, values_to = 'prob')
df.questions3 = df.demo.filt %>%
  select(subject, reported.h1.dichotomized, reported.h2.dichotomized, reported.h3.dichotomized) %>%
  pivot_longer(!subject, values_to = 'chosen')
df.questions4 = df.demo.filt %>%
  select(subject, actual.h1.prob.dichotomized, actual.h2.prob.dichotomized, actual.h3.prob.dichotomized) %>%
  pivot_longer(!subject, values_to = 'best.family')
df.questions5 = df.demo.filt %>%
  select(subject, process.accuracy.div.h1, process.accuracy.div.h2, process.accuracy.div.h3) %>%
  pivot_longer(!subject, values_to = 'process.accuracy.div')
df.questions = df.questions1
df.questions$prob = df.questions2$prob
df.questions$chosen = df.questions3$chosen
df.questions$best.family = df.questions4$best.family
df.questions$process.accuracy.div = df.questions5$process.accuracy.div
df.questions = df.questions %>%
  mutate(correct = chosen == best.family)

ggplot(df.questions, aes(x = prob, y = reported)) +
  geom_point(color = 'gray') +
  geom_smooth(method='lm', color = 'white') +
  theme_black() +
  labs(x = '\nProbability that participant\nused heuristic',
       y = 'Reported probability of\nusing heuristic')
m1 = lmer(reported ~ prob + (prob | subject), df.questions)
summary(rePCA(m1))
m2 = lmer(scale(reported) ~ scale(prob) + (1 | subject), df.questions)
summary(m2)

q.grouped = df.questions %>% group_by(best.family) %>%
  summarize(chosen.m = mean(chosen),
            chosen.se = se.prop(chosen),
            correct.m = mean(correct),
            correct.se = se.prop(correct),
            process.accuracy.div.m = mean(process.accuracy.div),
            process.accuracy.div.se = se(process.accuracy.div))
ggplot(q.grouped, aes(x = best.family, y = chosen.m)) +
  geom_col(color = 'white') +
  geom_errorbar(aes(ymin = chosen.m - chosen.se,
                    ymax = chosen.m + chosen.se),
                width = .2, color = 'white') +
  theme_black() +
  labs(x = 'Used heuristic',
       y = 'Reported using\nheuristic') 
m1 = glmer(chosen ~ best.family + (best.family | subject), df.questions, family = 'binomial')
summary(rePCA(m1))
m2 = glmer(chosen ~ best.family + (1 | subject), df.questions, family = 'binomial')
summary(m2)

ggplot(q.grouped, aes(x = best.family, y = correct.m)) +
  geom_col(color = 'white') +
  geom_errorbar(aes(ymin = correct.m - correct.se,
                    ymax = correct.m + correct.se),
                width = .2, color = 'white') +
  theme_black() +
  labs(x = 'Used heuristic',
       y = '% reporting correctly') 
ggplot(q.grouped, aes(x = best.family, y = process.accuracy.div.m)) +
  geom_col(color = 'white') +
  geom_errorbar(aes(ymin = process.accuracy.div.m - process.accuracy.div.se,
                    ymax = process.accuracy.div.m + process.accuracy.div.se),
                width = .2, color = 'white') +
  theme_black() +
  labs(x = 'Used heuristic',
       y = 'Heuristic error') 
m1 = glmer(correct ~ best.family + (best.family | subject), df.questions, family = 'binomial')
summary(rePCA(m1))
m2 = glmer(correct ~ best.family + (1 | subject), df.questions, family = 'binomial')
summary(m2)

## % chosen correctly
features.graph = df.demo.filt %>% select(subject,process.accuracy.h1.dichotomized, process.accuracy.h2.dichotomized, process.accuracy.h3.dichotomized, process.accuracy.model.dichotomized) %>%
  pivot_longer(!subject) %>%
  group_by(name) %>%
  summarize(val = mean(value,na.rm=T),
            val.se = se(value)) %>%
  mutate(name = factor(name, c('process.accuracy.h1.dichotomized', 'process.accuracy.h2.dichotomized', 'process.accuracy.h3.dichotomized', 'process.accuracy.model.dichotomized')))
ggplot(features.graph, aes(x = name, y = val)) +
  geom_col(color='white') +
  geom_errorbar(aes(ymin = val - val.se, ymax = val + val.se), width = .2, color = 'white') +
  geom_segment(x = 0, y = .5, xend = 3.5, yend = 0.5, color = 'white', linetype = 'dashed') +
  geom_segment(x = 3.5, y = 1/6, xend = 4.5, yend = 1/6, color = 'white', linetype = 'dashed') +
  theme_black() +
  labs(x = '', y = '% correct') +
  scale_x_discrete(labels = c('Q1\n(one vs.\nmultiple attributes)', 'Q2\n(binary vs.\ngraded weights)', 'Q3\n(binary vs.\ngraded attributes)', 'Overall model'))

get.ci.prop(df.demo.filt$process.accuracy.model.dichotomized)

get.ci.prop(df.demo.filt$process.accuracy.h1.dichotomized)
get.ci.prop(df.demo.filt$process.accuracy.h2.dichotomized)
get.ci.prop(df.demo.filt$process.accuracy.h3.dichotomized)

# by question
num.subj.filt = nrow(df.demo.filt)
h1.heatmap = df.demo.filt %>%
  group_by(actual.h1.prob.dichotomized, reported.h1.dichotomized) %>%
  summarize(num = n() / num.subj.filt * 100)
h2.heatmap = df.demo.filt %>%
  group_by(actual.h2.prob.dichotomized, reported.h2.dichotomized) %>%
  summarize(num = n() / num.subj.filt * 100)
h3.heatmap = df.demo.filt %>%
  group_by(actual.h3.prob.dichotomized, reported.h3.dichotomized) %>%
  summarize(num = n() / num.subj.filt * 100)

summary(lm(process.accuracy.div.h2 ~ process.accuracy.div.h1, df.demo.filt))
summary(lm(process.accuracy.div.h3 ~ process.accuracy.div.h1 + actual.h1.prob, df.demo.filt))
summary(lm(process.accuracy.div.h3 ~ process.accuracy.div.h2 + actual.h2.prob, df.demo.filt))

ggplot(h1.heatmap, aes(x = actual.h1.prob.dichotomized, y = reported.h1.dichotomized,
                       fill = num)) +
  geom_tile() +
  labs(x = 'Used heuristic', y = 'Reported using heuristic') +
  #scale_fill_brewer(palette = 'YlOrRd') +
  guides(fill = guide_colorbar(title = '% of subjects')) +
  geom_text(aes(label = paste0(round(num), '%')), color = "black", size = 4) +
  theme_black() +
  scale_fill_gradient(low = "white", high = "red")
ggplot(h2.heatmap, aes(x = actual.h2.prob.dichotomized, y = reported.h2.dichotomized,
                       fill = num)) +
  geom_tile() +
  labs(x = 'Used heuristic', y = 'Reported using heuristic') +
  #scale_fill_brewer(palette = 'YlOrRd') +
  guides(fill = guide_colorbar(title = '% of subjects')) +
  geom_text(aes(label = paste0(round(num), '%')), color = "black", size = 4) +
  theme_black() +
  scale_fill_gradient(low = "white", high = "red")
ggplot(h3.heatmap, aes(x = actual.h3.prob.dichotomized, y = reported.h3.dichotomized,
                       fill = num)) +
  geom_tile() +
  labs(x = 'Used heuristic', y = 'Reported using heuristic') +
  #scale_fill_brewer(palette = 'YlOrRd') +
  guides(fill = guide_colorbar(title = '% of subjects')) +
  geom_text(aes(label = paste0(round(num), '%')), color = "black", size = 4) +
  theme_black() +
  scale_fill_gradient(low = "white", high = "red")

## look at stuff by model
summary(lm(actual.model.loo.readable ~ actual.model.fac, data = df.demo.filt))
summary(lm(consistency1 ~ actual.model.fac, df.demo.filt)) 
summary(lm(consistency2 ~ actual.model.fac, df.demo.filt)) 

df.demo.filt.model = df.demo.filt %>% group_by(actual.model.fac) %>%
  summarize(actual.model.loo.readable.m = mean(actual.model.loo.readable),
            actual.model.loo.readable.se = se(actual.model.loo.readable),
            consistency1.m = mean(consistency1),
            consistency1.se = se(consistency1),
            consistency2.m = mean(consistency2),
            consistency2.se = se(consistency2),
            weight.accuracy.averaged.m = mean(weight.accuracy.averaged),
            weight.accuracy.averaged.se = se(weight.accuracy.averaged),
            process.accuracy.model.relprob.m = mean(process.accuracy.model.relprob),
            process.accuracy.model.relprob.se = se(process.accuracy.model.relprob),
            process.accuracy.model.dichotomized.m = mean(process.accuracy.model.dichotomized),
            process.accuracy.model.dichotomized.se = se.prop(process.accuracy.model.dichotomized))
ggplot(df.demo.filt.model, aes(x = actual.model.fac, y = actual.model.loo.readable.m)) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = actual.model.loo.readable.m - actual.model.loo.readable.se,
                    ymax = actual.model.loo.readable.m + actual.model.loo.readable.se))

# Weight awareness --------------------------------------------------------

# what weights did people report?
ggplot(df.attributes.filt, aes(x = reported.weight.graded.signed)) +
  geom_histogram(color='white') +
  theme_black() +
  labs(x = '\nReported attribute weight', y = '# of attributes')


# across all points
ggplot(df.attributes.filt, aes(x = fitted.weight.averaged, y = reported.weight.averaged.signed)) +
  geom_point(color = 'gray') +
  theme_black() +
  geom_smooth(method='lm',color = 'white')
labs(x = 'Self-reported weight', y = 'Fitted weight') +
  scale_x_continuous(breaks = c(-1, 0, 1))
m = lmer(scale(fitted.weight) ~ scale(rating.signed) + (scale(rating.signed) | attribute) + (0 + scale(rating.signed) | subject), data = df.attributes.filt)
summary(m)

sigmoid.model = nls(reported.weight.averaged.signed ~ 2/(1+exp(-b * fitted.weight.averaged))-1,
                    start = list(b = 10),
                    data = df.attributes.filt)
sigmoid.params = coef(sigmoid.model)
sigmoid = function(x) {
  2 / (1 + exp(-sigmoid.params[1] * x))-1
}
df.attributes.filt$sigmoid.output = sigmoid(df.attributes.filt$fitted.weight.averaged)
ggplot(df.attributes.filt, aes(x = fitted.weight.averaged, y = reported.weight.averaged.signed)) +
  geom_point(color = 'gray') +
  geom_line(aes(y = sigmoid.output), color = 'white') +
  theme_black() +
  #geom_smooth(method='lm',formula='y ~ sigmoid(x)',color = 'white') +
  labs(x = 'Fitted weights', y = 'Reported weights') +
  scale_x_continuous(breaks = c(-1, 0, 1)) +
  scale_y_continuous(breaks = c(-1, 0, 1))

test = df.attributes.filt %>% group_by(subject) %>% summarize(accuracy.sigmoid = cor(reported.weight.averaged.signed, sigmoid.output))

# plot subject-level accuracies
rand.errs = numeric(1000)
for (i in 1:length(rand.errs)) {
  test = numeric(nrow(df.demo.filt))
  for (j in 1:nrow(df.demo.filt)) {
    subj.rows = df.attributes$subject == df.demo.filt$subject[j]
    reported = sample(df.attributes$reported.weight.averaged.signed[subj.rows])
    fitted = df.attributes$fitted.weight.averaged[subj.rows]
    test[j] = cor(reported, fitted)
  }
  rand.errs[i] = mean(test, na.rm = T)
}

ggplot(df.demo.filt, aes(x = weight.accuracy.averaged)) +
  geom_histogram(color = 'gray') +
  #geom_vline(xintercept = mean(rand.errs), color = 'gray') +
  geom_vline(xintercept = mean(df.demo.filt$weight.accuracy.averaged, na.rm = T), color = 'red') +
  geom_vline(xintercept = mean(df.demo.filt$weight.accuracy.averaged, na.rm = T)+se(df.demo.filt$weight.accuracy.averaged), color = 'red', linetype = 'dashed') +
  geom_vline(xintercept = mean(df.demo.filt$weight.accuracy.averaged, na.rm = T)-se(df.demo.filt$weight.accuracy.averaged), color = 'red', linetype = 'dashed') +
  labs(x = 'Participant-level correlation between\nfitted and reported weights', y = '# of subjects') +
  #labs(x = '', y = '') +
  scale_y_continuous(breaks = NULL) +
  #scale_x_continuous(limits = c(-.2, 1.1), breaks = c(0, .5, 1)) +
  theme_black()

get.ci(df.demo.filt$weight.accuracy.averaged)
get.ci(df.demo.filt$weight.accuracy.best)

ggplot(df.demo.filt.multiatt, aes(x = weight.accuracy.averaged)) +
  geom_histogram(color = 'gray') +
  geom_vline(xintercept = mean(rand.errs), color = 'gray') +
  geom_vline(xintercept = mean(df.demo.filt$weight.accuracy.averaged, na.rm = T), color = 'red') +
  geom_vline(xintercept = mean(df.demo.filt$weight.accuracy.averaged, na.rm = T)+se(df.demo.filt$weight.accuracy.averaged), color = 'red', linetype = 'dashed') +
  geom_vline(xintercept = mean(df.demo.filt$weight.accuracy.averaged, na.rm = T)-se(df.demo.filt$weight.accuracy.averaged), color = 'red', linetype = 'dashed') +
  labs(x = 'Participant-level correlation between\nfitted and reported weights', y = '# of subjects') +
  #labs(x = '', y = '') +
  scale_y_continuous(breaks = NULL) +
  #scale_x_continuous(limits = c(-.2, 1.1), breaks = c(0, .5, 1)) +
  theme_black()
get.ci(df.demo.filt.multiatt$weight.accuracy.averaged)

# relation to model fit
ggplot(df.demo.filt, aes(x = actual.model.r2, y = process.accuracy.div)) +
  geom_point(color = 'gray') +
  theme_black() +
  geom_smooth(method='lm', color = 'white') +
  labs(x = 'Model fit', y = 'Heuristic error') +
  scale_x_continuous(breaks = c(-1, 0, 1))
ggplot(df.demo.filt, aes(x = actual.model.r2, y = weight.accuracy.averaged)) +
  geom_point(color = 'gray') +
  theme_black() +
  geom_smooth(method='lm', color = 'white') +
  labs(x = 'Model fit', y = 'Weight accuracy') +
  scale_x_continuous(breaks = c(-1, 0, 1))

test = df.demo.filt %>% group_by(actual.model.fac) %>%
  summarize(process.accuracy.model.relprob.m = mean(process.accuracy.model.relprob),
            process.accuracy.model.relprob.se = se(process.accuracy.model.relprob),
            process.accuracy.div.m = mean(process.accuracy.div),
            process.accuracy.div.se = se(process.accuracy.div),
            weight.accuracy.averaged.m = mean(weight.accuracy.averaged),
            weight.accuracy.averaged.se = se(weight.accuracy.averaged))
ggplot(test, aes(x = actual.model.fac, y = process.accuracy.model.relprob.m)) +
  geom_col(color = 'gray') +
  geom_errorbar(aes(ymin = process.accuracy.model.relprob.m - process.accuracy.model.relprob.se,
                    ymax = process.accuracy.model.relprob.m + process.accuracy.model.relprob.se), color = 'white', width = .2) +
  theme_black() +
  labs(x = 'Best-fitting model', y = 'Reported model BF') +
  scale_x_discrete(labels = c('Full', 'H3', 'H2', 'H2&3', 'H1', 'H1&3'))
ggplot(test, aes(x = actual.model.fac, y = weight.accuracy.averaged.m)) +
  geom_col(color = 'gray') +
  geom_errorbar(aes(ymin = weight.accuracy.averaged.m - weight.accuracy.averaged.se,
                    ymax = weight.accuracy.averaged.m + weight.accuracy.averaged.se), color = 'white', width = .2) +
  theme_black() +
  labs(x = 'Best-fitting model', y = 'Weight accuracy') +
  scale_x_discrete(labels = c('Full', 'H3', 'H2', 'H2&3', 'H1', 'H1&3'))

summary(lm(process.accuracy.div ~ actual.model.r2 + actual.model.fac + cc.correct.total, df.demo.filt))
summary(lm(weight.accuracy.averaged ~ actual.model.r2 + actual.model.fac + cc.correct.total, df.demo.filt))

test = df.demo.filt %>% group_by(actual.h1.prob.dichotomized2) %>%
  summarize(process.accuracy.div.m = mean(process.accuracy.div.h1),
            process.accuracy.div.se = se(process.accuracy.div.h1),
            process.accuracy.h1.dichotomized.m = mean(process.accuracy.h1.dichotomized),
            process.accuracy.h1.dichotomized.se = se.prop(process.accuracy.h1.dichotomized))
ggplot(test, aes(x = actual.h1.prob.dichotomized2, y = process.accuracy.div.m)) +
  geom_col(color = 'white') +
  geom_errorbar(aes(ymin = process.accuracy.div.m - process.accuracy.div.se,
                    ymax = process.accuracy.div.m + process.accuracy.div.se), color = 'white', width = .2) +
  theme_black() +
  labs(x = 'Used heuristic?', y = 'Heuristic error')
ggplot(test, aes(x = actual.h1.prob.dichotomized2, y = process.accuracy.h1.dichotomized.m)) +
  geom_col(color = 'white') +
  geom_errorbar(aes(ymin = process.accuracy.h1.dichotomized.m - process.accuracy.h1.dichotomized.se,
                    ymax = process.accuracy.h1.dichotomized.m + process.accuracy.h1.dichotomized.se), color = 'white', width = .2) +
  theme_black() +
  labs(x = 'Used heuristic', y = '% who got it right')

test2 = df.demo.filt %>% group_by(actual.h2.prob.dichotomized2) %>%
  summarize(process.accuracy.div.m = mean(process.accuracy.div.h2),
            process.accuracy.div.se = se(process.accuracy.div.h2),
            process.accuracy.h2.dichotomized.m = mean(process.accuracy.h2.dichotomized),
            process.accuracy.h2.dichotomized.se = se.prop(process.accuracy.h2.dichotomized))
ggplot(test2, aes(x = actual.h2.prob.dichotomized2, y = process.accuracy.div.m)) +
  geom_col(color = 'white') +
  geom_errorbar(aes(ymin = process.accuracy.div.m - process.accuracy.div.se,
                    ymax = process.accuracy.div.m + process.accuracy.div.se), color = 'white', width = .2) +
  theme_black() +
  labs(x = 'Used heuristic?', y = 'Heuristic error')
ggplot(test2, aes(x = actual.h2.prob.dichotomized2, y = process.accuracy.h2.dichotomized.m)) +
  geom_col(color = 'white') +
  geom_errorbar(aes(ymin = process.accuracy.h2.dichotomized.m - process.accuracy.h2.dichotomized.se,
                    ymax = process.accuracy.h2.dichotomized.m + process.accuracy.h2.dichotomized.se), color = 'white', width = .2) +
  theme_black() +
  labs(x = 'Used heuristic', y = '% who got it right')

test3 = df.demo.filt %>% group_by(actual.h3.prob.dichotomized2) %>%
  summarize(process.accuracy.div.m = mean(process.accuracy.div.h3),
            process.accuracy.div.se = se(process.accuracy.div.h3),
            process.accuracy.h3.dichotomized.m = mean(process.accuracy.h3.dichotomized),
            process.accuracy.h3.dichotomized.se = se.prop(process.accuracy.h3.dichotomized))
ggplot(test3, aes(x = actual.h3.prob.dichotomized2, y = process.accuracy.div.m)) +
  geom_col(color = 'white') +
  geom_errorbar(aes(ymin = process.accuracy.div.m - process.accuracy.div.se,
                    ymax = process.accuracy.div.m + process.accuracy.div.se), color = 'white', width = .2) +
  theme_black() +
  labs(x = 'Used heuristic?', y = 'Heuristic error')
ggplot(test3, aes(x = actual.h3.prob.dichotomized2, y = process.accuracy.h3.dichotomized.m)) +
  geom_col(color = 'white') +
  geom_errorbar(aes(ymin = process.accuracy.h3.dichotomized.m - process.accuracy.h3.dichotomized.se,
                    ymax = process.accuracy.h3.dichotomized.m + process.accuracy.h3.dichotomized.se), color = 'white', width = .2) +
  theme_black() +
  labs(x = 'Used heuristic', y = '% who got it right')
summary()

# Compare heuristic and weight awareness -------------------------------

df.demo.correct = df.demo.filt %>% group_by(process.accuracy.model.dichotomized) %>%
  summarize(accuracy.m = mean(weight.accuracy.averaged, na.rm = T), accuracy.se = se(weight.accuracy.averaged),
            accuracy.rank.m = mean(weight.accuracy.averaged.rank, na.rm = T), accuracy.rank.se = se(weight.accuracy.averaged.rank),
            accuracy.reported.m = mean(weight.accuracy.reported, na.rm = T), accuracy.reported.se = se(weight.accuracy.reported))
# averaged weights
ggplot(df.demo.filt, aes(x = process.accuracy.model.dichotomized, y = weight.accuracy.averaged.rank)) +
  geom_violin() +
  geom_point(data = df.demo.correct, aes(y = accuracy.rank.m)) +
  geom_errorbar(data = df.demo.correct,
                aes(y = accuracy.m, ymin = accuracy.m - accuracy.se,
                    ymax = accuracy.m + accuracy.se),
                width = .1) +
  theme_black() +
  labs(x = 'Chose correct model?', y = 'Weight accuracy\n(with model-averaged wts)') +
  scale_y_continuous(breaks = c(0,.5,1))
ggplot(df.demo.filt, aes(x = process.accuracy.div, y = weight.accuracy.averaged)) +
  geom_point(color = 'gray') +
  geom_smooth(method = 'lm', color = 'white') +
  theme_black() +
  labs(x = 'Reported model BF', y = 'Weight accuracy\n(with model-averaged wts)')
summary(lm(scale(weight.accuracy.averaged) ~ scale(process.accuracy.model.relprob) + actual.model.fac + actual.model.loo.readable, df.demo.filt))

# reported weights
ggplot(df.demo.filt, aes(x = process.accuracy.model.dichotomized, y = weight.accuracy.reported)) +
  geom_violin() +
  geom_point(data = df.demo.correct, aes(y = accuracy.reported.m)) +
  geom_errorbar(data = df.demo.correct,
                aes(y = accuracy.reported.m, ymin = accuracy.reported.m - accuracy.reported.se,
                    ymax = accuracy.reported.m + accuracy.reported.se),
                width = .1) +
  theme_black() +
  labs(x = 'Chose correct model?', y = 'Weight accuracy for reported model') +
  scale_y_continuous(breaks = c(0,.5,1))
ggplot(df.demo.filt, aes(x = process.accuracy.model.relprob, y = weight.accuracy.reported)) +
  geom_point(color = 'gray') +
  geom_smooth(method = 'lm', color =' white') +
  theme_black() +
  labs(x = 'Reported model BF', y = 'Weight accuracy\n(with weights of reported model)')
summary(lm(scale(weight.accuracy.reported) ~ scale(process.accuracy.model.relprob) + actual.model.fac + scale(actual.model.r2), df.demo.filt))

# kld / CEL

plot(df.demo.filt$process.accuracy.div, df.demo.filt$weight.accuracy.averaged)
cor.test(df.demo.filt$process.accuracy.div, df.demo.filt$weight.accuracy.averaged)
summary(lm(scale(weight.accuracy.averaged.rank) ~ scale(process.accuracy.div) + actual.model.fac + scale(actual.model.r2), df.demo.filt))

# Analyze potential correlates of awareness --------------------------------------------------------------

# comprehension checks
ggplot(df.demo.filt, aes(x = cc.correct.total, y = process.accuracy.model.prob)) +
  geom_point(color='gray') +
  geom_smooth(method='lm',color='gray') +
  theme_black() +
  labs(x = '\nTotal correct CC questions', y = 'Scaled likelihood of chosen model')
ggplot(df.demo.filt, aes(x = cc.correct.total, y = num.qs.correct)) +
  geom_point(color='gray') +
  geom_smooth(method='lm',color='gray') +
  theme_black() +
  labs(x = '\nTotal correct CC questions', y = 'Scaled likelihood of chosen model')
summary(lm(reported.model.ll ~ cc.correct.total, df.demo.filt))
summary(glm(chose.process.accuracy.model.dichotomized ~ cc.correct.total, df.demo.filt, family = 'binomial'))

# experience w/ choice domain
ggplot(df.demo.filt, aes(x = choice_domain, y = process.accuracy.div)) +
  geom_point(color='gray') +
  geom_smooth(method='lm', color = 'white') +
  theme_black() +
  labs(x = '\nExperience with choosing\nhomes to rent', y = 'Weight accuracy') +
  scale_x_continuous(breaks = c(0,50,100)) +
  scale_y_continuous(breaks = c(0,.5,1))
ggplot(df.demo.filt, aes(x = choice_domain, y = weight.accuracy.averaged)) +
  geom_point(color='gray') +
  geom_smooth(method='lm', color = 'white') +
  theme_black() +
  labs(x = '\nExperience with choosing\nhomes to rent', y = 'Weight accuracy') +
  scale_x_continuous(breaks = c(0,50,100)) +
  scale_y_continuous(breaks = c(0,.5,1))
summary(lm(scale(weight.accuracy.averaged) ~ scale(choice_domain) +
             actual.model.fac + actual.model.r2, data = df.demo.filt))
summary(lm(scale(process.accuracy.model.prob) ~ scale(choice_domain) +
             actual.model.fac + actual.model.r2, data = df.demo.filt))

# experience w/ tasks like this
ggplot(df.demo.filt, aes(x = choice_exp_num, y = weight.accuracy.averaged)) +
  geom_point() +
  geom_smooth(method = 'lm')
ggplot(df.demo.filt, aes(x = choice_exp_num, y = process.accuracy.model.prob)) +
  geom_point() +
  geom_smooth(method='lm')
summary(lm(scale(weight.accuracy.averaged) ~ scale(choice_exp_num) +
             actual.model.fac + actual.model.r2, data = df.demo.filt))
summary(lm(scale(process.accuracy.model.prob) ~ scale(choice_exp_num) +
             actual.model.fac + actual.model.r2, data = df.demo.filt))

# confidence
ggplot(df.demo.filt, aes(x = confidence, y = weight.accuracy.averaged)) +
  geom_point() +
  geom_smooth(method='lm')
summary(lm(scale(weight.accuracy.averaged) ~ scale(confidence) +
             actual.model.fac + actual.model.r2, data = df.demo.filt))
ggplot(df.demo.filt, aes(x = confidence, y = process.accuracy.model.relprob)) +
  geom_point() +
  geom_smooth(method='lm')
summary(lm(scale(process.accuracy.model.prob) ~ scale(confidence) +
             actual.model.fac + actual.model.r2, data = df.demo.filt))

# decision style
ggplot(df.demo.filt, aes(x = decisionstyle, y = process.accuracy.model.prob)) +
  geom_point() +
  geom_smooth(method='lm')
ggplot(df.demo.filt, aes(x = decisionstyle, y = weight.accuracy.averaged)) +
  geom_point() +
  geom_smooth(method='lm')

# mindfulness
ggplot(df.demo.filt, aes(x = mindfulness, y = weight.accuracy.averaged)) +
  geom_point() +
  geom_smooth(method='lm', color = 'brown') +
  labs(x = 'CAMS-R score', y = 'Parameter accuracy') +
  scale_x_continuous(breaks = c(2,10)) +
  scale_y_continuous(breaks = c(0, .5, 1))
ggplot(df.demo.filt, aes(x = mindfulness, y = process.accuracy.model.prob)) +
  geom_point() +
  geom_smooth(method='lm', color = 'brown') +
  labs(x = 'CAMS-R score', y = 'Parameter accuracy') +
  scale_x_continuous(breaks = c(2,10)) +
  scale_y_continuous(breaks = c(0, .5, 1))

# meditation
meditation.graph = df.demo.filt %>% group_by(meditation_exp1) %>%
  summarize(weight.accuracy.averaged.m = mean(weight.accuracy.averaged, na.rm = T),
            weight.accuracy.averaged.se = se(weight.accuracy.averaged),
            process.accuracy.qs.m = mean(process.accuracy.qs),
            process.accuracy.qs.se = se(process.accuracy.qs),
            process.accuracy.model.relprob.m = mean(process.accuracy.model.relprob, na.rm = T),
            process.accuracy.model.relprob.se = se(process.accuracy.model.relprob),
            process.accuracy.model.dichotomized.m = mean(process.accuracy.model.dichotomized, na.rm = T),
            process.accuracy.model.dichotomized.se = se.prop(process.accuracy.model.dichotomized))
ggplot(meditation.graph, aes(x = meditation_exp1, y = process.accuracy.model.relprob.m)) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = process.accuracy.model.relprob.m - process.accuracy.model.relprob.se,
                    ymax = process.accuracy.model.relprob.m + process.accuracy.model.relprob.se,
                    width = .2))
ggplot(meditation.graph, aes(x = meditation_exp1, y = process.accuracy.model.dichotomized.m)) +
  geom_point(size = 5, color = 'white') +
  geom_errorbar(aes(ymin = process.accuracy.model.dichotomized.m - process.accuracy.model.dichotomized.se,
                    ymax = process.accuracy.model.dichotomized.m + process.accuracy.model.dichotomized.se),
                width = .2, color = 'white') +
  #scale_y_continuous(limits = c(0.1, 0.4), breaks = c(0.1, 0.2, 0.3, 0.4), labels = c('10', '20', '30', '40')) +
  labs(x = '\nPrior meditation experience', y = '% reporting correct model') +
  theme_black()
ggplot(meditation.graph, aes(x = meditation_exp1, y = weight.accuracy.averaged.m)) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = weight.accuracy.averaged.m - weight.accuracy.averaged.se,
                    ymax = weight.accuracy.averaged.m + weight.accuracy.averaged.se,
                    width = .2))
ggplot(meditation.graph, aes(x = meditation_exp1, y = num.features.correct.m)) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = num.features.correct.m - num.features.correct.se,
                    ymax = num.features.correct.m + num.features.correct.se,
                    width = .2))
summary(lm(scale(weight.accuracy.averaged) ~ meditation_exp1 +
             actual.model.fac + actual.model.r2, data = df.demo.filt))
summary(lm(scale(process.accuracy.model.prob) ~ meditation_exp1 +
             actual.model.fac + actual.model.r2, data = df.demo.filt))

meditation.graph2 = df.demo.filt %>% filter(meditation_exp1 == 'Yes') %>%
  group_by(meditation_exp2) %>%
  summarize(weight.accuracy.averaged.m = mean(weight.accuracy.averaged, na.rm = T),
            weight.accuracy.averaged.se = se(weight.accuracy.averaged),
            process.accuracy.model.relprob.m = mean(process.accuracy.model.relprob, na.rm = T),
            process.accuracy.model.relprob.se = se(process.accuracy.model.relprob),
            process.accuracy.model.dichotomized.m = mean(process.accuracy.model.dichotomized, na.rm = T),
            process.accuracy.model.dichotomized.se = se.prop(process.accuracy.model.dichotomized),
            process.accuracy.qs.m = mean(process.accuracy.qs),
            process.accuracy.qs.se = se(process.accuracy.qs),)
ggplot(meditation.graph2, aes(x = meditation_exp2, y = process.accuracy.model.relprob.m)) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = process.accuracy.model.relprob.m - process.accuracy.model.relprob.se,
                    ymax = process.accuracy.model.relprob.m + process.accuracy.model.relprob.se,
                    width = .2))
ggplot(meditation.graph2, aes(x = meditation_exp2, y = process.accuracy.model.dichotomized.m)) +
  geom_point(size = 5, color = 'white') +
  geom_errorbar(aes(ymin = process.accuracy.model.dichotomized.m - process.accuracy.model.dichotomized.se,
                    ymax = process.accuracy.model.dichotomized.m + process.accuracy.model.dichotomized.se),
                width = .2, color = 'white') +
  #scale_y_continuous(limits = c(0.1, 0.4), breaks = c(0.1, 0.2, 0.3, 0.4), labels = c('10', '20', '30', '40')) +
  labs(x = '\nPrior meditation experience', y = '% reporting correct model') +
  theme_black()
ggplot(meditation.graph2, aes(x = meditation_exp2, y = weight.accuracy.averaged.m)) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = weight.accuracy.averaged.m - weight.accuracy.averaged.se,
                    ymax = weight.accuracy.averaged.m + weight.accuracy.averaged.se,
                    width = .2))

meditation.graph3 = df.demo.filt %>%
  group_by(meditation.over5) %>%
  summarize(weight.accuracy.averaged.m = mean(weight.accuracy.averaged, na.rm = T),
            weight.accuracy.averaged.se = se(weight.accuracy.averaged),
            process.accuracy.model.relprob.m = mean(process.accuracy.model.relprob, na.rm = T),
            process.accuracy.model.relprob.se = se(process.accuracy.model.relprob),
            process.accuracy.model.dichotomized.m = mean(process.accuracy.model.dichotomized, na.rm = T),
            process.accuracy.model.dichotomized.se = se.prop(process.accuracy.model.dichotomized),
            process.accuracy.qs.m = mean(process.accuracy.qs),
            process.accuracy.qs.se = se(process.accuracy.qs))
ggplot(meditation.graph3, aes(x = meditation.over5, y = process.accuracy.model.relprob.m)) +
  geom_point(size = 5, color = 'white') +
  geom_errorbar(aes(ymin = process.accuracy.model.relprob.m - process.accuracy.model.relprob.se,
                    ymax = process.accuracy.model.relprob.m + process.accuracy.model.relprob.se,
                    width = .2), color = 'white') +
  labs(x = '\n>5 yrs meditation experience', y = 'Posterior prob. of reported model') +
  theme_black()
ggplot(meditation.graph3, aes(x = meditation.over5, y = process.accuracy.model.dichotomized.m)) +
  geom_point(size = 5, color = 'white') +
  geom_errorbar(aes(ymin = process.accuracy.model.dichotomized.m - process.accuracy.model.dichotomized.se,
                    ymax = process.accuracy.model.dichotomized.m + process.accuracy.model.dichotomized.se),
                width = .2, color = 'white') +
  #scale_y_continuous(limits = c(0.1, 0.4), breaks = c(0.1, 0.2, 0.3, 0.4), labels = c('10', '20', '30', '40')) +
  labs(x = '\n>5 yrs meditation experience', y = '% reporting correct model') +
  theme_black()
ggplot(meditation.graph3, aes(x = meditation.over5, y = weight.accuracy.averaged.m)) +
  geom_point(size = 5, color = 'white') +
  geom_errorbar(aes(ymin = weight.accuracy.averaged.m - weight.accuracy.averaged.se,
                    ymax = weight.accuracy.averaged.m + weight.accuracy.averaged.se),
                width = .2, color = 'white') +
  labs(x = '\n>5 yrs meditation experience', y = 'Weight accuracy') +
  theme_black()

summary(lm(scale(weight.accuracy.averaged) ~ meditation.over5 +
             actual.model.fac + scale(actual.model.r2), data = df.demo.filt))
summary(lm(scale(process.accuracy.model.relprob) ~ meditation.over5 +
             actual.model.fac + scale(actual.model.r2), data = df.demo.filt))
summary(glm(process.accuracy.model.dichotomized ~ meditation.over5 +
              actual.model.fac + actual.model.r2, data = df.demo.filt, family = 'binomial'))

# norms
ggplot(df.demo.filt, aes(x = process.accuracy.model.relprob, y = norm.actual)) +
  geom_point(size = 1, color = 'gray') +
  geom_smooth(method='lm', color = 'white') +
  labs(x = 'Reported model BF', y = 'Prob. that they used\naligned model') +
  theme_black()
ggplot(df.demo.filt, aes(x = weight.accuracy.averaged, y = norm.actual)) +
  geom_point(size = 1, color = 'white') +
  labs(x = 'Reported correct\nmodel', y = 'Normative alignment of process') +
  theme_black()
norms.grouped = df.demo.filt %>% group_by(process.accuracy.model.dichotomized) %>%
  summarize(norm.actual.m = mean(norm.actual),
            norm.actual.se = se(norm.actual))
ggplot(norms.grouped, aes(x = process.accuracy.model.dichotomized, y = norm.actual.m)) +
  geom_point(size = 5, color = 'white') +
  geom_errorbar(aes(ymin = norm.actual.m - norm.actual.se,
                    ymax = norm.actual.m + norm.actual.se),
                width = .2, color = 'white') +
  scale_x_discrete(labels = c('No', 'Yes')) +
  labs(x = 'Reported correct\nmodel', y = 'Prob. that they used\naligned model') +
  theme_black()
summary(lm(scale(norm.actual) ~ process.accuracy.model.dichotomized + actual.model.fac + actual.model.r2, data = df.demo.filt))
summary(lm(scale(norm.actual) ~ scale(process.accuracy.model.relprob) + actual.model.fac + actual.model.r2, data = df.demo.filt))


test = df.demo.filt %>% group_by(h1.correct) %>%
  summarize(h1.norm.actual.m = mean(h1.norm.actual),
            h1.norm.actual.se = se.prop(h1.norm.actual))
ggplot(test, aes(x = h1.correct, y = h1.norm.actual.m)) +
  geom_point(size = 5, color = 'white') +
  geom_errorbar(aes(ymin = h1.norm.actual.m - h1.norm.actual.se,
                    ymax = h1.norm.actual.m + h1.norm.actual.se),
                width = .2, color = 'white') +
  scale_x_discrete(labels = c('No', 'Yes')) +
  labs(x = 'Reported correct\nmodel', y = '% who actually did') +
  theme_black()

test = df.demo.filt %>% group_by(h1.correct, norm.h1) %>%
  summarize(h1.norm.actual.m = mean(h1.norm.actual),
            h1.norm.actual.se = se.prop(h1.norm.actual))
ggplot(test, aes(x = h1.correct, y = h1.norm.actual.m, color = norm.h1, group = norm.h1)) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = h1.norm.actual.m - h1.norm.actual.se,
                    ymax = h1.norm.actual.m + h1.norm.actual.se),
                width = .2, color = 'white') +
  scale_x_discrete(labels = c('No', 'Yes')) +
  labs(x = 'Reported correct\nmodel', y = '% who actually did') +
  theme_black()

test = df.demo.filt.obs %>% group_by(chose.process.accuracy.model.dichotomized) %>%
  summarize(pct.features.norm.m = mean(pct.features.norm),
            pct.features.norm.se = se(pct.features.norm),
            pct.features.norm.reported.m = mean(pct.features.norm.reported),
            pct.features.norm.reported.se = se(pct.features.norm.reported),)
ggplot(test, aes(x = chose.process.accuracy.model.dichotomized, y = pct.features.norm.m)) +
  geom_point(size = 5, color = 'white') +
  geom_errorbar(aes(ymin = pct.features.norm.m - pct.features.norm.se,
                    ymax = pct.features.norm.m + pct.features.norm.se),
                width = .2, color = 'white') +
  scale_x_discrete(labels = c('No', 'Yes')) +
  labs(x = 'Reported correct\nmodel', y = '% of strategy properties\nthat were\nnormatively aligned') +
  theme_black()

#IQ
ggplot(df.demo.filt, aes(x = icar_num_correct, y = weight.accuracy.averaged)) +
  geom_point(color = 'gray') +
  geom_smooth(method='lm', color = 'white') +
  theme_black() +
  labs(x = 'IQ score', y = 'Weight accuracy')
ggplot(df.demo.filt, aes(x = icar_num_correct, y = process.accuracy.div)) +
  geom_jitter(color = 'gray') +
  geom_smooth(method='lm', color = 'white') +
  theme_black() +
  labs(x = 'IQ score', y = 'Reported model BF')
ggplot(df.demo.filt, aes(x = icar_num_correct, y = process.accuracy.qs)) +
  geom_point(color = 'gray') +
  geom_smooth(method='lm', color = 'white') +
  theme_black() +
  labs(x = 'IQ score', y = 'Process question accuracy')
summary(lm(scale(weight.accuracy.averaged) ~ scale(icar_num_correct) + actual.model.fac + actual.model.r2, df.demo.filt))
summary(lm(scale(process.accuracy.model.relprob) ~ scale(icar_num_correct) + actual.model.fac + scale(actual.model.r2), df.demo.filt))
summary(lm(scale(process.accuracy.model.relprob) ~ scale(icar_num_correct) + actual.model.fac + scale(actual.model.r2), df.demo.filt))

summary(glm(process.accuracy.model.dichotomized ~ scale(icar_num_correct) + actual.model.fac + actual.model.r2, df.demo.filt, family = 'binomial'))

# education
ggplot(df.demo.filt, aes(x = edu.num, y = weight.accuracy.averaged)) +
  geom_point(color = 'gray') +
  geom_smooth(method='lm', color = 'white') +
  theme_black() +
  labs(x = 'IQ score', y = 'Parameter accuracy')

# maia
ggplot(df.demo.filt, aes(x = maia, y = weight.accuracy.averaged)) +
  geom_point() +
  geom_smooth(method='lm')
ggplot(df.demo.filt, aes(x = maia, y = process.accuracy.model.prob)) +
  geom_point() +
  geom_smooth(method='lm')
summary(lm(weight.accuracy.best ~ maia + reported.model.fac + avg.model.ll, df.demo.filt))
summary(lm(reported.model.ll ~ maia + reported.model.fac + avg.model.ll, df.demo.filt))

# sris
ggplot(df.demo.filt, aes(x = sris.insight, y = weight.accuracy.averaged)) +
  geom_point() +
  geom_smooth(method='lm')
ggplot(df.demo.filt, aes(x = sris.insight, y = process.accuracy.model.prob)) +
  geom_point() +
  geom_smooth(method='lm')
summary(lm(weight.accuracy.best ~ sris.insight + best.model.fac + avg.model.ll, df.demo.filt))
summary(lm(reported.model.ll ~ sris.insight + best.model.fac + avg.model.ll, df.demo.filt))

ggplot(df.demo.filt, aes(x = sris.tendency, y = weight.accuracy.averaged)) +
  geom_point() +
  geom_smooth(method='lm')
ggplot(df.demo.filt, aes(x = sris.tendency, y = process.accuracy.model.prob)) +
  geom_point() +
  geom_smooth(method='lm')

# acs
ggplot(df.demo.filt, aes(x = acs.shifting, y = process.accuracy.div)) +
  geom_point(color = 'gray') +
  geom_smooth(method='lm', color = 'gray') +
  labs(x = "\nAttentional Control Scale\nscore", y = "Mean squared error") +
  theme_black()
#scale_x_continuous(breaks = c(20,100), limits = c(20,100))
#scale_y_continuous(breaks = c(0,1), limits = c(0,1))

# ANT stuff
ggplot(df.demo.filt, aes(x = orienting, y = weight.accuracy.averaged)) +
  geom_point() +
  geom_smooth(method='lm')
ggplot(df.demo.filt, aes(x = orienting, y = process.accuracy.model.prob)) +
  geom_point() +
  geom_smooth(method='lm')
summary(lm(scale(weight.accuracy.averaged) ~ scale(orienting) + actual.model.fac + actual.model.r2, df.demo.filt))
summary(lm(scale(process.accuracy.model.prob) ~ scale(orienting) + actual.model.fac + actual.model.r2, df.demo.filt))
summary(glm(process.accuracy.model.dichotomized ~ scale(orienting) + actual.model.fac + actual.model.r2, df.demo.filt, family = 'binomial'))

ggplot(df.demo.filt, aes(x = alerting, y = weight.accuracy.averaged)) +
  geom_point() +
  geom_smooth(method='lm')
ggplot(df.demo.filt, aes(x = alerting, y = process.accuracy.model.prob)) +
  geom_point() +
  geom_smooth(method='lm')

ggplot(df.demo.filt, aes(x = exec, y = weight.accuracy.averaged)) +
  geom_point() +
  geom_smooth(method='lm')
ggplot(df.demo.filt, aes(x = exec, y = process.accuracy.model.prob)) +
  geom_point() +
  geom_smooth(method='lm')

# satisfaction

ggplot(df.demo.filt, aes(x = weight.accuracy.averaged, y = satisfaction)) +
  geom_point() +
  geom_smooth(method='lm')
ggplot(df.demo.filt, aes(x = process.accuracy.model.prob, y = satisfaction)) +
  geom_point() +
  geom_smooth(method='lm')
ggplot(df.demo.filt %>% group_by(chose.process.accuracy.model.dichotomized) %>% summarize(s = mean(satisfaction)),
       aes(x = chose.process.accuracy.model.dichotomized, y = s)) +
  geom_col()

# power analysis
p_vals = numeric(100)
df.pwr = df.demo
for (i in 1:100) {
  print(i)
  resampled <-
    df.pwr %>%
    distinct(subject) %>%
    slice_sample(n = 250, replace = T) %>%
    group_by(subject) %>%
    mutate(instance = row_number()) %>%
    ungroup() %>%
    left_join(df.pwr,
              by = "subject") %>%
    mutate(subject = str_c(subject, instance))
  
  resampled_model <- glm(chose.process.accuracy.model.dichotomized ~ acs, family = 'binomial',
                         data = resampled)
  p_vals[i] <- summary(resampled_model)$coefficients[2,4]
}

m.mods = lm(weight.accuracy.averaged ~ decisionstyle + mindfulness + sris.tendency + sris.insight + maia + acs.shifting + acs.focusing +
              icar_num_correct + gender + age + edu.num + meditation_exp1 +
              alerting + orienting + exec +
              choice_domain + confidence + consistency1 + consistency2 + choice_exp_num +
              actual.model.fac + actual.model.r2,
            data = df.demo.filt)
summary(m.mods)

m.mods2 = lm(process.accuracy.div ~ decisionstyle + mindfulness + sris.tendency + sris.insight + maia + acs.shifting + acs.focusing +
               icar_num_correct + gender + age + edu.num + meditation_exp1 +
               alerting + orienting + exec +
               choice_domain + confidence + consistency1 + consistency2 + choice_exp_num +
               actual.model.fac + actual.model.r2,
             data = df.demo.filt)
summary(m.mods2)

# Save observer data for comparison with deciders ------------------------------------------------------
# Save observer data to the observer folder
if (participant_type == 'observer') {
  df.demo.filt.obs = df.demo.filt
  save(df.demo.filt.obs, file = paste0(filepath, 'observer_results.rdata'))
}

# split-half stuff
df.demo.filt.odd = df.demo.filt
df.attributes.filt.odd = df.attributes.filt
save(df.demo.filt.odd, df.attributes.filt.odd, file = paste0(filepath_modeloutput, '/dfs.rdata'))

# Compare deciders to observers ----------------------------------------------------
# This section should only be run if participant_type == 'decider'

load(paste0(filepath_observer, 'observer_results.rdata'))
df.demo.filt$type = 'Original'
df.demo.filt.obs$type = 'Observers'

df.demo.filt.both = full_join(df.demo.filt, df.demo.filt.obs)
df.demo.filt.both = df.demo.filt.both %>% mutate(type = factor(type, c('Original', 'Observers'), c('Original', 'Observers')))

df.demo.filt.both.multiatt = df.demo.filt.both %>%
  filter(actual.h1.prob.dichotomized == 'Multiple')

### process awareness
## heat map
df.demo.heat.obs = df.demo.filt.obs %>% group_by(reported.model.fac, actual.model.fac) %>%
  summarize(num.subj = n())
ggplot(df.demo.heat.obs, aes(x = actual.model.fac, y = reported.model.fac,
                             fill = num.subj)) +
  geom_tile() +
  labs(y = '\nSelf-reported model', x = 'Best-fitting model') +
  #scale_fill_brewer(palette = 'YlOrRd') +
  guides(fill = guide_colorbar(title = '# of subjects')) +
  theme_black()
df.demo.heat.normed.obs = df.demo.filt.obs %>% group_by(reported.model.fac, actual.model.fac) %>%
  summarize(num.subj = n()) %>%
  group_by(actual.model.fac) %>%
  mutate(num.subj.norm = num.subj / sum(num.subj))
ggplot(df.demo.heat.normed.obs, aes(x = actual.model.fac, y = reported.model.fac,
                                    fill = num.subj.norm)) +
  geom_tile() +
  geom_text(aes(label = round(num.subj.norm, 2))) +
  labs(y = '\nSelf-reported model', x = 'Best-fitting model') +
  #scale_fill_brewer(palette = 'YlOrRd') +
  scale_fill_continuous(limits = c(0,1), low = 'black', high = 'white') +
  guides(fill = guide_colorbar(title = '% of subjects')) +
  #theme_black() +
  scale_x_discrete(labels = c('Rational', 'Binary att vals', 'Binary wts', 'Binary wts +\nbinary att vals', 'Single att', 'Single Att +\nBinary att vals')) +
  scale_y_discrete(labels = c('Rational', 'Binary att vals', 'Binary wts', 'Binary wts +\nbinary att vals', 'Single att', 'Single Att +\nBinary att vals')) +
  theme_grey(base_size = 12) %+replace%
  theme(
    # Specify axis options
    axis.text.x = element_text(size = 10, color = "white", angle = 45, margin = margin(0, 10, 0, 0)),  
    axis.text.y = element_text(size = 10, color = "white", margin = margin(0, 10, 0, 0)),  
    axis.line = element_blank(),  
    axis.ticks = element_line(color = "white", size  =  0.2),  
    axis.title.x = element_text(size = 18, color = "white", margin = margin(0, 10, 0, 0)),  
    axis.title.y = element_text(size = 18, color = "white", angle = 90, margin = margin(0, 10, 0, 0)),  
    axis.ticks.length = unit(0.3, "lines"),   
    # Specify legend options
    legend.background = element_rect(color = NA, fill = "black"),  
    legend.key = element_rect(color = "white",  fill = "black"),  
    legend.key.size = unit(1.2, "lines"),  
    legend.key.height = NULL,  
    legend.key.width = NULL,      
    legend.text = element_text(size = 12*0.8, color = "white"),  
    legend.title = element_text(size = 12*0.8, face = "bold", hjust = 0, color = "white"),  
    legend.position = "right",  
    legend.text.align = NULL,  
    legend.title.align = NULL,  
    legend.direction = "vertical",  
    legend.box = NULL, 
    # Specify panel options
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, color = "white"),  
    # Specify facetting options
    strip.background = element_rect(fill = "grey30", color = "grey10"),  
    strip.text.x = element_text(size = 12*0.8, color = "white"),  
    strip.text.y = element_text(size = 12*0.8, color = "white",angle = -90),  
    # Specify plot options
    plot.background = element_rect(color = "black", fill = "black"),  
    plot.title = element_text(size = 12*1.2, color = "white"),  
    plot.margin = unit(rep(1, 4), "lines")
  )


## How good or bad were the models reported?

mean.orig = mean(df.demo.filt.both$process.accuracy.model.relprob[df.demo.filt.both$type == 'Original'], na.rm = T)
se.orig = se(df.demo.filt.both$process.accuracy.model.relprob[df.demo.filt.both$type == 'Original'])
mean.obs = mean(df.demo.filt.both$process.accuracy.model.relprob[df.demo.filt.both$type == 'Observers'], na.rm = T)
se.obs = se(df.demo.filt.both$process.accuracy.model.relprob[df.demo.filt.both$type == 'Observers'])
ggplot(df.demo.filt.both, aes(x = process.accuracy.model.relprob, fill = type, group = type)) +
  geom_histogram(color = 'white', alpha = .7, bins = 25, position = position_identity()) +
  labs(x = "\nReported model BF",
       y = "# of subjects") +
  scale_y_continuous(breaks = NULL) +
  theme_black() +
  geom_vline(xintercept = mean.orig, linetype = 1, color = '#e3211c') +
  geom_vline(xintercept = mean.orig - se.orig, linetype = 'dashed', color = '#e3211c') +
  geom_vline(xintercept = mean.orig + se.orig, linetype = 'dashed', color = '#e3211c') +
  geom_vline(xintercept = mean.obs, linetype = 1, color = '#377eb8') +
  geom_vline(xintercept = mean.obs - se.obs, linetype = 'dashed', color = '#377eb8') +
  geom_vline(xintercept = mean.obs + se.obs, linetype = 'dashed', color = '#377eb8') +
  scale_fill_brewer(palette = 'Set1') +
  theme(legend.title = element_blank())
summary(lm(scale(process.accuracy.model.relprob) ~ type, df.demo.filt.both))

mean.orig = mean(df.demo.filt.both$process.accuracy.div[df.demo.filt.both$type == 'Original'], na.rm = T)
se.orig = se(df.demo.filt.both$process.accuracy.div[df.demo.filt.both$type == 'Original'])
mean.obs = mean(df.demo.filt.both$process.accuracy.div[df.demo.filt.both$type == 'Observers'], na.rm = T)
se.obs = se(df.demo.filt.both$process.accuracy.div[df.demo.filt.both$type == 'Observers'])
ggplot(df.demo.filt.both, aes(x = process.accuracy.div, fill = type, group = type)) +
  geom_histogram(color = 'white', alpha = .7, bins = 25, position = position_identity()) +
  labs(x = "\nHeuristic error",
       y = "# of subjects") +
  scale_y_continuous(breaks = NULL) +
  theme_black() +
  geom_vline(xintercept = mean.orig, linetype = 1, color = '#e3211c') +
  geom_vline(xintercept = mean.orig - se.orig, linetype = 'dashed', color = '#e3211c') +
  geom_vline(xintercept = mean.orig + se.orig, linetype = 'dashed', color = '#e3211c') +
  geom_vline(xintercept = mean.obs, linetype = 1, color = '#377eb8') +
  geom_vline(xintercept = mean.obs - se.obs, linetype = 'dashed', color = '#377eb8') +
  geom_vline(xintercept = mean.obs + se.obs, linetype = 'dashed', color = '#377eb8') +
  scale_fill_brewer(palette = 'Set1') +
  theme(legend.title = element_blank())
summary(lm(scale(process.accuracy.div) ~ type, df.demo.filt.both))
t.test(df.demo.filt$process.accuracy.div,
       df.demo.filt.obs$process.accuracy.div)

## by feature

# continuous version
ggplot(df.demo.filt.both, aes(x = actual.h1.prob, y = reported.h1, color = type, group = type)) +
  geom_point() +
  geom_smooth(method='lm') +
  theme_black() +
  labs(x = '\nProcess Question 1:\nModel-fitted probability of using\nmultiple atts (vs. a single att)',
       y = 'Reported probability of using\nmultiple atts (vs. a single att)') +
  scale_color_brewer(palette = 'Set1')
summary(lm(scale(reported.h1) ~ scale(actual.h1.prob) * type, df.demo.filt.both))
ggplot(df.demo.filt.both, aes(x = actual.h2.prob, y = reported.h2, color = type, group = type)) +
  geom_point() +
  geom_smooth(method='lm') +
  theme_black() +
  labs(x = '\nProcess Question 2:\nModel-fitted probability of using\ngraded weights (vs. binary weights)',
       y = 'Reported probability of using\ngraded weights (vs. binary weights)') +
  scale_color_brewer(palette = 'Set1')
summary(lm(scale(reported.h2) ~ scale(actual.h2.prob) * type, df.demo.filt.both))
ggplot(df.demo.filt.both, aes(x = actual.h3.prob, y = reported.h3, color = type, group = type)) +
  geom_point() +
  geom_smooth(method='lm') +
  theme_black() +
  labs(x = '\nProcess Question 3:\nModel-fitted probability of using\ngraded att values (vs. binary att values)',
       y = 'Reported probability of using\ngraded att values (vs. binary att values)') +
  scale_color_brewer(palette = 'Set1')
summary(lm(scale(reported.h3) ~ scale(actual.h3.prob) * type, df.demo.filt.both))

# dichotomized version
h1.grouped = df.demo.filt.both %>% group_by(type, actual.h1.prob.dichotomized) %>%
  summarize(Multiple = mean(reported.h1.dichotomized == 'Multiple'),
            Multiple.se = se.prop(reported.h1.dichotomized == 'Multiple'))
ggplot(h1.grouped, aes(x = actual.h1.prob.dichotomized, y = Multiple, group = type, fill = type)) +
  geom_col(position = dodge, color = 'white') +
  geom_errorbar(aes(ymin = Multiple - Multiple.se,
                    ymax = Multiple + Multiple.se),
                width = .2, position = dodge, color = 'white') +
  theme_black() +
  labs(x = '\nProcess Question 1:\nModel-fitting said multiple atts?',
       y = '% reporting using multiple atts',
       fill = '') +
  geom_hline(yintercept = .5, color = 'red', linetype = 'dashed') +
  scale_fill_brewer(palette = 'Set1') +
  scale_x_discrete(labels = c('False', 'True'))
summary(glm(reported.h1.dichotomized ~ actual.h1.prob.dichotomized * type, df.demo.filt.both, family = 'binomial'))

h2.grouped = df.demo.filt.both %>% group_by(type, actual.h2.prob.dichotomized) %>%
  summarize(Multiple = mean(reported.h2.dichotomized == 'Graded'),
            Multiple.se = se.prop(reported.h2.dichotomized == 'Graded'))
ggplot(h2.grouped, aes(x = actual.h2.prob.dichotomized, fill = type, group = type, y = Multiple)) +
  geom_col(position = dodge, color = 'white') +
  geom_errorbar(aes(ymin = Multiple - Multiple.se,
                    ymax = Multiple + Multiple.se),
                width = .2, position = dodge, color = 'white') +
  theme_black() +
  labs(x = '\nProcess Question 2:\nModel-fitting said graded weights?',
       y = '% reporting using graded weights',
       fill = '') +
  geom_hline(yintercept = .5, color = 'red', linetype = 'dashed')+
  scale_x_discrete(labels = c('False', 'True')) +
  scale_fill_brewer(palette = 'Set1')
summary(glm(reported.h2.dichotomized ~ actual.h2.prob.dichotomized * type, df.demo.filt.both, family = 'binomial'))

h3.grouped = df.demo.filt.both %>% group_by(type, actual.h3.prob.dichotomized) %>%
  summarize(Multiple = mean(reported.h3.dichotomized == 'Graded'),
            Multiple.se = se.prop(reported.h3.dichotomized == 'Graded'))
ggplot(h3.grouped, aes(x = actual.h3.prob.dichotomized, fill = type, group = type, y = Multiple)) +
  geom_col(color = 'white', position = dodge) +
  geom_errorbar(aes(ymin = Multiple - Multiple.se,
                    ymax = Multiple + Multiple.se),
                width = .2, color = 'white', position = dodge) +
  theme_black() +
  labs(x = '\nProcess Question 3:\nModel-fitting said graded att values?',
       y = '% reporting using graded att values',
       fill = '') +
  geom_hline(yintercept = .5, color = 'red', linetype = 'dashed')+
  scale_x_discrete(labels = c('False', 'True'))+
  scale_fill_brewer(palette = 'Set1')
summary(glm(reported.h3.dichotomized ~ actual.h3.prob.dichotomized * type, df.demo.filt.both, family = 'binomial'))

# combine all 3 together
df.questions1 = df.demo.filt.both %>%
  select(subject, type, reported.h1, reported.h2, reported.h3) %>%
  rename(q1 = reported.h1, q2 = reported.h2, q3 = reported.h3) %>%
  pivot_longer(!c(subject, type), names_to = 'question', values_to = 'reported')
df.questions2 = df.demo.filt.both %>%
  select(subject, type, actual.h1.prob, actual.h2.prob, actual.h3.prob) %>%
  pivot_longer(!c(subject, type), values_to = 'prob')
df.questions3 = df.demo.filt.both %>%
  select(subject, type, reported.h1.dichotomized, reported.h2.dichotomized, reported.h3.dichotomized) %>%
  pivot_longer(!c(subject, type), values_to = 'chosen')
df.questions4 = df.demo.filt.both %>%
  select(subject, type, actual.h1.prob.dichotomized, actual.h2.prob.dichotomized, actual.h3.prob.dichotomized) %>%
  pivot_longer(!c(subject, type), values_to = 'best.family')
df.questions = df.questions1
df.questions$prob = df.questions2$prob
df.questions$chosen = df.questions3$chosen
df.questions$best.family = df.questions4$best.family

ggplot(df.questions, aes(x = prob, y = reported, color = type, group = type)) +
  geom_point() +
  geom_smooth(method='lm') +
  theme_black() +
  labs(x = '\nProb. of using\nheuristic',
       y = 'Reported probability of\nusing heuristic') +
  scale_color_brewer(palette = 'Set1')
m1 = lmer(scale(reported) ~ scale(prob)*type + (scale(prob) | subject), df.questions)
summary(rePCA(m1))
m2 = lmer(scale(reported) ~ scale(prob)*type + (1 | subject), df.questions)
summary(m2)

q.grouped.rev = df.questions %>% group_by(type, best.family) %>%
  summarize(chosen.m = mean(chosen),
            chosen.se = se.prop(chosen))
ggplot(q.grouped.rev, aes(x = best.family, y = chosen.m, fill = type, group = type)) +
  geom_col(color = 'white', position = dodge) +
  geom_errorbar(aes(ymin = chosen.m - chosen.se,
                    ymax = chosen.m + chosen.se),
                width = .2, color = 'white', position = dodge) +
  theme_black() +
  labs(x = 'Used heuristic',
       y = '% reporting\nusing heuristic')  +
  scale_fill_brewer(palette = 'Set1')+
  theme(legend.title = element_blank())
m1 = glmer(chosen ~ best.family*type + (best.family | subject), df.questions, family = 'binomial')
summary(rePCA(m1))
m2 = glmer(chosen ~ best.family*type + (1 | subject), df.questions, family = 'binomial')
summary(m2)


# % chosen correctly
features.graph = df.demo.filt.both %>% select(subject,type,process.accuracy.h1.dichotomized, process.accuracy.h2.dichotomized, process.accuracy.h3.dichotomized, process.accuracy.model.dichotomized) %>%
  pivot_longer(!c(subject,type)) %>%
  group_by(type, name) %>%
  summarize(val = mean(value,na.rm=T),
            val.se = se(value)) %>%
  mutate(name = factor(name, c('process.accuracy.h1.dichotomized', 'process.accuracy.h2.dichotomized', 'process.accuracy.h3.dichotomized', 'process.accuracy.model.dichotomized')))
ggplot(features.graph %>% filter(name == 'process.accuracy.model.dichotomized'), aes(x = name, y = val, fill = type, group = type)) +
  geom_col(color='white', position = dodge) +
  geom_errorbar(aes(ymin = val - val.se, ymax = val + val.se), width = .2, color = 'white', position = dodge) +
  geom_segment(x = 0, y = .5, xend = 3.5, yend = 0.5, color = 'white', linetype = 'dashed') +
  geom_segment(x = 3.5, y = 1/6, xend = 4.5, yend = 1/6, color = 'white', linetype = 'dashed') +
  theme_black() +
  labs(x = '', y = '% reporting correct\nset of heuristics') +
  #scale_x_discrete(labels = c('Q1\n(one vs.\nmultiple attributes)', 'Q2\n(binary vs.\ngraded weights)', 'Q3\n(binary vs.\ngraded attributes)', 'Overall model')) +
  scale_x_discrete(labels = c('')) +
  scale_fill_brewer(palette = 'Set1')

summary(glm(process.accuracy.h1.dichotomized ~ type, df.demo.filt.both, family = 'binomial'))
summary(glm(process.accuracy.h2.dichotomized ~ type, df.demo.filt.both, family = 'binomial'))
summary(glm(process.accuracy.h3.dichotomized ~ type, df.demo.filt.both, family = 'binomial'))
summary(glm(process.accuracy.model.dichotomized ~ type, df.demo.filt.both, family = 'binomial'))

## weight accuracy

mean.orig = mean(df.demo.filt.both$weight.accuracy.averaged[df.demo.filt.both$type == 'Original'], na.rm = T)
se.orig = se(df.demo.filt.both$weight.accuracy.averaged[df.demo.filt.both$type == 'Original'])
mean.obs = mean(df.demo.filt.both$weight.accuracy.averaged[df.demo.filt.both$type == 'Observers'], na.rm = T)
se.obs = se(df.demo.filt.both$weight.accuracy.averaged[df.demo.filt.both$type == 'Observers'])
ggplot(df.demo.filt.both, aes(x = weight.accuracy.averaged, fill = type, group = type)) +
  geom_histogram(color = 'white', alpha = .7, bins = 25, position = position_identity()) +
  geom_vline(xintercept = mean.orig, linetype = 1, color = '#e3211c') +
  geom_vline(xintercept = mean.orig - se.orig, linetype = 'dashed', color = '#e3211c') +
  geom_vline(xintercept = mean.orig + se.orig, linetype = 'dashed', color = '#e3211c') +
  geom_vline(xintercept = mean.obs, linetype = 1, color = '#377eb8') +
  geom_vline(xintercept = mean.obs - se.obs, linetype = 'dashed', color = '#377eb8') +
  geom_vline(xintercept = mean.obs + se.obs, linetype = 'dashed', color = '#377eb8') +
  labs(x = 'Weight accuracy', y = '# of subjects') +
  #labs(x = '', y = '') +
  scale_y_continuous(breaks = NULL) +
  #scale_x_continuous(limits = c(-.2, 1.1), breaks = c(0, .5, 1)) +
  theme_black() +
  scale_fill_brewer(palette = 'Set1') +
  theme(legend.title = element_blank())

summary(lm(scale(weight.accuracy.averaged) ~ type, df.demo.filt.both))

mean.orig.multiatt = mean(df.demo.filt.both.multiatt$weight.accuracy.averaged[df.demo.filt.both$type == 'Original'], na.rm = T)
se.orig.multiatt = se(df.demo.filt.both.multiatt$weight.accuracy.averaged[df.demo.filt.both$type == 'Original'])
mean.obs.multiatt = mean(df.demo.filt.both.multiatt$weight.accuracy.averaged[df.demo.filt.both$type == 'Observers'], na.rm = T)
se.obs.multiatt = se(df.demo.filt.both.multiatt$weight.accuracy.averaged[df.demo.filt.both$type == 'Observers'])
ggplot(df.demo.filt.both.multiatt, aes(x = weight.accuracy.averaged, fill = type, group = type)) +
  geom_histogram(color = 'white', alpha = .7, bins = 25, position = position_identity()) +
  geom_vline(xintercept = mean.orig, linetype = 1, color = '#e3211c') +
  geom_vline(xintercept = mean.orig - se.orig, linetype = 'dashed', color = '#e3211c') +
  geom_vline(xintercept = mean.orig + se.orig, linetype = 'dashed', color = '#e3211c') +
  geom_vline(xintercept = mean.obs, linetype = 1, color = '#377eb8') +
  geom_vline(xintercept = mean.obs - se.obs, linetype = 'dashed', color = '#377eb8') +
  geom_vline(xintercept = mean.obs + se.obs, linetype = 'dashed', color = '#377eb8') +
  labs(x = 'Weight accuracy', y = '# of subjects') +
  #labs(x = '', y = '') +
  scale_y_continuous(breaks = NULL) +
  #scale_x_continuous(limits = c(-.2, 1.1), breaks = c(0, .5, 1)) +
  theme_black() +
  scale_fill_brewer(palette = 'Set1') +
  theme(legend.title = element_blank())

## moderators
# IQ
ggplot(df.demo.filt.both, aes(x = icar_num_correct, y = weight.accuracy.averaged, color = type)) +
  geom_point() +
  geom_smooth(method='lm', color = 'white') +
  theme_black() +
  labs(x = 'IQ score', y = 'Weight accuracy') +
  facet_wrap(~type,ncol = 1) +
  scale_color_brewer(palette = 'Set1')
ggplot(df.demo.filt.both, aes(x = icar_num_correct, y = process.accuracy.model.relprob, color = type)) +
  geom_point() +
  geom_smooth(method='lm', color = 'white') +
  theme_black() +
  labs(x = 'IQ score', y = 'Posterior prob. of reported model') +
  facet_wrap(~type,ncol = 1) +
  scale_fill_brewer(palette = 'Set1')
summary(lm(scale(weight.accuracy.averaged) ~ scale(icar_num_correct) + actual.model.fac + scale(actual.model.r2), df.demo.filt.obs))
summary(lm(scale(weight.accuracy.averaged) ~ scale(icar_num_correct) * type + actual.model.fac + scale(actual.model.r2), df.demo.filt.both))

summary(lm(scale(process.accuracy.qs) ~ scale(icar_num_correct) + actual.model.fac + scale(actual.model.r2), df.demo.filt.obs))
summary(lm(scale(process.accuracy.model.relprob) ~ scale(icar_num_correct) + actual.model.fac + scale(actual.model.r2), df.demo.filt.obs))
summary(lm(scale(process.accuracy.model.relprob) ~ scale(icar_num_correct)*type + actual.model.fac + scale(actual.model.r2), df.demo.filt.both))
summary(glm(process.accuracy.model.dichotomized ~ scale(icar_num_correct) + actual.model.fac + scale(actual.model.r2), df.demo.filt.obs, family = 'binomial'))
summary(glm(process.accuracy.model.dichotomized ~ scale(icar_num_correct)*type + actual.model.fac + scale(actual.model.r2), df.demo.filt.both, family = 'binomial'))

# meditation
meditation.graph3 = df.demo.filt.both %>%
  group_by(type, meditation.over5) %>%
  summarize(weight.accuracy.averaged.m = mean(weight.accuracy.averaged, na.rm = T),
            weight.accuracy.averaged.se = se(weight.accuracy.averaged),
            process.accuracy.model.prob.m = mean(process.accuracy.model.prob, na.rm = T),
            process.accuracy.model.prob.se = se(process.accuracy.model.prob),
            process.accuracy.model.dichotomized.m = mean(process.accuracy.model.dichotomized, na.rm = T),
            process.accuracy.model.dichotomized.se = se.prop(process.accuracy.model.dichotomized),
            process.accuracy.qs.m = mean(process.accuracy.qs),
            process.accuracy.qs.se = se(process.accuracy.qs))
ggplot(meditation.graph3, aes(x = meditation.over5, y = process.accuracy.model.prob.m, group = type, color = type)) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = process.accuracy.model.prob.m - process.accuracy.model.prob.se,
                    ymax = process.accuracy.model.prob.m + process.accuracy.model.prob.se,
                    width = .2))
ggplot(meditation.graph3, aes(x = meditation.over5, y = process.accuracy.model.dichotomized.m, group = type, color = type)) +
  geom_point(size = 5, color = 'white') +
  geom_errorbar(aes(ymin = process.accuracy.model.dichotomized.m - process.accuracy.model.dichotomized.se,
                    ymax = process.accuracy.model.dichotomized.m + process.accuracy.model.dichotomized.se),
                width = .2, color = 'white') +
  #scale_y_continuous(limits = c(0.1, 0.4), breaks = c(0.1, 0.2, 0.3, 0.4), labels = c('10', '20', '30', '40')) +
  labs(x = '\n>5 yrs meditation experience', y = '% reporting correct model') +
  theme_black()
ggplot(meditation.graph3, aes(x = meditation.over5, y = weight.accuracy.averaged.m, group = type, color = type)) +
  geom_point(size = 5, color = 'white') +
  geom_errorbar(aes(ymin = weight.accuracy.averaged.m - weight.accuracy.averaged.se,
                    ymax = weight.accuracy.averaged.m + weight.accuracy.averaged.se),
                width = .2, color = 'white') +
  labs(x = '\n>5 yrs meditation experience', y = 'Weight accuracy') +
  theme_black()

summary(lm(scale(weight.accuracy.averaged) ~ meditation.over5 +
             actual.model.fac + actual.model.r2, data = df.demo.filt))
summary(glm(process.accuracy.model.dichotomized ~ meditation.over5 +
              actual.model.fac + actual.model.r2, data = df.demo.filt, family = 'binomial'))

# all
m.mods = lm(weight.accuracy.averaged ~ decisionstyle + mindfulness + sris.tendency + sris.insight + maia + acs.shifting + acs.focusing +
              icar_num_correct + gender + edu.num + meditation.over5 +
              alerting + orienting + exec +
              choice_domain + confidence + consistency1 + consistency2 + choice_exp_num +
              actual.model.fac + actual.model.r2 + actual.model.prob,
            data = df.demo.filt.obs)
summary(m.mods)

m.mods2 = lm(process.accuracy.model.prob ~ decisionstyle + mindfulness + sris.tendency + sris.insight + maia + acs.shifting + acs.focusing +
               icar_num_correct + gender + edu.num + meditation_exp1 +
               alerting + orienting + exec +
               choice_domain + confidence + consistency1 + consistency2 + choice_exp_num +
               actual.model.fac + actual.model.r2 + actual.model.prob,
             data = df.demo.filt)
summary(m.mods2)

m.mods3 = glm(process.accuracy.model.dichotomized ~ decisionstyle + mindfulness + sris.tendency + sris.insight + maia + acs.shifting + acs.focusing +
                icar_num_correct + gender + edu.num + meditation_exp1 +
                alerting + orienting + exec +
                choice_domain + confidence + consistency1 + consistency2 + choice_exp_num +
                actual.model.fac + actual.model.r2 + actual.model.prob,
              data = df.demo.filt, family = 'binomial')
summary(m.mods3)

# do matching shit
test = df.demo.filt.both %>% filter(!is.na(target_id_num)) %>% count(target_id_num)

df.demo.filt$h1.correct.obs = NA
for (i in 1:nrow(df.demo.filt)) {
  subj = df.demo.filt$subject[i]
  df.obs.cur = df.demo.filt.obs %>% filter(target_id == subj)
  if (nrow(df.obs.cur) > 0) {
    df.demo.filt$h1.correct.obs[i] = mean(df.obs.cur$h1.correct)
  }
}

mean(df.demo.filt$h1.correct)
mean(df.demo.filt$h1.correct.obs, na.rm = T)
wilcox.test(as.numeric(df.demo.filt$h1.correct), df.demo.filt$h1.correct.obs, paired = T)

# Do split-half analyses --------------------------------------------------

## half1 vs half2, or even vs. odd
analysis_types = c('split-half', 'even-odd')
analysis_type = analysis_types[2]

if (analysis_type == 'split-half') {
  load(paste0(filepath, '/modeling-output/half1/dfs.rdata'))
  load(paste0(filepath, '/modeling-output/half2/dfs.rdata'))
} else {
  load(paste0(filepath, '/modeling-output/even/dfs.rdata'))
  load(paste0(filepath, '/modeling-output/odd/dfs.rdata'))
  
  df.demo.filt.half1 = df.demo.filt.even
  df.demo.filt.half2 = df.demo.filt.odd
  df.attributes.filt.half1 = df.attributes.filt.even
  df.attributes.filt.half2 = df.attributes.filt.odd
}

df.demo.filt.half1$type = 'Half1'
df.demo.filt.half2$type = 'Half2'
df.attributes.filt.half1$type = 'Half1'
df.attributes.filt.half2$type = 'Half2'

df.demo.filt.splithalf = full_join(df.demo.filt.half1, df.demo.filt.half2)
df.demo.filt.splithalf = df.demo.filt.splithalf %>% mutate(type = factor(type, c('Half1', 'Half2')))
df.attributes.filt.splithalf = full_join(df.attributes.filt.half1, df.attributes.filt.half2)
df.attributes.filt.splithalf = df.attributes.filt.splithalf %>% mutate(type = factor(type, c('Half1', 'Half2')))

df.demo.filt.splithalf.grp = df.demo.filt.splithalf %>%
  group_by(type) %>%
  summarize(actual.h1.prob.m = mean(actual.h1.prob), actual.h1.prob.se = se(actual.h1.prob),
            actual.h2.prob.m = mean(actual.h2.prob), actual.h2.prob.se = se(actual.h2.prob),
            actual.h3.prob.m = mean(actual.h3.prob), actual.h3.prob.se = se(actual.h3.prob),
            actual.qs.prob.m = mean(c(actual.h1.prob.m, actual.h2.prob.m, actual.h3.prob.m)),
            process.accuracy.div.m = mean(process.accuracy.div), process.accuracy.div.se = se(process.accuracy.div),
            weight.accuracy.averaged.m = mean(weight.accuracy.averaged), weight.accuracy.averaged.se = se(weight.accuracy.averaged),
  )

# heuristic questions
splithalf.actual.h1.prob = df.demo.filt.splithalf %>%
  select(subject, type, actual.h1.prob) %>%
  pivot_wider(names_from = type,
            values_from = actual.h1.prob)
splithalf.actual.h2.prob = df.demo.filt.splithalf %>%
  select(subject, type, actual.h2.prob) %>%
  pivot_wider(names_from = type,
              values_from = actual.h2.prob)
splithalf.actual.h3.prob = df.demo.filt.splithalf %>%
  select(subject, type, actual.h3.prob) %>%
  pivot_wider(names_from = type,
              values_from = actual.h3.prob)
splithalf.actual.qs.prob = full_join(full_join(splithalf.actual.h1.prob,
                                     splithalf.actual.h2.prob),
                                     splithalf.actual.h3.prob)

ggplot(splithalf.actual.qs.prob, aes(x = Half1, y = Half2)) +
  geom_point(color = 'white') +
  geom_smooth(method = 'lm', color = 'white') +
  labs(x = '\nProbability of using heuristic\non even trials',
       y = 'Probability of using\nheuristic on odd trials') +
  theme_black()
cor.test(splithalf.actual.qs.prob$Half1, splithalf.actual.qs.prob$Half2)

# overall model
splithalf.actual.model = df.demo.filt.splithalf %>%
  select(subject, type, actual.model.fac) %>%
  pivot_wider(names_from = type,
              values_from = actual.model.fac)
splithalf.actual.model.heat.normed = splithalf.actual.model %>%
  group_by(Half1, Half2) %>%
  summarize(num.subj = n()) %>%
  group_by(Half1) %>%
  mutate(num.subj.normed = num.subj / sum(num.subj))
ggplot(splithalf.actual.model.heat.normed, aes(x = Half1, y = Half2,
                                    fill = num.subj.normed)) +
  geom_tile() +
  geom_text(aes(label = round(num.subj.normed, 2))) +
  labs(y = '\nBest-fitting model in half 2', x = 'Best-fitting model in half 1') +
  #scale_fill_brewer(palette = 'YlOrRd') +
  scale_fill_continuous(limits = c(0,1), low = 'black', high = 'white') +
  guides(fill = guide_colorbar(title = '% of subjects')) +
  theme_black()

# weights
splithalf.fitted.weight.Full = df.attributes.filt.splithalf %>%
  select(subject, attribute, type, fitted.weight.Full) %>%
  pivot_wider(names_from = type,
              values_from = fitted.weight.Full)
ggplot(splithalf.fitted.weight.Full, aes(x = Half1, y = Half2)) +
  geom_point(color = 'white') +
  geom_smooth(method = 'lm', color = 'white') +
  labs(x = '\nWeight on even trials',
       y = 'Weight on odd trials') +
  theme_black()
cor.test(splithalf.fitted.weight.Full$Half1, splithalf.fitted.weight.Full$Half2)

# weight accuracy
splithalf.weight.accuracy.averaged = df.demo.filt.splithalf %>%
  select(subject, type, weight.accuracy.averaged) %>%
  pivot_wider(names_from = type,
              values_from = weight.accuracy.averaged)
ggplot(splithalf.weight.accuracy.averaged, aes(x = Half1, y = Half2)) +
  geom_point(color = 'white') +
  geom_smooth(method = 'lm', color = 'white') +
  labs(x = '\nWeight accuracy in\nfirst half',
       y = 'Weight accuracy in\nsecond half') +
  theme_black()
cor.test(splithalf.weight.accuracy.averaged$Half1, splithalf.weight.accuracy.averaged$Half2)

# heuristic accuracy
splithalf.process.accuracy.div = df.demo.filt.splithalf %>%
  select(subject, type, process.accuracy.div) %>%
  pivot_wider(names_from = type,
              values_from = process.accuracy.div)
ggplot(splithalf.process.accuracy.div, aes(x = Half1, y = Half2)) +
  geom_point(color = 'white') +
  geom_smooth(method = 'lm', color = 'white') +
  labs(x = '\nHeuristic accuracy\nin first half',
       y = 'Heuristic accuracy\nin second half') +
  theme_black()
cor.test(splithalf.process.accuracy.div$Half1, splithalf.process.accuracy.div$Half2)

# Save data ---------------------------------------------------------------

save.image(paste0(filepath, 'analysis_output.rdata'))

# Test-retest -------------------------------------------------------------

# should include df.demo.movies and df.demo.rand
# df.demo.rand is df.demo.filt from rand, but without filtering people who did movies first
# df.demo.movies is df.demo.filt from movies (i.e., this script), but without filtering people who did rand first
load('testretest.rdata')

subj.both = intersect(df.demo.movies$subject, df.demo.rand$subject)
df.demo.movies.trt = df.demo.movies %>% filter(subject %in% subj.both) %>%
  dplyr::select(subject, weight.accuracy.averaged,
         process.accuracy.model.dichotomized, process.accuracy.model.relprob,
         process.accuracy.h1, process.accuracy.h2, process.accuracy.h3, process.accuracy.qs,
         which_version_first) %>%
  mutate(which_version = 'movies')
df.demo.rand.trt = df.demo.rand %>% filter(subject %in% subj.both) %>%
  dplyr::select(subject, weight.accuracy.averaged,
                process.accuracy.model.dichotomized, process.accuracy.model.relprob,
                process.accuracy.h1, process.accuracy.h2, process.accuracy.h3, process.accuracy.qs,
                which_version_first) %>%
  mutate(which_version = 'rand')

df.demo.trt = bind_rows(df.demo.movies.trt, df.demo.rand.trt)
df.demo.trt.wide = df.demo.trt %>% pivot_wider(names_from = which_version, values_from = c(weight.accuracy.averaged,
                                                                                           process.accuracy.model.dichotomized, process.accuracy.model.relprob,
                                                                                           process.accuracy.h1, process.accuracy.h2, process.accuracy.h3, process.accuracy.qs,))

cor.test(df.demo.movies.trt$weight.accuracy.averaged, df.demo.rand.trt$weight.accuracy.averaged)
ggplot(df.demo.trt.wide, aes(x = weight.accuracy.averaged_movies, y = weight.accuracy.averaged_rand)) +
  geom_point() +
  geom_smooth(method='lm')

summary(lm(weight.accuracy.averaged_rand ~ weight.accuracy.averaged_movies, data = df.demo.trt.wide))

with(df.demo.trt.wide, table(process.accuracy.model.dichotomized_movies, process.accuracy.model.dichotomized_rand))
ggplot(df.demo.trt.wide, aes(x = process.accuracy.model.relprob_movies, y = process.accuracy.model.relprob_rand)) +
  geom_point() +
  geom_smooth(method='lm')

ggplot(df.demo.trt.wide, aes(x = process.accuracy.qs_movies, y = process.accuracy.qs_rand)) +
  geom_point() +
  geom_smooth(method='lm')

# compare first to second
for (i in 1:nrow(df.demo.trt.wide)) {
  if (df.demo.trt.wide$which_version_first[i] == 'movies') {
    df.demo.trt.wide$weight.accuracy.averaged.first[i] = df.demo.trt.wide$weight.accuracy.averaged_movies[i]
    df.demo.trt.wide$weight.accuracy.averaged.second[i] = df.demo.trt.wide$weight.accuracy.averaged_rand[i]
    df.demo.trt.wide$process.accuracy.model.dichotomized.first[i] = df.demo.trt.wide$process.accuracy.model.dichotomized_movies[i]
    df.demo.trt.wide$process.accuracy.model.dichotomized.second[i] = df.demo.trt.wide$process.accuracy.model.dichotomized_rand[i]
    df.demo.trt.wide$process.accuracy.model.relprob.first[i] = df.demo.trt.wide$process.accuracy.model.relprob_movies[i]
    df.demo.trt.wide$process.accuracy.model.relprob.second[i] = df.demo.trt.wide$process.accuracy.model.relprob_rand[i]
  } else {
    df.demo.trt.wide$weight.accuracy.averaged.first[i] = df.demo.trt.wide$weight.accuracy.averaged_rand[i]
    df.demo.trt.wide$weight.accuracy.averaged.second[i] = df.demo.trt.wide$weight.accuracy.averaged_movies[i]
    df.demo.trt.wide$process.accuracy.model.dichotomized.first[i] = df.demo.trt.wide$process.accuracy.model.dichotomized_rand[i]
    df.demo.trt.wide$process.accuracy.model.dichotomized.second[i] = df.demo.trt.wide$process.accuracy.model.dichotomized_movies[i]
    df.demo.trt.wide$process.accuracy.model.relprob.first[i] = df.demo.trt.wide$process.accuracy.model.relprob_rand[i]
    df.demo.trt.wide$process.accuracy.model.relprob.second[i] = df.demo.trt.wide$process.accuracy.model.relprob_movies[i]
  }
}

df.demo.trt = df.demo.trt %>%
  mutate(which_order = factor(which_version_first == which_version, c(T,F), c('First', 'Second')))

df.demo.trt.graph = df.demo.trt %>%
  group_by(which_order, which_version_first) %>%
  summarize(weight.accuracy.averaged.m = mean(weight.accuracy.averaged),
            weight.accuracy.averaged.se = se(weight.accuracy.averaged),
            process.accuracy.model.dichotomized.m = mean(process.accuracy.model.dichotomized),
            process.accuracy.model.dichotomized.se = se.prop(process.accuracy.model.dichotomized),
            process.accuracy.model.relprob.m = mean(process.accuracy.model.relprob),
            process.accuracy.model.relprob.se = se(process.accuracy.model.relprob))

ggplot(df.demo.trt.graph, aes(x = which_order, y = weight.accuracy.averaged.m)) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = weight.accuracy.averaged.m - weight.accuracy.averaged.se,
                    ymax = weight.accuracy.averaged.m + weight.accuracy.averaged.se),
                width = .2) +
  labs(x = 'Which test', y = 'Weight accuracy') +
  facet_wrap(~which_version_first)

ggplot(df.demo.trt.graph, aes(x = which_order, y = process.accuracy.model.relprob.m)) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = process.accuracy.model.relprob.m - process.accuracy.model.relprob.se,
                    ymax = process.accuracy.model.relprob.m + process.accuracy.model.relprob.se),
                width = .2) +
  labs(x = 'Which test', y = 'Process accuracy') +
  facet_wrap(~which_version_first)