
# Only works in RStudio -- otherwise you have to set the path manually
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
 
require(dplyr)
require(lubridate)

df = read.csv('test.csv', col.names = c('movies','movies.compl','rand','rand.compl')) %>%
  mutate(movies.compl = mdy_hm(movies.compl), rand.compl = mdy_hm(rand.compl))

x = intersect(df$movies, df$rand)

movies.compl = numeric(length(x))
rand.compl = numeric(length(x))
for (i in 1:length(x)) {
  movies.compl[i] = df$movies.compl[df$movies == x[i]]
  rand.compl[i] = df$rand.compl[df$rand == x[i]]
}

movies.first = x[movies.compl < rand.compl]
rand.first = x[rand.compl < movies.compl]



demo_movies = read.csv('demo.csv') %>% mutate(stamp = ymd_hms(stamp))
demo_rand = read.csv('demo_rand.csv') %>% mutate(stamp = ymd_hms(stamp))
demo_movies_followup = read.csv('demo_followup.csv') %>% mutate(stamp = ymd_hms(stamp))
demo_rand_followup = read.csv('demo_followup_rand.csv') %>% mutate(stamp = ymd_hms(stamp))

demo_movies_subj = unique(demo_movies$subject)
demo_rand_subj = unique(demo_rand$subject)
subj.both = intersect(demo_movies_subj, demo_rand_subj)

demo_movies$did_followup = NA
demo_movies$did_followup_rand = NA
demo_movies$in_rand = NA
demo_movies$which_version_first = NA
for (i in 1:nrow(demo_movies)) {
  subj = demo_movies$subject[i]
  if (subj %in% demo_movies_followup$subject) {
    demo_movies$did_followup[i] = T
  } else {
    demo_movies$did_followup[i] = F
  }
  
  
  if (subj %in% demo_rand_followup$subject) {
    demo_movies$did_followup_rand[i] = T
  } else {
    demo_movies$did_followup_rand[i] = F
  }
  
  if (subj %in% demo_rand$subject) {
    demo_movies$in_rand[i] = T
    demo_stamp = demo_movies$stamp[i]
    rand_stamp = demo_rand$stamp[demo_rand$subject == subj]
    if (demo_stamp < rand_stamp) {
      demo_movies$which_version_first[i] = 'movies'
    } else {
      demo_movies$which_version_first[i] = 'rand'
    }
  } else {
    demo_movies$in_rand[i] = F
  }
}

# we want: ppl who only did movies, or who did movies before they did rand, but didn't do followup
followup_list2 = demo_movies$assignmentId[
  (!demo_movies$in_rand | demo_movies$which_version_first == 'movies') &
    !demo_movies$did_followup]
paste(followup_list2, collapse = ",")

sum(!demo_movies$in_rand | demo_movies$which_version_first == 'movies')

demo_rand$did_followup = NA
demo_rand$did_followup_movies = NA
demo_rand$in_movies = NA
demo_rand$which_version_first = NA
for (i in 1:nrow(demo_rand)) {
  subj = demo_rand$subject[i]
  if (subj %in% demo_rand_followup$subject) {
    demo_rand$did_followup[i] = T
  } else {
    demo_rand$did_followup[i] = F
  }
  
  
  if (subj %in% demo_movies_followup$subject) {
    demo_rand$did_followup_movies[i] = T
  } else {
    demo_rand$did_followup_movies[i] = F
  }
  
  if (subj %in% demo_movies$subject) {
    demo_rand$in_movies[i] = T
    rand_stamp = demo_rand$stamp[i]
    movies_stamp = demo_movies$stamp[demo_movies$subject == subj]
    if (rand_stamp < movies_stamp) {
      demo_rand$which_version_first[i] = 'rand'
    } else {
      demo_rand$which_version_first[i] = 'movies'
    }
  } else {
    demo_rand$in_movies[i] = F
  }
}

# we want: ppl who only did movies, or who did movies before they did rand, but didn't do either followup
followup_list2 = demo_rand$assignmentId[
  (!demo_rand$in_movies | demo_rand$which_version_first == 'rand') &
    !demo_rand$did_followup & !demo_rand$did_followup_movies]
paste(followup_list2, collapse = ",")

sum((!demo_rand$in_movies | demo_rand$which_version_first == 'rand'))



# we want: people who only did movies
retest_list_movies = demo_movies$assignmentId[!demo_movies$in_rand]
paste(retest_list_movies, collapse = ",")

# we want: people who only did rand
retest_list_rand = demo_rand$assignmentId[!demo_rand$in_movies]
paste(retest_list_rand, collapse = ",")
