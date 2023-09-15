# Script name: 01_import_data.R
# Project: Surprise
# Script purpose: import data
# @author: Corrado Caudek <corrado.caudek@unifi.it>
# Date Created: Tue Feb  9 13:03:02 2021
# Last Modified Date: Tue Sep 20 06:43:04 2022
# 
# Notes: This script must be run first. Then the resulting df dataframe must be 
#        processed by 10_tidying_data.R. The process must be repeated twice: 
#        once for the control group, once for the surprise group.


library("tidyverse")
library("here")

source(here("libraries", "fnct_surprise.R"))


# Data ingestion ----------------------------------------------------------

# Change exp_flag to analyze each of the two experiments!

exp_flag <- "SURPRISE" #"SURPRISE" # CONTROL

# Get file names and number of files.
if (exp_flag == "SURPRISE") {
  filenames <- Sys.glob(here("data", "surprise_videos", "*.dat"))
  n_files <- length(filenames)
} else if (exp_flag == "CONTROL") { # control experiment
  filenames <- Sys.glob(here("data", "no_surprise_videos", "*.dat"))
  # filenames <- Sys.glob(here("data", "giulia", "*.dat"))
  n_files <- length(filenames)
} else {
  message("\nWhich experiment do you want to analyze?\n")
}

# Get list of data.frames.
mylist <- list()
for (i in 1:n_files) {
  mylist[[i]] <- read.table(filenames[i], header = TRUE, row.names=NULL)
}

# Get info about number of columns for each file.
id <- rep(NA, n_files)
n_col <- rep(NA, n_files)
for (i in 1:length(mylist)) {
  id[i] <- as.character(unique(mylist[[i]]$subject_name))
  n_col[i] <- length(names(mylist[[i]]))
}

# Use files with 17 columns. Use only the files with 17 columns.  
# For the control group there are no problems.
# Ci sono file con 17, 19 e 20 colonne.  In alcuni file i nomi delle colonne
# sono sbagliati!
list_17 <- list()
index <- 1
for (i in 1:length(mylist)) {
  foo <- mylist[[i]]
  if (length(names(foo)) == 17) {
    list_17[[index]] <- foo
    index <- index + 1
  }
}
length(list_17)

new_list_17 <- list()
for (i in 1:length(list_17)) {
  
  foo <- list_17[[i]]
  n <- length(foo$subject_name)
  
  df <- data.frame(
    subj_name = rep(as.character(unique(foo$subject_name)), n))
  df$subject_number <- rep(unique(foo$subjectNumber), n)
  df$trial <- foo$iTrial
  df$target_or <- foo$targetOrientation
  df$flanker_or <- foo$flankersOrientation
  df$resp <- foo$response
  df$correct <- foo$correct
  df$movie_id <- foo$magicShows
  df$rt <- foo$rt
  df$date <- foo$date
  df$time_of_day <- foo$time_of_day
  delta_x <- foo$deltaX
  delta_y <- foo$deltaY
  df$is_surprise_clip <- foo$surpriseClip
  df$is_clip_trial <- foo$isClipTrial
  
  new_list_17[[i]] <- df
  rm(df)
}

tot_df <- do.call(what = rbind, args = new_list_17)

# Recover clip type, surprise / no surprise
baa <- as.character(tot_df$movie_id)
tot_df$movie_id2 <- sapply(strsplit(baa, "/stim/"), "[", 2)
summary(factor(tot_df$movie_id2))

# create numeric id index for each subject
tot_df$subj_id <- as.numeric(as.factor(tot_df$subj_name))
n_sub <- length(unique(tot_df$subj_id))

tot_df %>% 
  group_by(subj_id) %>% 
  summarise(
    n = n()
  ) %>% 
  as.data.frame()

# select the first 320 rows of each subject
foo <- tot_df %>% 
  arrange(desc(subj_id)) %>% 
  group_by(subj_id) %>% 
  slice(1:320)

if (exp_flag == "CONTROL") {
  # add remaining subjects
  additional_subjects <- add_control_subjects(tot_df)
  last <- bind_rows(list(foo, additional_subjects))
  # last <- temp[complete.cases(temp), ]
  length(unique(last$subj_name))
  
  # create new numeric id index for each subject
  last$subj_id <- as.numeric(as.factor(last$subj_name))
  n_sub <- length(unique(last$subj_id))
  n_sub
  
  last %>% 
    group_by(subj_id) %>% 
    summarise(
      n = n(),
      nn = length(unique(subject_number))
    ) %>% 
    as.data.frame()
  # some subjects have less than 320 rows!
  
  unique(last$subj_id)
} else {
  last <- foo
}

# add block
blk <- rep(1:4, each = 80)
ss <- list()
for (i in 1:length(unique(last$subj_id))) {
  one_subj <- last %>% 
    dplyr::filter(subj_id == i)
  one_subj$block <- blk[1:nrow(one_subj)]
  ss[[i]] <- one_subj
}

# create data.frame with all the individual data
df <- do.call(what = rbind, args = ss)

if (exp_flag == "CONTROL" ) {
  df_controls_giulia <- get_data_giulia()
  df <- bind_rows(df, df_controls_giulia)
}


# Exit message
n_subj <- length(unique(df$subj_name))
cat("\nDone! The data.frame df has been created.\nNumber of subjects: ", n_subj)



