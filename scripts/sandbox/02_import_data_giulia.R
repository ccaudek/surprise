# Script name: 02_import_data_giulia.R
# Project: Surprise
# Script purpose: import control data run by Giulia
# @author: Corrado Caudek <corrado.caudek@unifi.it>
# Date Created: Tue Feb  9 13:03:02 2021
# Last Modified Date: Tue Sep 20 06:43:04 2022
# 
# Notes: T


library("tidyverse")
library("here")

source(here("libraries", "fnct_surprise.R"))


get_data_giulia <- function() {
  
  # Get file names and number of files of the control experiment (Giulia)
  filenames <- Sys.glob(here("data", "no_surprise_video_giulia", "*.dat"))
  n_files <- length(filenames)
  
  # Get list of data.frames.
  mylist <- list()
  for (i in 1:n_files) {
    mylist[[i]] <- read.table(filenames[i], header = TRUE, row.names=NULL)
  }
  
  new_list <- list()
  for (i in 1:length(mylist)) {
    
    foo <- mylist[[i]]
    
    subj_name_part_1 <- unique(foo$row.names)
    subj_name_part_2 <- unique(foo$subject_name)
    subj_name_string <- paste0(subj_name_part_1, subj_name_part_2)
    
    n <- length(foo$subject_name)
    
    df <- data.frame(
      subj_name = rep(subj_name_string, n),
      subject_number = 1000, # TODO
      trial = foo$iTrial,
      target_or = foo$flankersOrientation,
      flanker_or = foo$response,
      resp = foo$correct,
      correct = foo$onscreen_stim_duration,
      movie_id = foo$rt,
      rt = foo$date,
      date = foo$time_of_day,
      time_of_day = foo$deltaX,
      is_surprise_clip = foo$isClipTrial,
      is_clip_trial = foo$isSurprise,
      movie_id2 = substring(foo$rt, 57),
      subj_id = 1000, # TODO
      block = foo$block
    )
    
    new_list[[i]] <- df
    rm(df)
  }
  
  tot_df <- do.call(what = rbind, args = new_list)
  
  tot_df$subject_number <- 5000 + as.numeric(as.factor(as.character(tot_df$subj_name)))
  tot_df$subj_id <- 5000 + as.numeric(as.factor(as.character(tot_df$subj_name)))
  
  tot_df
}




# Data ingestion ----------------------------------------------------------

# Get file names and number of files of the control experiment (Giulia)
filenames <- Sys.glob(here("data", "no_surprise_video_giulia", "*.dat"))
n_files <- length(filenames)


# Get list of data.frames.
mylist <- list()
for (i in 1:n_files) {
  mylist[[i]] <- read.table(filenames[i], header = TRUE, row.names=NULL)
}

# # Get info about number of columns for each file.
# id <- rep(NA, n_files)
# n_col <- rep(NA, n_files)
# for (i in 1:length(mylist)) {
#   id[i] <- as.character(unique(mylist[[i]]$subject_name))
#   n_col[i] <- length(names(mylist[[i]]))
# }
# 
# list_20 <- list()
# index <- 1
# for (i in 1:length(mylist)) {
#   foo <- mylist[[i]]
#   if (length(names(foo)) == 20) {
#     list_20[[index]] <- foo
#     index <- index + 1
#   }
# }
# length(list_20)

new_list <- list()
for (i in 1:length(mylist)) {
  
  foo <- mylist[[i]]
  
  subj_name_part_1 <- unique(foo$row.names)
  subj_name_part_2 <- unique(foo$subject_name)
  subj_name_string <- paste0(subj_name_part_1, subj_name_part_2)

  n <- length(foo$subject_name)
  
  df <- data.frame(
    subj_name = rep(subj_name_string, n),
    subject_number = 1000, # TODO
    trial = foo$iTrial,
    target_or = foo$flankersOrientation,
    flanker_or = foo$response,
    resp = foo$correct,
    correct = foo$onscreen_stim_duration,
    movie_id = foo$rt,
    rt = foo$date,
    date = foo$time_of_day,
    time_of_day = foo$deltaX,
    is_surprise_clip = foo$isClipTrial,
    is_clip_trial = foo$isSurprise,
    movie_id2 = substring(foo$rt, 57),
    subj_id = 1000, # TODO
    block = foo$block
  )
  
  new_list[[i]] <- df
  rm(df)
}

tot_df <- do.call(what = rbind, args = new_list)

tot_df$subject_number <- 5000 + as.numeric(as.factor(as.character(tot_df$subj_name)))
tot_df$subj_id <- 5000 + as.numeric(as.factor(as.character(tot_df$subj_name)))

# Exit message
n_subj <- length(unique(df$subj_name))
cat("\nDone! The data.frame df has been created.\nNumber of subjects: ", n_subj)


