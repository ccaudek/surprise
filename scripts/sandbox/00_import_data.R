#' Surprise Experiment
#' Reads all files in the data directory.
#'
#' File: 00_import_data.R
#'
#' This script was last modified on "Mon May 21 10:34:52 2018"




library("tidyverse")


# Data ingestion ----------------------------------------------------------


# Change exp_flag to analyze each of the two experiments!

exp_flag <- "SURPRISE"

# Get the list of files
if (exp_flag == "SURPRISE") {
  filenames <- Sys.glob("../data/surprise_videos/*.dat") # both Lucrezia's and Elena's data
  n_files <- length(filenames)
} else if (exp_flag == "CONTROL") { # control experiment
  filenames <- Sys.glob("../data/no_surprise_videos/*.dat")
  n_files <- length(filenames)
} else {
  message("\nWhich experiment do you want to analyze?\n")
}




mylist <- list()
for (i in 1:n_files) {
  mylist[[i]] <- read.table(filenames[i], header = TRUE, row.names=NULL)
}


# remove row names
mylist2 <- list()
for (j in 1:length(mylist)) {
  foo <- mylist[[j]]
  if (length(names(foo)) == 20) {
    mylist2[[j]] <- foo[-1]
    print(j, length(names(mylist[[j]])))
  } else {
    mylist2[[j]] <- foo
  }
}

# remove "block" and "isSurprise"
mylist3 <- list()
for (j in 1:length(mylist2)) {
  foo <- mylist2[[j]]
  if (length(names(foo)) == 19) {
    mylist3[[j]] <- foo[-c(3, 19)]
    print(j, length(names(mylist[[j]])))
  } else {
    mylist3[[j]] <- foo
  }
}

# check it out
for(i in 1:length(mylist3)) {
  print(length(names(mylist3[[i]])))
}


# create data.frame with all the individual data.
df <- do.call(what = rbind, args = mylist3)

message("\nNumber of subjects:\n")
# number of subjects
length(unique(df$subject_name))




# Exit message ------------------------------------------------------------


n_subj <- length(unique(df$subject_name))
cat("\nDone! The data.frame df has been created.\nNumber of subjects: ", n_subj)



