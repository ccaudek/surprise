#' ---
#' title: "Models' fit: control condition"
#' output:
#'   pdf_document:
#'     keep_tex: true
#' ---

# Script name: 20_variational_inference.R
# Project: surprise with flanker task
# Script purpose: brms analysis of the congruenty effect
# @author: Corrado Caudek <corrado.caudek@unifi.it>
# Date Created: Tue Sep 12 08:39:19 2023
# Last Modified Date: Fri Sep 15 08:10:18 2023


library("here")

suppressPackageStartupMessages({
    library("tidyverse")
    library("brms")
    library("cmdstanr")
    library("reshape")
    library("devtools")
    library("mice")
    library("tidybayes")
    library("emmeans")
    library("broom.mixed")
    library("rstanarm")
    library("patchwork")
    library("flankr")
})

theme_set(bayesplot::theme_default(base_family = "sans", base_size = 14))
set.seed(123)

source(here("libraries", "functions.R"))


#' Import complete data set

# Import the data that have been created by the previous scripts.
# The data have been created with the 01_, 10_, 11_ scripts in
# the present directory.
data <- get_data()

#' Tidy data
data_tidy <- tidy_flanker(data)

#' Perform some participants' flanker checks.
flanker_accuracy_overall <- get_flanker_accuracy(data_tidy, overall = TRUE)

# Get a list of participants who scored below 80% accuracy.
accuracy_removal <- flanker_accuracy_overall |>
    filter(accuracy < 0.80) |>
    pull(subj_id)

length(accuracy_removal)
# 1

# Remove the <80% accuracy participants from the flanker data.
flanker_data <- data_tidy |>
    filter(!subj_id %in% accuracy_removal)

flanker_data |>
    group_by(experiment, is_surprise_clip) |>
    summarize(
        n = n_distinct(subj_id)
    )

sort(unique(flanker_data$bysub_n))

flanker_data_clean <- flanker_data |>
    dplyr::filter(bysub_n > 240 & bysub_n < 321)


# Check number of subjects by condition.
flanker_data_clean |>
    group_by(experiment, is_surprise_clip) |>
    summarize(
        n = n_distinct(subj_id)
    )

#' Data wrangling for models' fit

hist(flanker_data_clean$rt)

#' select experiment

modeling_data_cntl <- flanker_data_clean |>
    dplyr::filter(experiment == "control")

modeling_data_surprise <- flanker_data_clean |>
    dplyr::filter(experiment == "surprise")

# unique(modeling_data_surprise$subj_id)


#' Check whether each subject has a sufficient number of trials in each block
#' min = 30
temp <- modeling_data_cntl |>
    group_by(subj_id, block) |>
    summarize(
        n = n_distinct(rt)
    ) |>
    as.data.frame()
hist(temp$n)

#' Check whether each subject has 4 blocks
temp <- modeling_data_cntl |>
    group_by(subj_id) |>
    summarize(
        n = n_distinct(block)
    ) |>
    as.data.frame()
temp


#' Data wrangling

modeling_data_cntl$rt <- modeling_data_cntl$rtTukey / 1000

modeling_data_cntl$subj_id <- factor(modeling_data_cntl$subj_id)
modeling_data_cntl$subject <- as.integer(modeling_data_cntl$subj_id)
sort(unique(modeling_data_cntl$subject))

# new names for required variables
modeling_data_cntl$congruency <- 
  ifelse(modeling_data_cntl$is_congruent_trial == "Congruent", 
         "congruent", "incongruent")
modeling_data_cntl$accuracy <- modeling_data_cntl$correct

mod_data_cntl <- modeling_data_cntl |>
  dplyr::select(subject, block, congruency, accuracy, rt)
mod_data_cntl$subj_name <- NULL

#' Remove NAs
complete_rows <- complete.cases(mod_data_cntl)
dd <- mod_data_cntl[complete_rows, ]



#' select block
# BLOCK <- 1
# modeling_data_cntl_blk <- modeling_data_cntl |>
#   dplyr::filter(block == BLOCK)

# Small test
temp <- dd[dd$subject %in% c(1, 2), ]
dim(temp)

#' set some parameters for the model fit routines

#' during the fit, how many sets of starting parameters should be explored?
n_start_parms <- 50

#' what should the variance across starting parameters be?
var_start_parms <- 20

#' how many trials to simulate during each iteration of the fit routine whilst
#' exploring multiple starting parameters?
n_first_pass <- 1000

#' how many trials to simulate during the final fit routine?
n_final_pass <- 50000

# pass the data to the model's fit functions

#' DSTP Model

# Data
dd <- temp

# Create an empty list to store results
res_dstp_control <- list()

# Get unique subjects and blocks
unique_subjects <- unique(dd$subject)
unique_blocks <- unique(dd$block)

# Loop through subjects and blocks
for (i_sub in unique_subjects) {
  for (i_block in unique_blocks) {
    # Select the data for one subject and one block
    one_subj_block <- dd %>%
      filter(subject == i_sub, block == i_block) %>%
      select(subject, congruency, accuracy, rt)
    
    # Fit the DSTP model
    m <- fitDSTP(data = one_subj_block)
    
    # Save results
    res_dstp_control[[paste("Subject", i_sub, "Block", i_block)]] <- m$bestParameters
    print(paste("Subject", i_sub, "Block", i_block))
  }
}


# First, create an empty dataframe with the desired column names
df <- data.frame(
  Subject = character(),
  Block = character(),
  Parameter1 = numeric(),
  Parameter2 = numeric(),
  Parameter3 = numeric(),
  Parameter4 = numeric(),
  Parameter5 = numeric(),
  Parameter6 = numeric(),
  Parameter7 = numeric(),
  stringsAsFactors = FALSE
)

# Loop through the list and fill the dataframe
for (key in names(res_dstp_control)) {
  subject_block <- unlist(strsplit(key, " "))
  subject <- subject_block[2]
  block <- subject_block[4]
  parameters <- unlist(res_dstp_control[key])
  
  # Add a new row to the dataframe
  new_row <- data.frame(
    Subject = subject,
    Block = block,
    Parameter1 = parameters[1],
    Parameter2 = parameters[2],
    Parameter3 = parameters[3],
    Parameter4 = parameters[4],
    Parameter5 = parameters[5],
    Parameter6 = parameters[6],
    Parameter7 = parameters[7]
  )
  
  df <- rbind(df, new_row)
}

# Reset row names to have consecutive integers
row.names(df) <- NULL

param_names <- 
  c("subject", "block", "A", "C", "mu_ta", "mu_fl", "mu_ss", "mu_rs2", "ter")
colnames(df) <- param_names

# Print the resulting dataframe
print(df)

rio::export(df, "control_DSTP_params.csv")

# eof ----





