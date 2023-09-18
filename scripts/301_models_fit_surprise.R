#' ---
#' title: "Models' fit: surprise condition"
#' output:
#'   pdf_document:
#'     keep_tex: true
#' ---


library("here")

suppressPackageStartupMessages({
    library("tidyverse")
    library("brms")
    library("cmdstanr")
    library("reshape")
    library("devtools")
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


#' Data wrangling for models' fit
#' Select experiment

modeling_data <- data_tidy |>
    dplyr::filter(experiment == "surprise") # surprise

# unique(modeling_data_surprise$subj_id)


#' Data wrangling

modeling_data$rt <- modeling_data$rt / 1000

modeling_data$subj_id <- factor(modeling_data$subj_id)
modeling_data$subject <- as.integer(modeling_data$subj_id)

# new names for required variables
modeling_data$congruency <-
  ifelse(modeling_data$is_congruent_trial == "Congruent", 
         "congruent", "incongruent")
modeling_data$accuracy <- modeling_data$correct

input_data <- modeling_data |>
  dplyr::select(subject, block, congruency, accuracy, rt) |> 
  ungroup()


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

# Sorted data
dd <- input_data %>%
  arrange(subject, block)

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





