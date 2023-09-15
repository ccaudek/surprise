#' Surprise Experiment
#' Pre-processing.
#'
#' File: 10_tidying_data.R
#'
#' This file was last modified on "Sun May  6 06:06:34 2018"



source("../lib/helpers.R")


df$rt <- ifelse(df$resp == "NORESP", NA, df$rt)
df$rt <- ifelse(df$rt < 200, NA, df$rt)



# Change column names to snake_case 
old_cols <- colnames(df)
new_cols <- to_snake_case(names(df))
colnames(df) <- new_cols


# Add block number 
create_block_num <- function(mydf) {
  blk <- rep(NA, nrow(mydf))
  for (i in 1:nrow(mydf)) {
    if (i == 1) {
      blk[i] <- 1
    } else {
      if (mydf$subject_number[i - 1] != mydf$subject_number[i])
        blk[i] <- blk[i - 1] + 1
      else
        blk[i] <- blk[i - 1] 
    }
  }
  blk
}

block_names_df <- df %>%
  group_by(subject_name) %>%
  do(data.frame(block = create_block_num(.)))

df$block <- block_names_df$block



# Transform and create variables ------------------------------------------


df$id <- as.numeric(factor(df$subject_name))

df$is_congruent_trial <- ifelse(
  df$target_orientation == df$flankers_orientation, 1, 0
)
df$is_congruent_trial <- factor(df$is_congruent_trial)
df <- df %>%
  mutate(is_congruent_trial = fct_recode(
    is_congruent_trial,
    "No" = "0",
    "Yes" = "1"
  ))

df <- rename(df, trial_within_blk = i_trial)
df <- rename(df, video_id = magic_shows)

df$is_clip_trial <- factor(df$is_clip_trial)
df <- df %>%
  mutate(is_clip_trial = fct_recode(is_clip_trial,
    "No" = "0",
    "Yes" = "1"
  ))
df$is_clip_trial <- factor(df$is_clip_trial, levels = c("Yes", "No"))


# How many trials are there after the video clip? 
get_trial_n_after_clip <- function(mydf) {
  
  trials_after <- rep(NA, nrow(mydf))
  for (i in 1:nrow(mydf)) {
    if (mydf$is_clip_trial[i] == "Yes") {
      trials_after[i] <- 0
    } else if (mydf$is_clip_trial[i - 1] == "Yes") {
      trials_after[i] <- 1
    } else if (mydf$is_clip_trial[i - 2] == "Yes") {
      trials_after[i] <- 2
    } else if (mydf$is_clip_trial[i - 3] == "Yes") {
      trials_after[i] <- 3
    } else if (mydf$is_clip_trial[i - 4] == "Yes") {
      trials_after[i] <- 4
    } else {
      trials_after[i] <- 999
    }
  }
  trials_after
  
}

trials_after_df <- df %>%
  group_by(subject_name) %>%
  do(data.frame(trials_after = get_trial_n_after_clip(.)))

df$trials_after_clip <- trials_after_df$trials_after

# Remember: only in the first 3 trials of each block, the third trial
# is also a clip trial! From trial 4 on, the third trial is not a clip 
# trial. 
# df$trials_after_clip[1:16]
# [1] 0 1 2 0 1 2 3 0 1 2 3 0 1 2 3 0
# The last trial of the block is a non-clip trial, even if is preceeded
# by 3 no-clip trials.

# in about the 1% of the cases there are 5 trials after the clip
# table(df$trials_after_clip)
# 593 / nrow(df)


# First trial after the video clip 
df$is_first_trial <- ifelse(df$trials_after_clip == 0, 1, 0)
df$is_first_trial <- factor(df$is_first_trial)
df <- df %>%
  mutate(is_first_trial = fct_recode(is_first_trial,
    "No" = "0",
    "Yes" = "1"
  ))
df$is_first_trial <- factor(df$is_first_trial, levels = c("Yes", "No"))
# table(df$is_first_trial)


# Get trial number for subject
bysubj_n <- df %>%
  group_by(subject_name) %>%
  summarise(
    n = n()
  )

if (0) ggplot(data=bysubj_n, aes(bysubj_n$n)) + 
  geom_histogram()

temp <- bysubj_n$n
ntrial <- NULL
for (i in 1:length(temp)) {
  v <- 1:temp[i]
  ntrial <- append(ntrial, v)
}
rm(temp)

df$ntrial <- ntrial

if (0) {
  bysub_n <- df %>% 
    group_by(subject_name) %>% 
    summarise(
      n = n()
    )
  data.frame(bysub_n)
}

# Subject 01pm23_f25 has only 72 trials.

# Subject 24sm21_f21 has 1122 trials!
# Remove last blocks of subjects 24sm21_f21

if (exp_flag == "SURPRISE") {
  temp1 <- df %>%
    dplyr::filter(subject_name == "24sm21_f21")
  temp2 <- temp1[1:345, ]
  
  temp3 <- df %>%
    dplyr::filter(subject_name != "24sm21_f21")
  
  df <- rbind(temp3, temp2)
  rm(temp1, temp2, temp3)
}



# Get trial number after the clip
if (exp_flag == "CONTROL") {
  df <- rename(df, is_surprise_clip = surprise_clip)
  df$is_surprise_clip <- factor(df$is_surprise_clip)
  df$is_surprise_clip <- "No"
} else {
  df <- rename(df, is_surprise_clip = surprise_clip)
  df$is_surprise_clip <- factor(df$is_surprise_clip)
  df <- df %>%
    mutate(is_surprise_clip = fct_recode(is_surprise_clip,
      "No" = "0",
      "Yes" = "1"
    ))
  df$is_surprise_clip <- factor(
    df$is_surprise_clip,
    levels = c("Yes", "No")
  )
}


# Remove 345 NAs ----------------------------------------------------------

# if (exp_flag == "CONTROL") {
#   
#   n_subjects <- length(unique(df$subject_name))
#   
#   out <- df %>%
#     group_by(subject_name) %>%
#     summarise(
#       y = mean(rt)
#     )
#   
#   good <- out$subject_name[1:n_subjects]
#   df <- df[df$subject_name %in% good, ]
#   
# }


# The coding is wrong for the trials that follow the trial immediately 
# after the presentation of the video clip! It will be corrected below.

for (i in 1:nrow(df)) {
  if (df$trials_after_clip[i] == 1) {
    df$is_surprise_clip[i] <- df$is_surprise_clip[i - 1]
  } else if (df$trials_after_clip[i] == 2) {
    df$is_surprise_clip[i] <- df$is_surprise_clip[i - 2]
  } else if (df$trials_after_clip[i] == 3) {
    df$is_surprise_clip[i] <- df$is_surprise_clip[i - 3]
  } else if (df$trials_after_clip[i] == 4) {
    df$is_surprise_clip[i] <- df$is_surprise_clip[i - 4]
  } 
}


# # Remove trials with no responses 
# df <- df %>%
#   dplyr::filter(response != "NORESP")
# df$response <- fct_drop(df$response)

df$rt <- ifelse(df$response == "NORESP", 0, df$rt)


# Outliers handling -------------------------------------------------------


if (0) plot_rt_bysub(clean_df, 121)

# RT and accuracy scores for interference and facilitation effects were calculated separately
# for both thepositive-target condition (PTC) and the negative-target condition(NTC). For each
# flanker condition of every target valence, we applied the data preparation procedure
# of Friedman and Miyake(2004)for RT measures. Only correct trials were included in the RT analyses
# (M = 92%, SD = 5%). Upper and lower criteria were determined through visual inspection of the
# overall RT distributions; values under 200 ms and over 1500 ms were eliminated. Foreach
# participant, RTs more than 3SDfrom the participant’s meanfor each condition were replaced with
# values that were 3 SD fromthe participant’s mean for that condition. Between-subjects RT
# distributions were then examined for each condition, and scoresabove or below 3 SD from
# the group mean were replaced withvalues that were plus or minus 3SDfrom the mean, respectively.

# NAs for RTs with responses too long or too short
df$trt <- ifelse(df$rt < 200 | df$rt > 2000, NA, df$rt)

# Responses outside the interval of the individual median +/- 3 × the
# interquartile range are considered as outliers (Tukey 1977)
trimIqr <- function(x) {
  nIQR <- 3
  iqr <- quantile(x, 3 / 4, na.rm = TRUE) - quantile(x, 1 / 4, na.rm = TRUE)
  x1 <- ifelse(
    (x > (quantile(x, 2 / 4, na.rm = TRUE) + 3 * iqr)) |
      (x < (quantile(x, 2 / 4, na.rm = TRUE) - 3 * iqr)),
    NA, x
  )
  return(x1)
}

df <- transform(df, rtTukey = ave(trt, id, FUN = trimIqr))

# qqPlot(log(df$trt))
# qqPlot(1/df$trt)
# qqPlot(df$trt)
# The log transformation works best


# Remove noisy subjects ---------------------------------------------------

if (0) {
  bysubj_acc <- df %>%
    group_by(subject_name) %>%
    summarise(
      mean(correct)
    )
  data.frame(bysubj_acc)
}


# Response key reversal
for (i in 1:nrow(df)) {
  if (df$subject_name[i] %in% c("89lv24_f20", "74np03_f25", 
                                "30ts30_f19", "69iv18_m19") ) {
    df$correct[i] <- ifelse(df$correct[i] == 0, 1, 0) 
  } 
}

if (0) data.frame(bysubj_acc)

# subject with accuracy smaller than 80%
bad <- c(
  #"89lv24_f20",
  #"74np03_f25",
  #"30ts30_f19",
  "48gp30_m19",
  "88fs29_m20",
  "65vz30_f20",
  "45ln03_f20",
  "40mz27_f19",
  "76fr07_f19",
  "19il25_f19",
  "67cp27_f19"
)


# Only good subjects
mydata <- subset(
  df,
  !subject_name %in% bad
)

mydata <- mydata %>%
  mutate(is_congruent_trial = fct_recode(is_congruent_trial,
    "Incongruent" = "No",
    "Congruent" = "Yes"
  ))

mydata <- mydata %>%
  mutate(
    is_surprise_clip = fct_recode(
      is_surprise_clip,
      "No Surprise" = "No",
      "Surprise" = "Yes"
    )
  )


# Remove unnecessary columns
mydata$unshuffled_trial <- NULL
mydata$onscreen_stim_duration <- NULL



# Write CSV in R ----------------------------------------------------------

# write.csv(mydata, file = "surprise_cntr.csv", row.names = FALSE)
# write.csv(mydata, file = "surprise_expr.csv", row.names = FALSE)



# Exit message ------------------------------------------------------------


n_subj <- length(unique(mydata$subject_name))
cat("\nPre-processing done!\nNumber of subjects: ", n_subj)

