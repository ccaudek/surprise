# Script name: 10_tidying_data.R
# Project: Surprise
# Script purpose: pre-processing
# @author: Corrado Caudek <corrado.caudek@unifi.it>
# Date Created: Sun May  6 06:06:34 2018
# Last Modified Date: Tue Sep 20 06:47:39 2022
# 
# Notes: 

library("here")

source(here("libraries", "helpers.R"))

length(unique(df$subj_name))


# filter wrong RTs 
summary(df$resp)
summary(df$rt)

df$rt <- ifelse(df$resp == "NORESP" | df$rt <= 150 | df$rt > 2500, NA, df$rt)
hist(df$rt)

# Transform and create variables
df$is_congruent_trial <- ifelse(df$target_or == df$flanker_or, 1, 0)
df$is_congruent_trial <- factor(df$is_congruent_trial)
df <- df %>%
  mutate(is_congruent_trial = fct_recode(
    is_congruent_trial,
    "No"  = "0",
    "Yes" = "1"
  )
)

df$is_clip_trial <- factor(df$is_clip_trial)
df <- df %>%
  mutate(is_clip_trial = fct_recode(is_clip_trial,
    "No"  = "0",
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
  group_by(subj_name) %>%
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
  group_by(subj_name) %>%
  summarise(
    n = n()
  )

if (0) ggplot(data=bysubj_n, aes(bysubj_n$n)) + 
  geom_histogram()

# index trial number for each subject
temp <- bysubj_n$n
ntrial <- NULL
for (i in 1:length(temp)) {
  v <- 1:temp[i]
  ntrial <- append(ntrial, v)
}
rm(temp)

df$ntrial <- ntrial

if (exp_flag == "CONTROL") {
  df$is_surprise_clip <- factor(df$is_surprise_clip)
  df$is_surprise_clip <- "No"
} else {
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
df$is_surprise_clip <- factor(df$is_surprise_clip)


# The coding is wrong for the trials that follow the trial immediately 
# after the presentation of the video clip! It will be corrected below.
# for (i in 1:nrow(df)) {
#   if (df$trials_after_clip[i] == 1) {
#     df$is_surprise_clip[i] <- df$is_surprise_clip[i - 1]
#   } else if (df$trials_after_clip[i] == 2) {
#     df$is_surprise_clip[i] <- df$is_surprise_clip[i - 2]
#   } else if (df$trials_after_clip[i] == 3) {
#     df$is_surprise_clip[i] <- df$is_surprise_clip[i - 3]
#   } else if (df$trials_after_clip[i] == 4) {
#     df$is_surprise_clip[i] <- df$is_surprise_clip[i - 4]
#   } 
# }


# # Remove trials with no responses 
# df <- df %>%
#   dplyr::filter(response != "NORESP")
# df$response <- fct_drop(df$response)


# Outliers handling -------------------------------------------------------

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
df1 <- df
df1$trt <- ifelse(df1$rt < 150 | df1$rt > 2500, NA, df1$rt)  # 200, 2000

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

df1 <- transform(df1, rtTukey = ave(trt, subj_id, FUN = trimIqr))

if (0) plot_rt_bysub(df, 16)

out <- df1 %>% 
  group_by(subj_name) %>% 
  summarise(
    p = mean(correct, na.rm = TRUE)
  )
hist(out$p)

# accuracty too low; subjects with p < 0.3 switched the response categories.
bad_subjs <- out %>%
  dplyr::filter(p < 0.75 & p > 0.3)

# remove this subject
df2 <- df1[df1$subj_name != "el_gr_12_09_94_912_f", ]

out <- df2 %>% 
  group_by(subj_name) %>% 
  summarise(
    p = mean(correct, na.rm = TRUE)
  )
hist(out$p)

reversed_resp_subjs <- out %>% 
  dplyr::filter(p < 0.3)

temp <- df2 %>% 
  dplyr::filter(subj_name != "69iv18_m19")

temp1 <- df2 %>% 
  dplyr::filter(subj_name == "69iv18_m19")

temp1$correct <- ifelse(temp1$correct == 0, 1, 0)

df2 <- rbind(temp, temp1)


out <- df2 %>% 
  group_by(subj_name) %>% 
  summarise(
    p = mean(correct, na.rm = TRUE)
  )
hist(out$p)

bysubj_acc <- df2 %>%
  group_by(subj_name) %>%
  summarise(
    p = mean(correct)
  )
data.frame(bysubj_acc)


good_subj <- bysubj_acc %>%
  dplyr::filter(p > 0.75)


# Only good subjects
mydata <- subset(
  df2,
  subj_name %in% good_subj$subj_name
)

out <- mydata %>% 
  group_by(subj_name) %>% 
  summarise(
    p = mean(correct, na.rm = TRUE)
  )
hist(out$p)

mydata %>% 
  summarise(
    n = n_distinct(subj_name)
  )


mydata <- mydata %>%
  mutate(is_congruent_trial = fct_recode(is_congruent_trial,
    "Incongruent" = "No",
    "Congruent" = "Yes"
  )
)

length(unique(mydata$subj_name))


if (exp_flag == "SURPRISE") {
  mydata <- mydata %>%
    mutate(
      is_surprise_clip = fct_recode(
        is_surprise_clip,
        "No Surprise" = "No",
        "Surprise" = "Yes"
      )
    )
}


# Write CSV in R
if (exp_flag == "CONTROL") {
  write.csv(mydata, here("data", "processed", "control_2021b.csv"), row.names = FALSE)
}

if (exp_flag == "SURPRISE") {
  write.csv(mydata, here("data", "processed", "surprise_2021b.csv"), row.names = FALSE)
}


# Exit message
n_subj <- length(unique(mydata$subj_name))
cat("\nPre-processing done!\nNumber of subjects: ", n_subj)

