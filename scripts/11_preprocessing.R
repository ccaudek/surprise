#' Surprise Experiment
#' 
#' Merge the csv files with the experimental and control conditions
#'
#' File: 11_preprocessing.R
#'
#' This file was last modified on ""Tue Jun 18 13:42:47 2019"


library("here")
suppressPackageStartupMessages(library("tidyverse")) 

exp_df  <- read_csv(here("data", "processed", "surprise_2021b.csv"))
cntr_df <- read_csv(here("data", "processed", "control_2021b.csv"))

names(exp_df)
names(cntr_df)

# baa <- as.character(cntr_df$video_id)
# cntr_df$video_id2 <- sapply(strsplit(baa, "/stim/"), "[", 2)

# temp_exp <- data.frame(
#   subj_name = exp_df$subj_name, trial = exp_df$trial, target_or = exp_df$target_or,        
#   flanker_or = exp_df$flanker_or, resp = exp_df$resp, correct = exp_df$correct, 
#   movie_id = exp_df$movie_id2, rt = exp_df$rt,  is_surprise_clip = exp_df$is_surprise_clip,  
#   is_clip_trial = exp_df$is_clip_trial, block = exp_df$block,
#   is_congruent_trial = exp_df$is_congruent_trial, trials_after_clip = exp_df$trials_after_clip, 
#   is_first_trial = exp_df$is_first_trial, ntrial = exp_df$ntrial, experiment = "surprise"
# )
exp_df$ID <- as.numeric(factor(exp_df$subj_name))
exp_df$experiment <- "surprise"

# temp_contr <- data.frame(
#   subj_name = cntr_df$subj_name, trial = cntr_df$trial, target_or = cntr_df$target_or,        
#   flanker_or = cntr_df$flanker_or, resp = cntr_df$resp, correct = cntr_df$correct, 
#   movie_id = cntr_df$video_id2, rt = cntr_df$rt,  is_surprise_clip = cntr_df$is_surprise_clip,  
#   is_clip_trial = cntr_df$is_clip_trial, block = cntr_df$block,
#   is_congruent_trial = cntr_df$is_congruent_trial, trials_after_clip = cntr_df$trials_after_clip, 
#   is_first_trial = cntr_df$is_first_trial, ntrial = cntr_df$ntrial, experiment = "control"
# )

cntr_df$ID <- max(exp_df$ID) + as.numeric(factor(cntr_df$subj_name))
cntr_df$experiment <- "control"

both_exps_df <- rbind(exp_df, cntr_df)
table(both_exps_df$ID, both_exps_df$experiment)

# The correct subject identifier is ID.
length(unique(both_exps_df$ID))
length(unique(both_exps_df$subj_id))

both_exps_df$subj_id <- NULL
both_exps_df <- both_exps_df |> 
  dplyr::rename(subj_id = ID)

length(unique(both_exps_df$subj_id))

write.csv(
  both_exps_df, 
  here("data", "processed", "surprise_control_exps.csv"), 
  row.names = FALSE
)


# eof ----








# Check whether there are outlying participants ---------------------------

if (0) {
  
  bysubraw <- mydata %>%
    group_by(subj_name, is_congruent_trial) %>%
    summarise(
      mrt = mean(rtTukey, trim = 0.1, na.rm = TRUE)
    )
  
  bysub_con <- dplyr::filter(bysubraw, is_congruent_trial == "Congruent")
  bysub_inc <- dplyr::filter(bysubraw, is_congruent_trial == "Incongruent")
  
  bysub_df <- data.frame(subj_name = bysub_con$subj_name)
  
  bysub_df$interference_eff <- 
    filter(bysubraw, is_congruent_trial == "Incongruent")$mrt - 
    filter(bysubraw, is_congruent_trial == "Congruent")$mrt
  
  hist(bysub_df$interference_eff)
  
  summary(bysub_df$interference_eff)
  #     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
  # -184.677  -30.657    6.733    2.882   50.141  128.259 
  
  temp <- mydata %>%
    group_by(subj_name) %>%
    summarise(
      acc = mean(correct, na.rm = TRUE),
      med_rt = median(rtTukey, na.rm = TRUE)
    )
  
  dd <- merge(bysub_df, temp, by = "subj_name")
  
  bysub_df$acc <- temp$acc
  bysub_df$med_rt <- temp$med_rt
  bysub_df
  
  mean(dd$med_rt) + c(-1, 1) * 2.5 * sd(dd$med_rt)
  
}
# Nothing strange! Five subjects are too slow and have almost perfect accuracy.  
# They could be deleted, because they did not perform the task "as quickly as possible".

