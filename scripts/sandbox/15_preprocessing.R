#' Surprise Experiment
#' Preparation for RT measures: procedure of Pe, Vandekerckhove, Kuppens (2013).
#'
#' File: 15_preprocessing.R
#'
#' This file was last modified on "Mon May  7 15:25:46 2018"




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

