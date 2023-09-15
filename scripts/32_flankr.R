#' Modeling of the Flanker Task Performance in the Surprise Experiment
#'
#' File: 31_flankr.R
#'
#' See: James A. Grange (2015). flankr: An R package implementing computational
#' models of attentional selectivity
#' 
#' devtools::install_github("JimGrange/flankr")
#'
#' This script was last modified on "Sat Sep 23 09:32:33 2017"


library("flankr")
library("tidyverse")
library("forcats")
library("car")


df1 <- read.csv(here("data", "processed", "control_2021b.csv"))
df1$experiment <- "control"
df1$experiment <- factor(df1$experiment)

df2 <- read.csv(here("data", "processed", "surprise_2021b.csv"))
df2$experiment <- "surprise"
df2$subj_id <- df2$subj_id + max(df1$subj_id)
df2$experiment <- factor(df2$experiment)

df <- bind_rows(list(df1, df2))

df %>% glimpse()

length(unique(df$subj_id))

length(unique(df[df$experiment == "control", ]$subj_id))
length(unique(df[df$experiment == "surprise", ]$subj_id))

df$rts <- df$rt / 1000

# exclude the first trial after the clip because it behaves differently than the three
# following trials
no_first_trial_dat <- df %>% 
  dplyr::filter(trials_after_clip != 0 & trials_after_clip != 4)


dat <- with(no_first_trial_dat,
  data.frame(subject = subj_id, condition = experiment,
             congruency = is_congruent_trial, accuracy = correct, rt = rts))

dat <- dat %>%
  mutate(congruency = fct_recode(
    congruency,
    "incongruent" = "Incongruent",
    "congruent" = "Congruent"
    )
)

# remove extreme RTs
dat <- dat[dat$rt > 0.2 & dat$rt < 1.5, ]
length(unique(dat$subject))

dd <- dat %>% 
  dplyr::filter(condition == "surprise" & !is.na(condition))  # control

n_sub <- length(unique(dd$subject))
n_sub

dd <- transform(dd, id=match(subject, unique(subject)))
sort(unique(dd$id))
length(unique(dd$id))

res_surprise <- list()

for (i_sub in 1:n_sub) {
  
  # select the data of one subject
  one_subj <- dd %>% 
    dplyr::filter(
      id == i_sub
    )
  # fit the model
  m <- fitDSTP(data = one_subj)
  # save results
  res_surprise[[i_sub]] <- m$bestParameters
  print(i_sub)
  
}

param_names <- c("A", "C", "mu_ta", "mu_fl", "mu_ss", "mu_rs2", "ter")

# surpr_exp_df <- do.call(rbind.data.frame, res)
# colnames(surpr_exp_df) <- param_names

# params_control <- data.frame(t(sapply(res, c)))
# colnames(params_control) <- param_names
# saveRDS(params_control, "dstp_control.Rds")

params_surprise <- data.frame(t(sapply(res_surprise, c)))
colnames(params_surprise) <- param_names
saveRDS(params_surprise, "dstp_surprise.Rds")

saveRDS(params, "dstp_params.Rds")







boxplot(params_surprise$A)
boxplot(params_surprise$C)
boxplot(params_surprise$mu_ta)
boxplot(params_surprise$mu_fl)
boxplot(params_surprise$mu_ss)
boxplot(params_surprise$mu_rs2)
boxplot(params_surprise$ter)

# ggplot(
#   data=params_control, aes(x=block, y=ce)) +
#   facet_wrap(~ exp) +
#   geom_violin(aes(fill=exp, color=exp)) +
#   geom_boxplot(width=.4, outlier.shape=NA) +
#   geom_hline(yintercept=0, linetype="dashed", color = "gray") +
#   ylim(-250, 250)


dim(params_control)

params_control$exp <- "control"
params_surprise$exp <- "surprise"

param_names <- c("A", "C", "mu_ta", "mu_fl", "mu_ss", "mu_rs2", "ter")

params <- rbind(params_control, params_surprise)







# df$rt_grp <- cut(df$rts,
#      quantile(df$rts, c(0, 1/5, 2/5, 3/5, 4/5, 1)),
#      include.lowest = TRUE,
#      labels = c('1', '2', '3', '4', '5'))
# 
# out <- df %>%
#   group_by(is_congruent_trial, rt_grp) %>%
#   summarise(
#     p = mean(correct, na.rm = TRUE)
#   )
# 
# out_rt <- df %>%
#   group_by(is_congruent_trial, rt_grp) %>%
#   summarise(
#     mrt = mean(rts, na.rm = TRUE)
#   )
# 
# out$mrt <- out_rt$mrt
# 
# ggplot(aes(x = mrt, y = p,
#            group = is_congruent_trial,
#            color = is_congruent_trial),
#       data = out) +
#       geom_point(size = 4) +
#       geom_line() +
#       coord_cartesian(ylim=c(0.5, 1)) +
#       labs(colour = "Stimulus Congruency",
#            x = "Response Time (s)",
#            y = "Accuracy") 
