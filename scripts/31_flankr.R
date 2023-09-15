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
  dplyr::filter(condition == "control")  # control

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





params_control <- data.frame(t(sapply(res, c)))
colnames(params_control) <- param_names
saveRDS(params_control, "dstp_control.Rds")


boxplot(params_control$A)
boxplot(params_control$C)
boxplot(params_control$mu_ta)
boxplot(params_control$mu_fl)
boxplot(params_control$mu_ss)
boxplot(params_control$mu_rs2)
boxplot(params_control$ter)

# ggplot(
#   data=params_control, aes(x=block, y=ce)) +
#   facet_wrap(~ exp) +
#   geom_violin(aes(fill=exp, color=exp)) +
#   geom_boxplot(width=.4, outlier.shape=NA) +
#   geom_hline(yintercept=0, linetype="dashed", color = "gray") +
#   ylim(-250, 250)




tot_df <- rbind(control_df, surpr_exp_df, nosurpr_exp_df)
tot_df$experiment <- c(rep("control", 10), rep("exp_surprise", 20), rep("exp_nosurprise", 20))

tot_df$exp <- ifelse(tot_df$experiment == "control", 0, 1)



fm <- lm(cbind(mu_ta, mu_fl, mu_ss) ~ exp, data = tot_df)
Anova(fm)
summary(fm)

fm <- lm(mu_fl ~ exp, data = tot_df)
Anova(fm)
summary(fm)

fm <- lm(mu_ss ~ exp, data = tot_df)
Anova(fm)
summary(fm)

fm <- lm(C ~ exp, data = tot_df)
Anova(fm)
summary(fm)


t.test(mu_fl ~ experiment, data = tot_df)
t.test(mu_ss ~ experiment, data = tot_df)
t.test(ter ~ experiment, data = tot_df)
t.test(A ~ experiment, data = tot_df)
t.test(C ~ experiment, data = tot_df)








# select the data of one subject
one_subj_dat <- dat %>% 
  dplyr::filter(
    subject == 2
  )






m1 <- fitDSTP(data = df_surprise)

m2 <- fitDSTP(data = df_nosurprise)

plot1 <- plotFitDSTP(modelFit = m1, data = df_surprise)
plot2 <- plotFitDSTP(modelFit = m2, data = df_nosurprise)



plot <- plotFitDSTP(modelFit = fit5, data = supr_df)




fit3 <- fitDSTP(data = dd, conditionName = "No Surprise")
# control condition (no surprise)
# (1) mu_ta (2) mu_fl (3) A (4) mu_ss (5) C (6) mu_rs2 (7) t_er
# mu = rate or drift of a given diffusion process and condition; 
# ta = target; 
# fl = flanker; 
# A = criterion for response selection; 
# SS = stimulus selection process; 
# C = criterion for stimulus selection; 
# RS2 = response selection process in Phase 2; 
# t_er = nondecisional parameter that represents time used for stimulus encoding, 
# response execution, and so forth.
# $bestParameters
# [1] 0.130 0.111 0.000 0.049 0.370 1.192 0.190
# 
# $g2
# [1] 46.63694
# 
# $bBIC
# [1] 1977.794

# surprise condition in the 2-video condition
# $bestParameters
# [1] 0.160 0.131 0.005 0.053 0.281 1.369 0.119
# 
# $g2
# [1] 24.97205
# 
# $bBIC
# [1] 1938.025

# no surprise condition in the 2-video condition
# $bestParameters
# [1] 0.215 0.142 0.096 0.131 0.300 1.373 0.096
# 
# $g2
# [1] 118.0751
# 
# $bBIC
# [1] 2045.489

# Height of response selection boundary
# Height of stimulus selection boundary
# Drift rate for central target during response selection stage 1 
# Drift rate for flankers during response selection stage 1 
# Drift rate for stimulus selection
# Drift rate for stage 2 of response selection
# Non-decision time


# $bestParameters
# [1] 0.111 0.083 0.016 0.073 0.338 1.272 0.261
#
# $g2
# [1] 39.5643
#
# $bBIC
# [1] 1959.925

# surprise NO
# $bestParameters
# [1] 0.133 0.109 0.005 0.065 0.322 1.367 0.218
#
# $g2
# [1] 48.37825
#
# $bBIC
# [1] 1985.64

plot <- plotFitDSTP(modelFit = fit2, data = dd)

#' Conclusion: The diffusion model is adequate for the control condition (no
#' surprise in any trial), but not for the experimental conditions, in which
#' in half of the trials the videos show a surprising event. The data will not
#' be analyzed with a DM approach, because the process is more complex than a
#' single decision task.


nBootstraps <- 1000

mean(tapply(mydata$rtTukey, mydata$subject_name, length))

nTrials <- 140

bestParameters <- c(0.130, 0.111, 0.000, 0.049, 0.370, 1.192, 0.190)

bootData <- matrix(0, nrow = nBootstraps, ncol = length(bestParameters))

for (i in 1:nBootstraps) {
  
  simData <- simulateDSTP(parms = bestParameters, nTrials = nTrials)
  
  fit <- fitDSTP(data = simData, 
                 parms = bestParameters,
                 nTrials = 10000,
                 multipleSubjects = FALSE)
  
  bootData[i, ] <- fit$bestParameters
}

apply(bootData, 2, sd) / sqrt(dim(bootData)[1])







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
