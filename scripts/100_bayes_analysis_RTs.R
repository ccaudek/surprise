#' Surprise Experiment
#' Exploratory analyses.
#'
#' File: 100_bayes_analysis_RTs.R
#'
#' This script was last modified on "Thu Dec 12 11:38:51 2019"


library("here")
suppressPackageStartupMessages(library("tidyverse")) 
library("brms")
# library("cowplot")
# library("reshape")
library("devtools")
suppressPackageStartupMessages(library("rethinking")) 

source(here("libraries", "fnct_surprise.R"))
source(here("libraries", "helpers.R"))

options(mc.cores = parallel::detectCores())

## use sum-to-zero contrasts
options(contrasts=c("contr.sum", "contr.poly"))



# Read data ---------------------------------------------------------------


d1 <- read.csv(here("data", "processed", "surprise_control_exps.csv"))

d1 %>% glimpse()

# Number of subjects per esperiment.
d1 %>%
  group_by(experiment) %>% 
  summarise(
    n_distinct(subj_name)
  )



# Data wrangling ----------------------------------------------------------


# Remove blocks with less than 10 trials 
out <- d1 %>% 
  group_by(subj_name, block) %>% 
  summarise(
    ntrials_per_block = n()
  )
d2 <- left_join(d1, out)
d3 <- d2[d2$ntrials_per_block > 10, ]

# add number of trials for each subject
d3 <- d3 %>% 
  group_by(subj_name) %>% 
  mutate(
    bysub_n = length(rtTukey)
  )


# remove 5th trials in the last sequence of trials and blocks after the fourth
d4 <- d3 %>%
  dplyr::filter(
    trials_after_clip != "4" & block < 5
  )

d4$experiment <- factor(d4$experiment)
d4$subj_name <- factor(d4$subj_name)
d4$resp <- factor(d4$resp)
d4$movie_id <- factor(d4$movie_id)
d4$is_surprise_clip <- factor(d4$is_surprise_clip)
d4$is_clip_trial <- factor(d4$is_clip_trial)



# Contrasts ---------------------------------------------------------------


contrasts(d4$experiment) <- contr.sum(2)
contrasts(d4$experiment)

contrasts(d4$is_surprise_clip) <- contr.sum(2)
contrasts(d4$is_surprise_clip)

contrasts(d4$is_congruent_trial) <- contr.sum(2)
contrasts(d4$is_congruent_trial)



# Correct responses only --------------------------------------------------


# Overall proportion of errors
table(d4$correct)[1] / sum(table(d4$correct))

# model er ~ rt to eyeball where performance becomes better than chance

library("mgcv")

foo <- data.frame(
  er = abs(1 - d4$correct),
  rt = d4$rt
)

foo <- foo[complete.cases(foo), ]

y = foo$er
x = foo$rt

fit = mgcv::bam(
  formula = y ~ s(x, bs='ts', k=10), 
  family = binomial
)

preds = data.frame(
  x = seq(min(x), max(x), length.out = 1e3)
)

preds$value = plogis(predict(fit, newdata = preds))

plot(preds, type = 'l')
abline(h = 0.5)
abline(v = 250)

d4$RT <- d4$rt / 1000

d4$toss <- d4$RT < 0.250
d4$toss <- ifelse(is.na(d4$RT), TRUE, d4$toss)
d4$toss <- ifelse(is.na(d4$toss), TRUE, d4$toss)
d4$toss <- ifelse(d4$resp == "NORESP", TRUE, d4$toss)
d4$toss <- ifelse(d4$RT > 2.5, TRUE, d4$toss)

d5 <- d4[!d4$toss, ]
ds <- d5[!is.na(d5$subj_name), ]

rm(d1, d2, d3, d4, d5)


ds$error <- ifelse(ds$correct == 1, 0, 1)
ds$id <- as.numeric(factor(ds$subj_name))

rand_sample_subjs <- sample(unique(ds$id), 9)
temp <- ds[ds$id %in% rand_sample_subjs, ]

# viz of individual Ss' rt distribtions
ggplot(
  data = temp[temp$correct == 1, ], 
  mapping = aes(
    x = RT, 
    # fill = error, 
    y = ..count..)
  ) +
  facet_wrap(
    ~ id, 
    scale = 'free_y') +
  geom_histogram(
    position = 'identity', 
    alpha = .5
   # , binwidth = 50
   ) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.title = element_blank(), 
    axis.text = element_blank(), 
    axis.ticks = element_blank()
  )



# melsm: mixed effects location scale model -------------------------------


# me: mixed effects model
fit_1 <- brm(
  bf(
    RT ~ 1 + (1 | p | ID),
    sigma ~ 1 + (1 | p | ID)
  ),
  data = dt_cor %>% filter(is_congruent_trial == "Incongruent", experiment == "control"),
  cores = parallel::detectCores(),
  init_r = 0.05,
  chains = 2,
  iter = 2000,
  warmup = 1000
)


# remove scale model
fit_2 <- brm(
  bf(RT ~ 1 + (1 | p | ID)), 
  data = dt_cor %>% filter(is_congruent_trial == "Incongruent", experiment == "control"),
  cores = parallel::detectCores(),
  init_r = 0.05,
  chains = 2,
  iter = 2000,
  warmup = 1000
)

# Compare with WAIC.
incong_waic <- waic(fit_1, fit_2)





ds_cor <- ds %>% 
  dplyr::filter(correct == 1)




form <- brmsformula(
  RT ~ is_congruent_trial + (is_congruent_trial | c | ID), 
  sigma  ~ is_congruent_trial + (is_congruent_trial | c | ID)
)



priors <- c(
  # correlations
  set_prior("lkj(2)", class = "cor"),
  # location random SD
  # set_prior("normal(0,0.25)", 
  #           class = "sd", 
  #           coef = "congruencyincongruent", 
  #           group = "ID"),
  # location random SD
  set_prior("normal(0, 0.25)",
            class = "sd", 
            coef = "Intercept",
            group = "ID"),
  # scale random SD
  # set_prior("normal(0,1)", 
  #           class = "sd",
  #           coef = "congruencyincongruent", 
  #           group = "ID", 
  #           dpar = "sigma"), 
  # scale random SD
  set_prior("normal(0, 1)", 
            class = "sd", 
            coef = "Intercept", 
            group = "id", 
            dpar = "sigma"), 
  # fixed effect priors
  set_prior("normal(0, 5)", 
            class = "b"),
  set_prior("normal(0, 5)", 
            class = "b", 
            dpar = "sigma"),
  set_prior("normal(0, 5)", 
            class = "Intercept"),
  set_prior("normal(0, 5)", 
            class = "Intercept", 
            dpar = "sigma"))




fit_flank <- brm(
  form, 
  data = ds_cor %>% filter(experiment == "control"), 
  init_r = 0.05, 
  cores = parallel::detectCores(),
  chains = 2, 
  iter = 1000, 
  warmup = 500,  
  prior = priors
)



dt_cor$trials_after <- factor(dt_cor$trials_after_clip)



# Experimental effect -----------------------------------------------------



ds_cor <- ds %>% 
  dplyr::filter(correct == 1)

ds_cor$trials_after_clip <- factor(ds_cor$trials_after_clip)


bysubj <- ds_cor %>% 
  group_by(experiment, id, block, trials_after_clip, is_congruent_trial) %>% 
  summarise(
    RT = median(RT, na.rm = TRUE)
  )


ds_cor_surprise <- 
  bysubj[bysubj$experiment == "surprise" & bysubj$block == 1, ]


bf1 <- bf(
  RT ~ is_congruent_trial * trials_after_clip + 
    (is_congruent_trial * trials_after_clip | p | id)
)


priors_1 <- c(
  set_prior("student_t(5, 0, 10.0)", class = "b"),
  set_prior("student_t(5, 0, 10.0)", class = "Intercept"),
  set_prior("cauchy(0, 10.0)", class = "sd")
)



fit1_exp <- brm(
  bf1, 
  data = ds_cor_surprise, #  %>% dplyr::filter(experiment == "surprise"), 
  family = exgaussian(),
  init_r = 0.05, 
  chains = 4,
  cores = parallel::detectCores(),
  iter = 4000,
  warmup = 1000, 
  prior = priors_1
)


pp_check(fit1_exp)

print(fit1_exp, 4)

conditional_effects(
  fit1_exp, 
  "is_congruent_trial:trials_after_clip"
)






ds_cor_control <- 
  bysubj[bysubj$experiment == "control" & bysubj$block == 1, ]


fit1_cnt <- brm(
  bf1, 
  data = ds_cor_control,  
  family = exgaussian(),
  init_r = 0.05, 
  chains = 4,
  cores = parallel::detectCores(),
  iter = 4000,
  warmup = 1000, 
  prior = priors_1
)


pp_check(fit1_cnt)

print(fit1_cnt, 4)

conditional_effects(
  fit1_cnt, 
  "is_congruent_trial:trials_after_clip"
)










plot(pc)[[1]] + ggplot2::ylim(0.50, 1.0)


fit3_surprise <- brm(
  form, 
  data = dt_cor %>% filter(experiment == "surprise", block == 1), 
  family = exgaussian(),
  inits = 0, 
  cores = 4, 
  chains = 2, 
  iter = 1500, 
  warmup = 1000,  
  prior = priors_1
)

pp_check(fit2_surprise)
print(fit3_surprise, 4)
ps <- marginal_effects(
  fit3_surprise, 
  "trials_after:is_congruent_trial"
)
plot(ps)[[1]] + ggplot2::ylim(0.50, 1.0)




#' computing the size of the congruency effect for each
#' subject and trial number after the clip

bysub_df <- mydata %>% 
  group_by(
    subj_name, is_surprise_clip, trials_after_clip, is_congruent_trial
  ) %>% 
  summarise(
    y = mean(rtTukey, na.rm = TRUE, trim = 0.1)
  )

inc <- bysub_df %>% 
  dplyr::filter(is_congruent_trial == "Incongruent")

con <- bysub_df %>% 
  dplyr::filter(is_congruent_trial == "Congruent")

congr_eff_df <- inc
congr_eff_df$is_congruent_trial <- NULL
congr_eff_df$y <- NULL

congr_eff_df$congr_eff <- inc$y - con$y

hist(congr_eff_df$congr_eff)

plot_df <- congr_eff_df %>% 
  group_by(is_surprise_clip, trials_after_clip) %>% 
  summarise(
    y = mean(congr_eff, na.rm = TRUE),
    stderr = sqrt(var(congr_eff) / n())
  )
plot_df$lower <- plot_df$y - plot_df$stderr
plot_df$upper <- plot_df$y + plot_df$stderr

pd <- position_dodge(0.2)

ggplot(plot_df, aes(x=trials_after_clip, y=y)) +
  geom_pointrange(aes(ymin=lower, ymax=upper), lwd=1.05, position=pd) +
  geom_line(position=pd, lwd=1.05) +
  geom_point(position=pd, size=5) +
  facet_wrap( ~ is_surprise_clip) +
  ylim(-14, 65) +
  xlab("Trials after the video clip") +
  ylab("Congruency Effect (ms)") +
  # ggtitle(title)
  geom_hline(yintercept=0, lty=2) +
  theme(legend.position="bottom")







# Effect of the experiment on the RTs -------------------------------------

# Are RTs longer, on average, when the video clips signal a 
# surprising event?

hist(final$rtTukey)
rtdata <- complete.cases(final$rtTukey)
retimes::mexgauss(rtdata)

# Parameter of the ex-gaussian distribution in the two conditions
# control condition
#         mu      sigma        tau 
# 0.93562532 0.04419088 0.05892118 

# surprise condition
#         mu      sigma        tau 
# 0.89897785 0.06645159 0.08860212 
# The difference is in the parameters sigma and tau.
# Better do this analysis with a Bayesian hierarchic model!



# Congruency effect for each subject and condition ------------------------

bysubraw <- final %>%
  group_by(
    subject_name, block, is_surprise_clip, trials_after_clip, 
    is_congruent_trial
  ) %>%
  summarise(
    mrt = mean(rtTukey, trim = 0.1, na.rm = TRUE)
  )

bysub_con <- bysubraw %>% dplyr::filter(
  is_congruent_trial == "Congruent"
)
bysub_inc <- bysubraw %>% dplyr::filter(
  is_congruent_trial == "Incongruent"
)

bysub_df <- merge(
  bysub_con,
  bysub_inc,
  by = c("subject_name", "block", "is_surprise_clip", 
         "trials_after_clip")
)

bysub_df <- bysub_df %>% 
  mutate(
    interference_eff = mrt.y - mrt.x
  )




plot_df <- bysub_df %>%
    group_by(is_surprise_clip, block, trials_after_clip) %>%
    summarise(
      y = mean(interference_eff, trim = .1, na.rm = TRUE),
      stderr = sqrt(var(interference_eff, na.rm = TRUE) / n()),
      n()
    )
plot_df$lower <- plot_df$y - plot_df$stderr
plot_df$upper <- plot_df$y + plot_df$stderr
  
pd <- position_dodge(0.2)
  
ggplot(plot_df, aes(x=trials_after_clip, y=y)) +
  geom_pointrange(aes(ymin=lower, ymax=upper), lwd=1.05, position=pd) +
  geom_line(position=pd, lwd=1.05) +
  geom_point(position=pd, size=5) +
  facet_wrap(block ~is_surprise_clip, ncol = 2) +
  xlab("Trials after the video clip") +
  ylab("Congruency Effect (ms)") +
  # ggtitle(title)
  geom_hline(yintercept=0, lty=2) +
  theme(legend.position="bottom")


plot_df <- bysub_df %>%
  group_by(is_surprise_clip, trials_after_clip) %>%
  summarise(
    y = mean(interference_eff, trim = .1, na.rm = TRUE),
    stderr = sqrt(var(interference_eff, na.rm = TRUE) / n()),
    n()
  )
plot_df$lower <- plot_df$y - plot_df$stderr
plot_df$upper <- plot_df$y + plot_df$stderr

pd <- position_dodge(0.2)

ggplot(plot_df, aes(x=trials_after_clip, y=y)) +
  geom_pointrange(aes(ymin=lower, ymax=upper), lwd=1.05, position=pd) +
  geom_line(position=pd, lwd=1.05) +
  geom_point(position=pd, size=5) +
  facet_wrap(~ is_surprise_clip, ncol = 2) +
  xlab("Trials after the video clip") +
  ylab("Congruency Effect (ms)") +
  # ggtitle(title)
  geom_hline(yintercept=0, lty=2) +
  theme(legend.position="bottom")


plot_df <- bysub_df %>%
  group_by(subject_name, is_surprise_clip) %>%
  summarise(
    y = mean(interference_eff, trim = .1, na.rm = TRUE)
  )
t.test(y ~ is_surprise_clip, data = plot_df)




plot_df <- bysub_df %>%
  
  group_by(block) %>%
  summarise(
    y = mean(interference_eff, trim = .1, na.rm = TRUE),
    stderr = sqrt(var(interference_eff, na.rm = TRUE) / n()),
    n()
  )
plot_df$lower <- plot_df$y - plot_df$stderr
plot_df$upper <- plot_df$y + plot_df$stderr


pd <- position_dodge(0.2)

ggplot(plot_df, aes(x=block, y=y)) +
  geom_pointrange(aes(ymin=lower, ymax=upper), lwd=1.05, position=pd) +
  geom_line(position=pd, lwd=1.05) +
  geom_point(position=pd, size=5) +
  xlab("Trials after the video clip") +
  ylab("Congruency Effect (ms)") +
  # ggtitle(title)
  geom_hline(yintercept=0, lty=2) +
  theme(legend.position="bottom")


t.test

BESTout <- BESTmcmc(y1, y2,  parallel=TRUE) 









bysub_df %>%
  group_by(block, is_clip_trial) %>%
  summarise(
    y = mean(interference_eff, na.rm = TRUE)
  )

length(unique(final$subject_name))


dd <- merge(bysub_df, temp, by = "subject_name")

bysub_df$acc <- temp$acc
bysub_df$med_rt <- temp$med_rt
bysub_df

mean(dd$med_rt) + c(-1, 1) * 2.5 * sd(dd$med_rt)
























# Stimulus sequence
with(
  df,
  data.frame(
    trial_within_blk, is_surprise_clip, is_clip_trial, trials_after_clip
  )[1:50, ]
)


temp <- mydata %>%
  dplyr::filter(
    block < 5
  )

ggerrorplot(temp,
            x = "block", 
            y = "correct",
            desc_stat = "mean_se",
            color = "is_congruent_trial", 
            palette = "jco",
            facet.by = "is_surprise_clip",
            ylab = "Accuracy",
            xlab = "Block",
            position = position_dodge(0.3) # Adjust the space between bars
)



ggerrorplot(temp,
            x = "block", 
            y = "rtTukey",
            desc_stat = "mean_se",
            color = "is_congruent_trial", 
            palette = "jco",
            facet.by = "is_surprise_clip",
            ylab = "Reaction Times (ms)",
            xlab = "Block",
            position = position_dodge(0.3) # Adjust the space between bars
)


ggerrorplot(temp,
            x = "is_clip_trial", 
            y = "rtTukey",
            desc_stat = "mean_se",
            color = "is_congruent_trial", 
            palette = "jco",
            facet.by = "is_surprise_clip",
            ylab = "Reaction Times (ms)",
            xlab = "Block",
            position = position_dodge(0.3) # Adjust the space between bars
)





# Histogram of accuracy
out <- mydata %>%
  group_by(subject_name, is_surprise_clip) %>%
  summarise(
    acc = mean(correct, na.rm = TRUE)
  )
# data.frame(out)

p <- ggdensity(out,
  x = "acc",
  add = "mean",
  rug = FALSE,
  color = "is_surprise_clip",
  fill = "is_surprise_clip",
  palette = c("#00AFBB", "#E7B800", "#FC4E07")
)

p2 <- ggpar(p,
  title = "Density Plot for Accuracy",
  subtitle = "Flanker task after video clip",
  caption = "Surprise project",
  xlab = "Accuracy",
  ylab = "Density",
  legend.title = "Surprise Video"
)
p2


# Histogram of median RTs as a function of surprise
out <- mydata %>%
  dplyr::filter(correct == 1) %>%
  group_by(subject_name, is_surprise_clip) %>%
  summarise(
    median_rt = median(rt)
  )
# data.frame(out)

p <- ggdensity(out,
  x = "median_rt",
  add = "mean", rug = FALSE,
  color = "is_surprise_clip", fill = "is_surprise_clip",
  palette = c("#00AFBB", "#E7B800", "#FC4E07")
)

p2 <- ggpar(p,
  title = "Density Plot for Response Latencies",
  subtitle = "Flanker task after video clip",
  caption = "Surprise project",
  xlab = "Median Reaction Times (ms)",
  ylab = "Density",
  legend.title = "Surprise Video"
)
p2


# Histogram of median RTs as a function of number of trials after the clip

temp <- mydata
temp$f_n_after <- factor(temp$trials_after_clip)
temp <- temp[temp$trials_after_clip != 4, ]

ggplot() + 
  geom_density(data=temp, aes(x=lrt, group=f_n_after, fill=f_n_after),alpha=0.5, adjust=2) + 
  xlab("Log reaction time") +
  ylab("Density")


# -------------------------------------------------------------------
# Log transformation of RTs works best
par(mfrow = c(1, 3))
qqPlot(mydata$trt, ylab = "Response Latencies (s)", main = "Reaction Times")
qqPlot(log(mydata$trt), ylab = "Log Response Latencies (s)", main = "Log Reaction Times")
qqPlot(1 / mydata$trt, ylab = "1 / RT", main = "Reciprocal of Reaction Times")
par(mfrow = c(1, 1))


ggerrorplot(final,
  x = "trials_after_clip", 
  y = "lrt",
  desc_stat = "mean_se",
  color = "is_congruent_trial", 
  palette = "jco",
  # facet.by = "block",
  ylab = "Log RTs",
  xlab = "Trial after the Video Clip",
  position = position_dodge(0.3) # Adjust the space between bars
)



ggerrorplot(four_trials_df,
  x = "trials_after_clip", 
  y = "lrt",
  desc_stat = "mean_se",
  color = "is_congruent_trial", 
  palette = "jco",
  facet.by = "is_surprise_clip",
  ylab = "Log RTs",
  xlab = "Trial after the Video Clip",
  position = position_dodge(0.3) # Adjust the space between bars
)


b1 <- final %>%
  dplyr::filter(block == 1)

ggerrorplot(b1,
  x = "trials_after_clip",
  y = "lrt",
  desc_stat = "mean_se",
  color = "is_congruent_trial", 
  palette = "jco",
  facet.by = "is_surprise_clip",
  ylab = "Log RTs",
  xlab = "Trials after the clip",
  position = position_dodge(0.3) # Adjust the space between bars
)


# -------------------------------------------------------------------
#  table of congruency effect and surprise
out <- final %>%
  dplyr::filter(correct == 1) %>%
  # group_by(block, is_surprise_clip, is_congruent_trial) %>%
  group_by(is_first_trial, is_congruent_trial) %>%
  summarise(
    m = mean(rt, na.rm = TRUE, trim = 0.1),
    se = sqrt(var(rt, na.rm = TRUE) / n()),
    n = n()
  )
data.frame(out)


# -------------------------------------------------------------------
# Plot reaction times
temp <- final %>%
  filter(trials_after_clip == 1)

ggerrorplot(final,
  x = "block",
  y = "lrt",
  desc_stat = "mean_se",
  color = "is_congruent_trial",
  palette = "jco",
  facet.by = "is_surprise_clip",
  ylab = "Log Reaction Times",
  xlab = "Block of Trials",
  position = position_dodge(0.3) # Adjust the space between bars
)


ggerrorplot(final,
            x = "block",
            y = "lrt",
            desc_stat = "mean_se",
            color = "is_surprise_clip",
            palette = "jco",
            facet.by = "is_congruent_trial",
            ylab = "Log Reaction Times",
            xlab = "Block of Trials",
            position = position_dodge(0.3) # Adjust the space between bars
)


bysub_df <- final %>%
  group_by(id, block, is_congruent_trial, is_surprise_clip) %>%
  summarise(
    mrt = mean(lrt, na.rm = TRUE), 
    n = n()
  )

bysub_df$mrt <- ifelse(bysub_df$n < 5, NA, bysub_df$mrt)

con_df <- dplyr::filter(bysub_df, is_congruent_trial == "Congruent")
con_df$is_congruent_trial <- NULL
inc_df <- dplyr::filter(bysub_df, is_congruent_trial == "Incongruent")
inc_df$is_congruent_trial <- NULL

con_df$n <- NULL
inc_df$n <- NULL

int_eff_df <- merge(con_df, inc_df, 
                    by = c("id", "block", "is_surprise_clip"))

int_eff_df$interf_eff <- int_eff_df$mrt.y - int_eff_df$mrt.x


ggerrorplot(int_eff_df,
            x = "block", 
            y = "interf_eff",
            desc_stat = "mean_se",
            color = "is_surprise_clip", 
            palette = "jco",
            ylab = "Log Reaction Times",
            xlab = "Block",
            position = position_dodge(0.3) # Adjust the space between bars
)






b1 <- final %>%
  filter(trials_after_clip == 1)

ggerrorplot(b1,
  x = "block", 
  y = "lrt",
  desc_stat = "mean_se",
  color = "is_congruent_trial", 
  palette = "jco",
  facet.by = c("is_surprise_clip"),
  ylab = "Log Reaction Times",
  xlab = "Block",
  position = position_dodge(0.3) # Adjust the space between bars
)


ggerrorplot(final,
  x = "is_surprise_clip", 
  y = "rt",
  desc_stat = "mean_se",
  palette = "jco",
  # facet.by = "is_surprise_clip",
  ylab = "Reaction Time (ms)",
  xlab = "Condition",
  position = position_dodge(0.3) # Adjust the space between bars
)


ggerrorplot(mydata,
  x = "is_surprise_clip", 
  y = "correct",
  desc_stat = "mean_se",
  palette = "jco",
  # facet.by = "is_surprise_clip",
  ylab = "Accuracy",
  xlab = "Surprise",
  position = position_dodge(0.3) # Adjust the space between bars
)


# -------------------------------------------------------------------
# Plot accuracy
all_trials_df <- df

all_trials_df <- all_trials_df %>%
  mutate(is_surprise_clip = fct_recode(is_surprise_clip,
    "No Surprise" = "No",
    "Surprise" = "Yes"
  ))

all_trials_df <- rename(all_trials_df, Condition = is_congruent_trial)
all_trials_df <- all_trials_df %>%
  mutate(Condition = fct_recode(Condition,
    "Incongruent" = "No",
    "Congruent" = "Yes"
  ))

ggerrorplot(all_trials_df,
  x = "block", y = "correct",
  desc_stat = "mean_se",
  color = "Condition", palette = "jco",
  facet.by = "is_surprise_clip",
  ylab = "Accuracy",
  xlab = "Block of Trials",
  ylim = c(0.75, 1.0),
  position = position_dodge(0.3) # Adjust the space between bars
)

ggerrorplot(all_trials_df,
  x = "is_surprise_clip", y = "correct",
  desc_stat = "mean_se",
  color = "condition", palette = "jco",
  facet.by = "block",
  ylab = "Accuracy",
  xlab = "Video Clip",
  ylim = c(0.75, 1.0),
  position = position_dodge(0.3) # Adjust the space between bars
)

all_trials_df <- all_trials_df[!is.na(all_trials_df$is_clip_prev_1) &
  all_trials_df$block < 5, ]

ggerrorplot(all_trials_df,
  x = "is_clip_prev_1", y = "correct",
  desc_stat = "mean_se",
  color = "condition", palette = "jco",
  facet.by = c("block"),
  ylab = "Accuracy",
  xlab = "First Trial After the Video Clip",
  ylim = c(0.75, 1.0),
  position = position_dodge(0.3) # Adjust the space between bars
)

# -------------------------------------------------------------------

names_df <- c(
  "id", "vel", "block", "i_trial", "video_id", "is_surprise_clip",
  "is_congruent_trial", "trials_after_clip"
)

dd <- dplyr::select(correct_df, names_df)

block_1_df <- dd %>%
  dplyr::filter(block == 4)

subject_means <- block_1_df %>%
  group_by(id, is_congruent_trial) %>%
  summarize(vel = mean(vel, na.rm = TRUE))
subject_means

# barplot <- ggplot(subject_means, aes(x = is_congruent_trial, y = vel)) +
#     stat_summary(
#     geom = "bar",
#     fun.y = "mean",
#     col = "black",
#     fill = "gray70"
#     )
# barplot

subject_means_wide <- spread(subject_means,
  key = is_congruent_trial,
  value = vel,
  sep = "_"
)
subject_means_wide

lims <- c(min(correct_df$vel, na.rm = TRUE), max(correct_df$vel, na.rm = TRUE))
wsplot <-
  ggplot(subject_means_wide, aes(
    x = is_congruent_trial_No,
    y = is_congruent_trial_Yes
  )) +
  geom_point() +
  geom_abline() +
  scale_x_continuous("Incongruent", limits = lims) +
  scale_y_continuous("Congruent", limits = lims) +
  theme(aspect.ratio = 1)
wsplot

# -------------------------------------------------------------------

message("\nDescriptive step: done!")
