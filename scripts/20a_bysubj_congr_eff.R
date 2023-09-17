#' ---
#' title: "Variational Inference analysis on raw data"
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
#
# Variational Inference with brms for the raw data

#' TODO: The results improve when I remove the subjects with only 1 or two blocks.
#' But there are also subjects who completed the task multiple times, more than 4 blocks.
#' I would be nice to keep only the first 4 blocks of them.

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

# temp <- flanker_data |>
#   dplyr::filter(bysub_n > 320)
# GwUYsO9QRc 8Pw2D108S6 TaAaeeUTFe CDvLepjySu pbXAI7g4iB ZNP649z7Mu

# one_subj <- flanker_data[flanker_data$subj_id == "ZNP649z7Mu", ]
# dim(one_subj)

sort(unique(flanker_data$bysub_n))

flanker_data_clean <- flanker_data |>
  dplyr::filter(bysub_n > 240 & bysub_n < 321)


# Check number of subjects by condition.
flanker_data_clean |>
  group_by(experiment, is_surprise_clip) |>
  summarize(
    n = n_distinct(subj_id)
  )

#' Select correct trials
dt_cor <- flanker_data_clean |>
  dplyr::filter(correct == 1)

nrow_total <- nrow(dt_cor)

# remove missing data on rt Tukey.
dt_cor <- dt_cor[!is.na(dt_cor$rtTukey), ]
nrow_na_removed <- nrow(dt_cor)

# percent removed
(1 - nrow_na_removed / nrow_total) * 100

nrow_total - nrow_na_removed
nrow_na_removed


#' Select correct trials by experiment

# Select correct trials of the surprise experiment
surprise_cor_df <- dt_cor[dt_cor$experiment == "surprise", ]

# Select correct trials of the control experiment
control_cor_df <- dt_cor[dt_cor$experiment == "control", ]


#' Control: Compute the congruency effect for each subject

control_cor_df$blk <- factor(control_cor_df$block)

# Removing the first trial does not improve the results
# control_cor_df <- control_cor_df[control_cor_df$is_clip_trial == "No", ]

cntl_rt_diff_df <- control_cor_df %>%
  group_by(subj_id, blk, is_congruent_trial) %>%
  summarise(mean_log_rt = mean(log(rtTukey), na.rm = TRUE)) %>%
  ungroup() %>%
  spread(key = is_congruent_trial, value = mean_log_rt) %>%
  mutate(
    mean_rt_congruent = exp(Congruent),
    mean_rt_incongruent = exp(Incongruent),
    ce = mean_rt_incongruent - mean_rt_congruent
  ) %>%
  arrange(subj_id, blk)


cntl_rt_diff_df |>
  group_by(blk) |>
  summarize(
    rt_diff = mean(ce),
    stderr = sqrt(var(ce) / n())
  )

plot(density(cntl_rt_diff_df$ce))

mod_cntl_1 <- brm(
  bf(ce ~ blk + (1 + blk | subj_id)),
  family = student(),
  backend = "cmdstanr",
  init = 0.1,
  data = cntl_rt_diff_df
)

pp_check(mod_cntl_1)

loo_mod_cntl_1 <- loo(mod_cntl_1)
print(loo_mod_cntl_1)

summary(mod_cntl_1)

conditional_effects(mod_cntl_1, "blk")

c_eff <- conditional_effects(mod_cntl_1, "blk")
p_c <- plot(c_eff, plot = FALSE)[[1]] +
  theme(legend.position = "bottom") +
  labs(
    y = "Congruency Effect (ms)",
    x = "Block of Trials"
  ) +
  ylim(-30, 60) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "darkgray")
print(p_c)




# get the adjusted means
em <- emmeans(mod_cntl_1, ~blk)
em

# get all possible contrasts
cont <- contrast(em, "tukey")
cont

hyp1 <- c(hyp1 = "Intercept > 0")
hypothesis(mod_cntl_1, hyp1)

hyp2 <- c(hyp1 = "Intercept + blk2 > 0")
hypothesis(mod_cntl_1, hyp2)



#' Surprise: Compute the congruency effect for each subject

surprise_cor_df$blk <- factor(surprise_cor_df$block)

rt_diff_df <- surprise_cor_df |>
  group_by(subj_id, blk, is_surprise_clip, is_congruent_trial) |>
  summarise(mean_log_rt = mean(log(rtTukey), na.rm = TRUE)) |>
  ungroup() |>
  spread(key = is_congruent_trial, value = mean_log_rt) |>
  mutate(
    mean_rt_congruent = exp(Congruent),
    mean_rt_incongruent = exp(Incongruent),
    ce = mean_rt_incongruent - mean_rt_congruent
  ) %>%
  arrange(subj_id, blk, is_surprise_clip) |>
  ungroup()

rt_diff_df |>
  group_by(is_surprise_clip, blk) |>
  summarize(
    rt_diff = mean(ce),
    stderr = sqrt(var(ce) / n())
  )

plot(density(rt_diff_df$ce))

rt_diff_df$blk <- factor(rt_diff_df$block)

rt_diff_df$zce <- rt_diff_df$ce |>
  scale() |>
  as.numeric()

mod_1 <- brm(
  bf(ce ~ is_surprise_clip * blk +
    (1 + is_surprise_clip * blk | subj_id)),
  family = student(),
  backend = "cmdstanr",
  init = 0.1,
  data = rt_diff_df
)

pp_check(mod_1)

loo_mod_1 <- loo(mod_1)
print(loo_mod_1)
# Computed from 4000 by 958 log-likelihood matrix
#
# Estimate   SE
# elpd_loo  -1007.9 28.1
# p_loo       219.0  6.2
# looic      2015.7 56.2
# ------
#   Monte Carlo SE of elpd_loo is 0.4.
#
# Pareto k diagnostic values:
#   Count Pct.    Min. n_eff
# (-Inf, 0.5]   (good)     951   99.3%   132
#  (0.5, 0.7]   (ok)         7    0.7%   697
#    (0.7, 1]   (bad)        0    0.0%   <NA>
#    (1, Inf)   (very bad)   0    0.0%   <NA>
#
#   All Pareto k estimates are ok (k < 0.7).
# See help('pareto-k-diagnostic') for details.

c_eff <- conditional_effects(mod_1, "blk:is_surprise_clip")
p_s <- plot(c_eff, plot = FALSE)[[1]] +
  theme(legend.position = "bottom") +
  labs(
    y = "Congruency Effect (ms)",
    x = "Block of Trials"
  ) +
  ylim(-30, 60) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "darkgray")
print(p_s)

p <- p_c + p_s
print(p)
ggsave("congr_effect.pdf", width = 10, height = 6)






#' both data set

rt_diff_df$exp <- "surprise"
cntl_rt_diff_df$exp <- "control"
cntl_rt_diff_df$is_surprise_clip <- "No Surprise"

dt_tot <- bind_rows(rt_diff_df, cntl_rt_diff_df) |>
  dplyr::select(exp, subj_id, is_surprise_clip, blk, ce)

dt_tot$exp <- factor(dt_tot$exp)
dt_tot$is_surprise_clip <- factor(dt_tot$is_surprise_clip)

summary(dt_tot)

m_1 <- brm(
  bf(ce ~ exp * blk + (1 * blk | subj_id)),
  family = student(),
  backend = "cmdstanr",
  init = 0.01,
  data = dt_tot
)

loo_m_1 <- loo(m_1)
print(loo_m_1)

m_2 <- brm(
  bf(ce ~ exp + blk + (1 + blk | subj_id)),
  family = student(),
  backend = "cmdstanr",
  init = 0.01,
  data = dt_tot
)

loo_m_2 <- loo(m_2)
print(loo_m_2)

comp <- loo_compare(loo_m_1, loo_m_2)
comp

summary(m_2)

m_3 <- brm(
  bf(ce ~ exp + (1 | subj_id)),
  family = student(),
  backend = "cmdstanr",
  init = 0.01,
  data = dt_tot
)
loo_m_3 <- loo(m_3)

comp <- loo_compare(loo_m_2, loo_m_3)
comp

m_4 <- brm(
  bf(ce ~ blk + (1 + blk | subj_id)),
  family = student(),
  backend = "cmdstanr",
  init = 0.01,
  data = dt_tot
)
loo_m_4 <- loo(m_4)

comp <- loo_compare(loo_m_1, loo_m_4)
comp



# make_stancode(
#   ce ~ exp * is_surprise_clip * blk +
#     (1 + is_surprise_clip * blk | subj_id),
#   family = student()
# )

pp_check(m_1)

loo_m_1 <- loo(m_1)
print(loo_m_1)


conditions <- make_conditions(m_1, "exp")
conditional_effects(m_1, "blk:exp", conditions = conditions)

conditional_effects(m_1, "blk:exp")
conditional_effects(m_1, "exp")

summary(m_1)

# Overlay data points from rt_diff_df


s_p <- for_plot |>
  dplyr::filter(is_surprise_clip == "Surprise")
s_p$blk <- as.numeric(s_p$blk)

p + geom_point(aes(x = s_p$blk, y = s_p$rt_diff), color = "black", shape = 3, size = 3)


comp <- loo_compare(loo_m0, loo_m1)
print(comp, digits = 2)



#' No Surprise: Compute the congruency effect for each subject

cntl_rt_diff_df <- control_cor_df %>%
  group_by(subj_id, blk, is_congruent_trial) %>%
  summarise(mean_log_rt = mean(log(rtTukey), na.rm = TRUE)) %>%
  ungroup() %>%
  spread(key = is_congruent_trial, value = mean_log_rt) %>%
  mutate(
    mean_rt_congruent = exp(Congruent),
    mean_rt_incongruent = exp(Incongruent),
    ce = mean_rt_incongruent - mean_rt_congruent
  ) %>%
  arrange(subj_id, blk)

# cntl_rt_diff_df <- control_cor_df |>
#   group_by(subj_id, blk, is_congruent_trial) |>
#   summarise(mean_rt = median(m_rt, na.rm = TRUE)) |>
#   spread(key = is_congruent_trial, value = mean_rt) |>
#   mutate(ce = Incongruent - Congruent) |>
#   arrange(subj_id, blk) |>
#   ungroup()

cntl_rt_diff_df$zce <- scale(cntl_rt_diff_df$ce) |> as.numeric()

cntl_rt_diff_df |>
  group_by(blk) |>
  summarize(
    rt_diff = mean(ce),
    stderr = sqrt(var(ce) / n()),
    n = n_distinct(subj_id),
    zrt_ce = mean(zce)
  )

plot(density(cntl_rt_diff_df$zce))

cntl_rt_diff_df$blk <- factor(cntl_rt_diff_df$block)

prior <- prior_string("normal(0, 1)", class = "b")

cmod_1 <- brm(
  bf(ce ~ blk +
    (1 + blk | subj_id)),
  family = student(),
  prior = prior,
  backend = "cmdstanr",
  init = 0.1,
  data = cntl_rt_diff_df
)

pp_check(cmod_1)


loo_cmod_1 <- loo(cmod_1)
print(loo_cmod_1)
# Computed from 4000 by 958 log-likelihood matrix
#
# Estimate   SE
# elpd_loo  -5440.8 29.6
# p_loo       189.7  3.6
# looic     10881.6 59.2
# ------
#   Monte Carlo SE of elpd_loo is 0.3.
#
# All Pareto k estimates are good (k < 0.5).
# See help('pareto-k-diagnostic') for details.

c_eff <- conditional_effects(cmod_1, "blk")

p <- plot(c_eff, plot = FALSE)[[1]] +
  theme(legend.position = "bottom") +
  labs(
    y = "Congruency Effect (ms)",
    x = "Block of Trials"
  )
print(p)






















# ------------------------------------------------------------------------------
# Compute the congruency effect for each subject for both experiments
# ------------------------------------------------------------------------------

dt_cor$blk <- factor(dt_cor$block)

rt_diff_df <- dt_cor |>
  group_by(experiment, subj_id, blk, is_surprise_clip, is_congruent_trial) |>
  summarise(
    m_rt = median(rt, na.rm = TRUE),
    n = n()
  ) %>%
  group_by(experiment, subj_id, blk, is_surprise_clip, is_congruent_trial) |>
  summarise(mean_rt = mean(m_rt, na.rm = TRUE)) |>
  spread(key = is_congruent_trial, value = mean_rt) |>
  mutate(ce = Incongruent - Congruent) |>
  arrange(experiment, subj_id, blk, is_surprise_clip) |>
  ungroup()

rt_diff_df |>
  group_by(experiment, blk, is_surprise_clip) |>
  summarize(
    rt_diff = mean(ce),
    stderr = sqrt(var(ce) / n())
  )

plot(density(rt_diff_df$ce))

rt_diff_df$blk <- factor(rt_diff_df$block)

mod1 <- brm(
  bf(ce ~ experiment * is_surprise_clip * blk +
    (1 + is_surprise_clip * blk | subj_id)),
  family = student(),
  backend = "cmdstanr",
  init = 0.1,
  data = rt_diff_df
)

pp_check(mod1)

loo_mod1 <- loo(mod1)
print(loo_mod1)


# Three-way interaction
conditions <- make_conditions(mod1, "experiment")
conditional_effects(mod1, "blk:is_surprise_clip", conditions = conditions)















# Step 1: Log-transform the rt values
surprise_cor_df <- surprise_cor_df |>
  mutate(log_rt = log(rt))

# Step 2 & 3: Calculate the mean log_rt for each condition and each subj_id,
# and then find the difference between the conditions
rt_diff <- surprise_cor_df |>
  group_by(subj_id, block, is_surprise_clip, is_congruent_trial) |>
  summarise(
    mean_log_rt = mean(log_rt, na.rm = TRUE),
    n = n()
  ) %>%
  group_by(subj_id, block, is_surprise_clip, is_congruent_trial) |>
  summarise(mean_log_rt = mean(mean_log_rt, na.rm = TRUE)) |>
  spread(key = is_congruent_trial, value = mean_log_rt) |>
  mutate(rt_diff_log_scale = Congruent - Incongruent) |>
  arrange(subj_id, block, is_surprise_clip) |>
  ungroup()

# Step 4: Exponentiate the differences to get back to the original scale
rt_diff <- rt_diff |>
  mutate(congr_eff = exp(rt_diff_log_scale))


rt_diff |>
  group_by(block, is_surprise_clip) |>
  summarize(
    m = mean(congr_eff)
  )




# ------------------------------------------------------------------------------
# Data wrangling
# ------------------------------------------------------------------------------

surprise_cor_df$BL <- surprise_cor_df$block
surprise_cor_df$blk <- factor(surprise_cor_df$block)
surprise_cor_df$BF <- surprise_cor_df$blk

surprise_cor_df$zrt <- scale(surprise_cor_df$rt) |> as.numeric()

surprise_cor_df$CT <- surprise_cor_df$is_congruent_trial |>
  as.factor()

surprise_cor_df$SC <- surprise_cor_df$is_surprise_clip |>
  as.factor()

surprise_cor_df$movie_id <- factor(surprise_cor_df$movie_id)

# I tried to remove first trial after the video-clip, but the results
# are worse.

d <- surprise_cor_df |>
  dplyr::select(rt, zrt, CT, SC, BL, BF, subj_id, movie_id)

# rio::export(d, "surprise_correct_data.csv")


# ------------------------------------------------------------------------------
# brm() analysis
# ------------------------------------------------------------------------------


# ---- Baseline model ----

m0 <- brm(
  bf(
    zrt ~ 1 + (1 | subj_id) + (1 | movie_id)
  ),
  algorithm = "meanfield",
  family = asym_laplace(),
  iter = 20000, # Increase the number of iterations
  data = d
)
loo_m0 <- loo(m0)


# ---- Model without random slopes ----

m1 <- brm(
  bf(zrt ~ CT * SC * BL +
    (1 | subj_id) + (1 | movie_id)),
  algorithm = "meanfield",
  family = asym_laplace(),
  iter = 20000, # Increase the number of iterations
  data = d
)
loo_m1 <- loo(m1)

comp <- loo_compare(loo_m0, loo_m1)
print(comp, digits = 2)
#    elpd_diff se_diff
# m1    0.00      0.00 <- preferred model
# m0 -897.70     51.26


# ---- Model with block as factor ----

m2 <- brm(
  bf(zrt ~ CT * SC * BF +
    (1 + CT * SC * BF | subj_id) + (1 | movie_id)),
  algorithm = "meanfield",
  family = asym_laplace(),
  iter = 20000, # Increase the number of iterations
  data = d
)
loo_m2 <- loo(m2)

comp <- loo_compare(loo_m1, loo_m2)
print(comp, digits = 2)
#    elpd_diff se_diff
# m2    0.00      0.00 <- preferred model
# m1 -414.50     62.35

print(loo_m2)
# Computed from 1000 by 36032 log-likelihood matrix
#
# Estimate    SE
# elpd_loo -33961.8 186.7
# p_loo      1735.5  11.4
# looic     67923.7 373.4
# ------
#   Monte Carlo SE of elpd_loo is 1.5.
#
# Pareto k diagnostic values:
#   Count Pct.    Min. n_eff
# (-Inf, 0.5]   (good)     36002 99.9%   481
#  (0.5, 0.7]   (ok)          30  0.1%   750
#    (0.7, 1]   (bad)          0  0.0%   <NA>
#    (1, Inf)   (very bad)     0  0.0%   <NA>
#
#   All Pareto k estimates are ok (k < 0.7).

pp_check(m2)

# ---- Model with block as numeric ----

m3 <- brm(
  bf(zrt ~ CT * SC * BL +
    (1 + CT * SC * BL | subj_id) + (1 | movie_id)),
  algorithm = "meanfield",
  family = asym_laplace(),
  iter = 20000, # Increase the number of iterations
  init = 0.01,
  data = d
)
loo_m3 <- loo(m3)

comp <- loo_compare(loo_m3, loo_m2)
print(comp, digits = 2)
#    elpd_diff se_diff
# m3    0.00      0.00
# m2 -479.28     44.14

pp_check(m3)
# fit is ok.

print(loo_m3)
# Computed from 1000 by 36032 log-likelihood matrix
#
# Estimate    SE
# elpd_loo -33482.6 188.7
# p_loo      1321.7   9.6
# looic     66965.1 377.3
# ------
#   Monte Carlo SE of elpd_loo is 1.3.
#
# Pareto k diagnostic values:
#   Count Pct.    Min. n_eff
# (-Inf, 0.5]   (good)     36003 99.9%   363
#  (0.5, 0.7]   (ok)          29  0.1%   805
#    (0.7, 1]   (bad)          0  0.0%   <NA>
#    (1, Inf)   (very bad)     0  0.0%   <NA>
#
#   All Pareto k estimates are ok (k < 0.7).
# See help('pareto-k-diagnostic') for details.

# ------------------------------------------------------------------------------

# Test whether the kind of video is important within the surprise experiment.

m4 <- brm(
  bf(zrt ~ CT * BL +
    (1 + CT * BL | subj_id) + (1 | movie_id)),
  algorithm = "meanfield",
  family = asym_laplace(),
  iter = 20000, # Increase the number of iterations
  init = 0.01,
  data = d
)
loo_m4 <- loo(m4)


comp <- loo_compare(loo_m3, loo_m4)
print(comp, digits = 2)
#     elpd_diff se_diff
# m3     0.00      0.00
# m4 -1752.48     81.93

# Interpretation:
# The variable coding the kind of video within the surprise experiment
# is necessary in the model.

# ------------------------------------------------------------------------------

# Remove the three-way interaction

m5 <- brm(
  bf(zrt ~ CT * SC + CT * BL + SC * BL +
    (1 + CT * SC + CT * BL + SC * BL | subj_id) + (1 | movie_id)),
  algorithm = "meanfield",
  family = brms::asym_laplace(),
  iter = 20000, # Increase the number of iterations
  init = 0.01,
  data = d
)
loo_m5 <- loo(m5)

# Test of the three-way interaction
comp <- loo_compare(loo_m3, loo_m5)
print(comp, digits = 2)
#     elpd_diff se_diff
# m3     0.00      0.00
# m5 -2696.98     80.65


# ------------------------------------------------------------------------------

# Conditional plot.

mod <- m1

# Three-way interaction
conditions <- make_conditions(mod, "SC")
c_eff <- conditional_effects(mod, "BL:CT", conditions = conditions)
plot(c_eff, plot = FALSE)[[1]] +
  theme(legend.position = "bottom") +
  labs(
    y = "Reaction Times (standardized)",
    x = "Block of Trials"
  )


# ------------------------------------------------------------------------------

# Plot of the raw data means

plot(density(log(surprise_cor_df$rt)))

surprise_cor_df$lrt <- log(surprise_cor_df$rt)

# Calculate the within-subject mean and standard error for each condition
subject_summary <- surprise_cor_df %>%
  group_by(subj_id, SC, CT, BL) %>%
  summarize(
    subj_mean = mean(lrt, na.rm = TRUE),
    .groups = "drop"
  )

# Calculate the overall mean and within-subject standard error for each condition
plot_df <- subject_summary %>%
  group_by(SC, CT, BL) %>%
  summarize(
    m = mean(subj_mean, na.rm = TRUE),
    stderr = sd(subj_mean, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Create the lower and upper bounds for the error bars
plot_df$lower <- plot_df$m - plot_df$stderr
plot_df$upper <- plot_df$m + plot_df$stderr

plot_df$m <- exp(plot_df$m)
plot_df$stderr <- exp(plot_df$stderr)
plot_df$lower <- exp(plot_df$lower)
plot_df$upper <- exp(plot_df$upper)

# Create the plot
pd <- position_dodge(0.5)

ggplot(plot_df, aes(x = BL, y = m, color = CT)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), lwd = 1.05, position = pd) +
  geom_line(position = pd, lwd = 1.05) +
  geom_point(position = pd, size = 5) +
  facet_grid(~SC) +
  xlab("Block") +
  ylab("RT (ms)") +
  theme(legend.position = "bottom")

# ------------------------------------------------------------------------------

# Compute contrasts.

mod <- m2 # Use block as factor

# get the adjusted means
em <- emmeans(mod, ~ BF:CT | SC)
em

# get all possible contrasts
cont <- contrast(em, "tukey")
cont
# SC = No Surprise:
#   contrast                          estimate lower.HPD upper.HPD
# BF1 Congruent - BF2 Congruent      0.32631   0.29430   0.35578
# BF1 Congruent - BF3 Congruent      0.44675   0.41391   0.47644
# BF1 Congruent - BF4 Congruent      0.49250   0.45221   0.53113
# BF1 Congruent - BF1 Incongruent    0.19668   0.17075   0.22125 <--
# BF1 Congruent - BF2 Incongruent    0.24445   0.19228   0.30565
# BF1 Congruent - BF3 Incongruent    0.36000   0.31146   0.42427
# BF1 Congruent - BF4 Incongruent    0.45139   0.38309   0.51873
# BF2 Congruent - BF3 Congruent      0.12009   0.07644   0.16764
# BF2 Congruent - BF4 Congruent      0.16695   0.11851   0.22396
# BF2 Congruent - BF1 Incongruent   -0.12938  -0.17373  -0.08870
# BF2 Congruent - BF2 Incongruent   -0.08079  -0.12907  -0.03532
# BF2 Congruent - BF3 Incongruent    0.03611  -0.02806   0.10334
# BF2 Congruent - BF4 Incongruent    0.12492   0.04511   0.19469
# BF3 Congruent - BF4 Congruent      0.04716  -0.00305   0.09602
# BF3 Congruent - BF1 Incongruent   -0.24895  -0.29246  -0.20968
# BF3 Congruent - BF2 Incongruent   -0.20214  -0.27587  -0.13989
# BF3 Congruent - BF3 Incongruent   -0.08467  -0.12743  -0.03278
# BF3 Congruent - BF4 Incongruent    0.00452  -0.06401   0.08327
# BF4 Congruent - BF1 Incongruent   -0.29522  -0.34456  -0.24918
# BF4 Congruent - BF2 Incongruent   -0.24843  -0.31901  -0.17724
# BF4 Congruent - BF3 Incongruent   -0.13031  -0.19320  -0.05411
# BF4 Congruent - BF4 Incongruent   -0.04127  -0.08958   0.01604
# BF1 Incongruent - BF2 Incongruent  0.04865  -0.00251   0.09942
# BF1 Incongruent - BF3 Incongruent  0.16402   0.11625   0.21584
# BF1 Incongruent - BF4 Incongruent  0.25500   0.19797   0.32092
# BF2 Incongruent - BF3 Incongruent  0.11774   0.04505   0.18150
# BF2 Incongruent - BF4 Incongruent  0.20672   0.12766   0.29031
# BF3 Incongruent - BF4 Incongruent  0.08991   0.00297   0.16201
#
# SC = Surprise:
#   contrast                          estimate lower.HPD upper.HPD
# BF1 Congruent - BF2 Congruent      0.19315   0.14267   0.24610
# BF1 Congruent - BF3 Congruent      0.35223   0.30370   0.40369
# BF1 Congruent - BF4 Congruent      0.36052   0.29377   0.42024
# BF1 Congruent - BF1 Incongruent   -0.01154  -0.05246   0.03207
# BF1 Congruent - BF2 Incongruent    0.25841   0.16531   0.34981
# BF1 Congruent - BF3 Incongruent    0.41177   0.32316   0.51147
# BF1 Congruent - BF4 Incongruent    0.48783   0.38370   0.58670
# BF2 Congruent - BF3 Congruent      0.15877   0.08179   0.22615
# BF2 Congruent - BF4 Congruent      0.16804   0.08894   0.24930
# BF2 Congruent - BF1 Incongruent   -0.20424  -0.26618  -0.13863
# BF2 Congruent - BF2 Incongruent    0.06509  -0.01224   0.13943
# BF2 Congruent - BF3 Incongruent    0.21839   0.12111   0.33398
# BF2 Congruent - BF4 Incongruent    0.29485   0.17489   0.40322
# BF3 Congruent - BF4 Congruent      0.00858  -0.07214   0.08905
# BF3 Congruent - BF1 Incongruent   -0.36370  -0.43144  -0.30122
# BF3 Congruent - BF2 Incongruent   -0.09618  -0.20223   0.01206
# BF3 Congruent - BF3 Incongruent    0.05873  -0.02254   0.13822
# BF3 Congruent - BF4 Incongruent    0.13316   0.02524   0.25172
# BF4 Congruent - BF1 Incongruent   -0.37217  -0.44555  -0.29326
# BF4 Congruent - BF2 Incongruent   -0.10417  -0.21060   0.00588
# BF4 Congruent - BF3 Incongruent    0.05144  -0.05794   0.16684
# BF4 Congruent - BF4 Incongruent    0.12594   0.04637   0.21173
# BF1 Incongruent - BF2 Incongruent  0.26826   0.18546   0.35157
# BF1 Incongruent - BF3 Incongruent  0.42335   0.33276   0.50296
# BF1 Incongruent - BF4 Incongruent  0.49816   0.41057   0.60008
# BF2 Incongruent - BF3 Incongruent  0.15391   0.02994   0.26440
# BF2 Incongruent - BF4 Incongruent  0.22878   0.10895   0.35214
# BF3 Incongruent - BF4 Incongruent  0.07672  -0.04414   0.20851
#
# Point estimate displayed: median
# HPD interval probability: 0.95

# get the posterior draws from the contrasts
cont_posterior <- gather_emmeans_draws(cont)

# plot
ggplot(
  cont_posterior,
  aes(y = contrast, x = .value)
) +
  stat_halfeye() +
  facet_wrap(~SC) +
  geom_vline(xintercept = 0, color = "red", lty = 2)


# ------------------------------------------------------------------------------
# Control experiment
# ------------------------------------------------------------------------------

control_cor_df$blk <- factor(control_cor_df$block)

control_cor_df$zrt <- scale(control_cor_df$rt) |>
  as.numeric()

control_cor_df$CT <- control_cor_df$is_congruent_trial |>
  as.factor()

control_cor_df$BF <- control_cor_df$blk
control_cor_df$BL <- control_cor_df$block

control_cor_df$movie_id <- factor(control_cor_df$movie_id)
control_cor_df$subj_id <- factor(control_cor_df$subj_id)


dc <- control_cor_df |>
  dplyr::select(zrt, CT, BL, BF, subj_id, movie_id)


# ---- Model with block as numeric ----

c3 <- brm(
  bf(zrt ~ CT * BL +
    (1 + CT * BL | subj_id) + (1 | movie_id)),
  algorithm = "meanfield",
  family = brms::asym_laplace(),
  iter = 20000, # Increase the number of iterations
  init = 0.01,
  data = dc
)
loo_c3 <- loo(c3)

print(loo_c3)
# Computed from 1000 by 23818 log-likelihood matrix
#
# Estimate    SE
# elpd_loo -25123.0 154.2
# p_loo       728.5   6.1
# looic     50246.0 308.4
# ------
#   Monte Carlo SE of elpd_loo is 0.9.
#
# Pareto k diagnostic values:
#   Count Pct.    Min. n_eff
# (-Inf, 0.5]   (good)     23814 100.0%  494
# (0.5, 0.7]   (ok)           4   0.0%  907
# (0.7, 1]   (bad)          0   0.0%  <NA>
#   (1, Inf)   (very bad)     0   0.0%  <NA>
#
#   All Pareto k estimates are ok (k < 0.7).
# See help('pareto-k-diagnostic') for details.

summary(c3)
# Population-Level Effects:
#                  Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept           -0.30      0.04    -0.38    -0.23 1.01      847      790
# CTIncongruent        0.17      0.03     0.11     0.23 1.00      729      831
# BL                  -0.10      0.01    -0.12    -0.08 1.00     1082     1070
# CTIncongruent:BL     0.04      0.01     0.02     0.06 1.00      979      952

mod_c <- c3

conditional_effects(mod_c, "BL:CT")

# Two-way interaction
c_eff <- conditional_effects(mod_c, "BL:CT")
plot(c_eff, plot = FALSE)[[1]] +
  theme(legend.position = "bottom") +
  labs(
    y = "Reaction Times (standardized)",
    x = "Block of Trials"
  )


# ------------------------------------------------------------------------------

# No interaction
c4 <- brm(
  bf(zrt ~ CT + BL +
    (1 + CT + BL | subj_id) + (1 | movie_id)),
  algorithm = "meanfield",
  family = brms::asym_laplace(),
  iter = 20000, # Increase the number of iterations
  init = 0.1,
  data = dc
)
loo_c4 <- loo(c4)

# Test of the interaction
comp <- loo_compare(loo_c3, loo_c4)
print(comp, digits = 2)
#     elpd_diff se_diff
# c3     0.00      0.00
# c4 -1166.27     63.11


message("\n20_variational_inference.R: done!")

# eof ----





flanker_data_clean |>
  group_by(experiment) |>
  summarize(
    n = n_distinct(subj_id)
  )

#' models fit

modeling_data_cntl <- flanker_data_clean |>
  dplyr::filter(experiment == "control")

unique(modeling_data_cntl$subj_id)

subset_ss <- c("EotVPAq8iy", "xmbeMys6sI")

temp <- modeling_data_cntl[modeling_data_cntl$subj_id %in% subset_ss, ]
temp$subj_id <- factor(temp$subj_id)

temp1 <- temp[temp$block == 1, ]

temp1$rt <- temp1$rtTukey / 1000

temp1$subject <- as.integer(temp1$subj_id)
sort(unique(temp1$subject))

temp1$congruency <- ifelse(temp1$is_congruent_trial == "Congruent", "congruent", "incongruent")
temp1$accuracy <- temp1$correct

temp2 <- temp1 |>
  dplyr::select(subject, congruency, accuracy, rt)

complete_rows <- complete.cases(temp2)

# Subset the data frame to keep only complete rows
dd <- temp2[complete_rows, ]

dd$subj_name <- NULL

# ======================================
# set some parameters for the model
# fit routines
# ======================================

# during the fit, how many sets of starting parameters should be explored?
n_start_parms <- 50

# what should the variance across starting parameters be?
var_start_parms <- 20

# how many trials to simulate during each iteration of the fit routine whilst
# exploring multiple starting parameters?
n_first_pass <- 1000

# how many trials to simulate during the final fit routine?
n_final_pass <- 50000

# pass the data to the model_fit functions. This is the funcion that fits the
# SSP and DSTP model. Note there are two functions being called:
# 1. fit_[model]_group: fit each model to group-level data. This is used to
#    assess on a group level which model fits the data better. [TODO]
# 2. fit_[model]_individual: fit each model to individual participants' data.
#    This is used to get parameters per participant, used for the correlations
#    to depression & pleasure scales.

library(flankr)

# fit each model to individual subjects



res_control <- list()

n_sub <- 2


for (i_sub in 1:n_sub) {
  # select the data of one subject
  one_subj <- dd %>%
    dplyr::filter(
      subject == i_sub
    )
  # fit the DSTP model
  m <- fitDSTP(data = one_subj)
  # save results
  res_control[[i_sub]] <- m$bestParameters
  print(i_sub)
}

param_names <- c("A", "C", "mu_ta", "mu_fl", "mu_ss", "mu_rs2", "ter")

cntl_exp_df <- do.call(rbind.data.frame, res_control)
colnames(cntl_exp_df) <- param_names


#' SSP
res_ssp_control <- list()

n_sub <- 2

for (i_sub in 1:n_sub) {
  # select the data of one subject
  one_subj <- dd %>%
    dplyr::filter(
      subject == i_sub
    )
  # fit the SSP model
  m <- fitSSP(data = one_subj)
  # save results
  res_ssp_control[[i_sub]] <- m$bestParameters
  print(i_sub)
}

param_names <- c("A", "ter", "p", "rd", "sda")

res_ssp_control <- do.call(rbind.data.frame, res_ssp_control)
colnames(res_ssp_control) <- param_names
