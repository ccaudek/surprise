#' ---
#' title: "Descriptive statistics"
#' output:
#'   pdf_document:
#'     keep_tex: true
#' ---

#' Analysis of the congruency effect as a functio of experiment, subject, 
#' block, and is_surprise_clip.

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


data_tidy |> 
  group_by(experiment, is_surprise_clip) |> 
  summarize(
    n = n_distinct(subj_id)
  )

#' Select correct trials by experiment

dt_cor <- data_tidy |> 
  dplyr::filter(correct == 1)

# Select correct trials of the surprise experiment
surprise_cor_df <- dt_cor[dt_cor$experiment == "surprise", ]

# Select correct trials of the control experiment
control_cor_df <- dt_cor[dt_cor$experiment == "control", ]

#' Control: Compute the congruency effect for each subject

control_cor_df$blk <- factor(control_cor_df$block)

# Removing the first trial does not improve the results
# control_cor_df <- control_cor_df[control_cor_df$is_clip_trial == "No", ]

cntl_rt_diff_df <- control_cor_df %>%
  dplyr::filter(trials_after_clip != 0) |> 
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
  bf(ce ~ blk + (1 | subj_id:blk)),
  family = asym_laplace(),
  backend = "cmdstanr",
  init = 0.1,
  data = cntl_rt_diff_df
  #control = list(max_treedepth = 15, adapt_delta = 0.95)
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
  ylim(-35, 60) +
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
  dplyr::filter(trials_after_clip != 0) |> 
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


mod_1 <- brm(
  bf(ce ~ is_surprise_clip * blk +
    (1 + is_surprise_clip * blk | subj_id)),
  family = asym_laplace(),
  backend = "cmdstanr",
  init = 0.1,
  # algorithm = "meanfield",
  data = rt_diff_df
  #control = list(max_treedepth = 15, adapt_delta = 0.95)
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
  ylim(-35, 60) +
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
  family = asym_laplace(),
  backend = "cmdstanr",
  init = 0.1,
  data = dt_tot
)

loo_m_1 <- loo(m_1)
print(loo_m_1)

conditional_effects(m_1, "blk:exp")


# get the adjusted means
em <- emmeans(m_1, ~exp | blk)
em

# get all possible contrasts
cont <- contrast(em, "tukey")
cont



