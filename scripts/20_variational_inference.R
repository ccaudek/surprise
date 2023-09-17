#' ---
#' title: "Variational Inference analysis on raw data"
#' output:
#'   pdf_document:
#'     keep_tex: true
#' ---

library("here")

suppressPackageStartupMessages(
  {
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
  }
)

theme_set(bayesplot::theme_default(base_family = "sans", base_size = 14))
set.seed(123)

source(here("libraries", "functions.R"))


#' Import complete data set


#' Import the data that have been created by the previous scripts. The data 
#' have been created with the 01_, 10_, 11_ scripts in the present directory.

data <- get_data()


#' Tidy data

data_tidy <- tidy_flanker(data)

#' Perform some participants' flanker checks.

flanker_accuracy_overall <- get_flanker_accuracy(data_tidy, overall = TRUE)

#' Get a list of participants who scored below 80% accuracy.

accuracy_removal <- flanker_accuracy_overall |> 
  filter(accuracy < 0.80) |> 
  pull(subj_id)

length(accuracy_removal)

#' Remove the <80% accuracy participants from the flanker data.

flanker_data <- data_tidy |> 
  filter(!subj_id %in% accuracy_removal)

#' Check number of subjects by condition.

flanker_data |>
  group_by(experiment, is_surprise_clip) |>
  summarize(
    n = n_distinct(subj_id)
  ) 


#' Select correct trials only

dt_cor <- flanker_data |> 
  dplyr::filter(correct == 1)

nrow_total <- nrow(dt_cor)

#' Remove missing data on rt.

dt_cor <- dt_cor[!is.na(dt_cor$rt), ]

nrow_na_removed <- nrow(dt_cor)

# percent removed
(1 - nrow_na_removed / nrow_total) * 100

nrow_total - nrow_na_removed  

nrow_total


#' Select correct trials by experiment

#' Select correct trials of the surprise experiment 
surprise_cor_df <- dt_cor[dt_cor$experiment == "surprise", ]

#' Select correct trials of the control experiment 
control_cor_df <- dt_cor[dt_cor$experiment == "control", ]


#' Data wrangling

surprise_cor_df$BL <- surprise_cor_df$block
surprise_cor_df$blk <- factor(surprise_cor_df$block)
surprise_cor_df$BF <- surprise_cor_df$blk

surprise_cor_df$zrt <- scale(surprise_cor_df$rt) |> as.numeric()

surprise_cor_df$CT <- surprise_cor_df$is_congruent_trial |>
  as.factor()

surprise_cor_df$SC <- surprise_cor_df$is_surprise_clip |>
  as.factor()

surprise_cor_df$movie_id <- factor(surprise_cor_df$movie_id)

d <- surprise_cor_df |>
  dplyr::select(rt, zrt, CT, SC, BL, BF, subj_id, movie_id) 


#' Behavioural analysis of the flanker effect
#' Response time

#' Baseline model
#' ex-gaussian produces a better fit that asym-laplace

m0 <- brm(
  bf(
    zrt ~ 1 + (1 | subj_id) + (1 | movie_id)
  ), 
  algorithm = "meanfield",
  family = exgaussian(),
  iter = 20000, # Increase the number of iterations
  data = d
)
loo_m0 <- loo(m0)


#' Model without random slopes 

m1 <- brm(
  bf(zrt ~ CT * SC * BL +
       (1 | subj_id) + (1 | movie_id)
  ), 
  algorithm = "meanfield",
  family = exgaussian(),
  iter = 20000, 
  data = d
)
loo_m1 <- loo(m1)

comp <- loo_compare(loo_m0, loo_m1)
print(comp, digits = 2)


#' Model with block as factor 

m2 <- brm(
  bf(zrt ~ CT * SC * BF +
       (1 + CT * SC * BF | subj_id) + (1 + CT * SC * BF | movie_id)
  ), 
  algorithm = "meanfield",
  family = asym_laplace(),
  iter = 20000, 
  init = 0.1,
  data = d
)
loo_m2 <- loo(m2)

comp <- loo_compare(loo_m1, loo_m2)
print(comp, digits = 2)


print(loo_m2)

pp_check(m2)


#' Model with block as numeric 

m3 <- brm(
  bf(zrt ~ CT * SC * BL +
       (1 + CT * SC * BL | subj_id) + (1 + CT * SC * BL | movie_id)
  ), 
  algorithm = "meanfield",
  family = asym_laplace(),
  iter = 25000, 
  init = 0.1,
  data = d
)
loo_m3 <- loo(m3)

comp <- loo_compare(loo_m3, loo_m2)
print(comp, digits = 2)

pp_check(m3)


print(loo_m3)


#' Test whether the kind of video is important within the surprise experiment.

m4 <- brm(
  bf(zrt ~ CT * BL +
       (1 + CT * BL | subj_id) + (1 | movie_id)
  ), 
  algorithm = "meanfield",
  family = exgaussian(),
  iter = 20000, 
  init = 0.1,
  data = d
)
loo_m4 <- loo(m4)


comp <- loo_compare(loo_m3, loo_m4)
print(comp, digits = 2)


#' Remove the three-way interaction

m5 <- brm(
  bf(zrt ~ CT * SC + CT * BL + SC * BL + 
       (1 + CT * SC + CT * BL + SC * BL | subj_id) + (1 | movie_id)
  ), 
  algorithm = "meanfield",
  family = brms::exgaussian(),
  iter = 20000, # Increase the number of iterations
  init = 0.1,
  data = d
)
loo_m5 <- loo(m5)

# Test of the three-way interaction
comp <- loo_compare(loo_m3, loo_m5)
print(comp, digits = 2)


#' Conditional plot.

mod <- m2

#' Three-way interaction
conditions <- make_conditions(mod, "SC")
c_eff <- conditional_effects(mod, "BF:CT", conditions=conditions) 
plot(c_eff, plot = FALSE)[[1]] +
  theme(legend.position = "bottom") +
  labs(
    y = "Reaction Times (standardized)",
    x = "Block of Trials"
  )


#' Plot of the raw data means

plot(density(log(surprise_cor_df$rt)))

surprise_cor_df$lrt <- log(surprise_cor_df$rt)

# Calculate the within-subject mean and standard error for each condition
subject_summary <- surprise_cor_df %>%
  group_by(subj_id, SC, CT, BL) %>%
  summarize(
    subj_mean = mean(lrt, na.rm = TRUE),
    .groups = 'drop'
  )

# Calculate the overall mean and within-subject standard error for each condition
plot_df <- subject_summary %>%
  group_by(SC, CT, BL) %>%
  summarize(
    m = mean(subj_mean, na.rm = TRUE),
    stderr = sd(subj_mean, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
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


#' Compute contrasts.

mod <- m3  # Use block as factor

# get the adjusted means
em <- emmeans (mod, ~ BF:CT | SC)
em

# get all possible contrasts
cont <- contrast(em, "tukey")
cont

# get the posterior draws from the contrasts
cont_posterior <- gather_emmeans_draws(cont)

# plot
ggplot(cont_posterior,
       aes(y = contrast, x = .value)) +
  stat_halfeye() +
  facet_wrap(~SC) +
  geom_vline(xintercept = 0, color = "red", lty = 2)


#' Control experiment

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


#' Control Model with block as numeric

c3 <- brm(
  bf(zrt ~ CT * BF +
       (1 + CT * BF | subj_id) + (1 + + CT * BF | movie_id)
  ), 
  algorithm = "meanfield",
  family = asym_laplace(),
  iter = 20000, 
  init = 0.1,
  data = dc
)
loo_c3 <- loo(c3)

print(loo_c3)

summary(c3)

mod_c <- c3

conditional_effects(mod_c, "BF:CT")

# Two-way interaction
c_eff <- conditional_effects(mod_c, "BL:CT") 
plot(c_eff, plot = FALSE)[[1]] +
  theme(legend.position = "bottom") +
  labs(
    y = "Reaction Times (standardized)",
    x = "Block of Trials"
  )


#' Model with no interaction

c4 <- brm(
  bf(zrt ~ CT + BL +
       (1 + CT + BL | subj_id) + (1 | movie_id)
  ), 
  algorithm = "meanfield",
  family = brms::exgaussian(),
  iter = 20000, 
  init = 0.1,
  data = dc
)
loo_c4 <- loo(c4)

#' Test of the interaction

comp <- loo_compare(loo_c3, loo_c4)
print(comp, digits = 2)


#' Accuracy modeling

# TODO

# Get the accuracy split by congruency.
flanker_accuracy <- get_flanker_accuracy(flanker_data, overall = FALSE)

flanker_accuracy |> 
  group_by(experiment, is_congruent_trial) |> 
  summarize(
    acc = mean(accuracy)
  )


message("\n20_variational_inference.R: done!")

# eof ----




