# Script name: 20_variational_inference.R
# Project: surprise with flanker task
# Script purpose: brms analysis of the congruenty effect
# @author: Corrado Caudek <corrado.caudek@unifi.it>
# Date Created: Tue Sep 12 08:39:19 2023
# Last Modified Date: Fri Sep 15 08:10:18 2023
#
# 👉 Variational Inference with brms for the raw data

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

# source(here("libraries", "fnct_surprise.R"))
# source(here("libraries", "helpers.R"))
source(here("libraries", "functions.R"))


# ------------------------------------------------------------------------------
# Import complete data set
# ------------------------------------------------------------------------------

# Import the data that have been created by the previous scripts.
# The data have been created with the 01_, 10_, 11_ scripts in 
# the present directory.
data <- get_data()


# ------------------------------------------------------------------------------
# Tidy data
# ------------------------------------------------------------------------------

# tidy the data frame
data_tidy <- tidy_flanker(data)

# Perform some participants' flanker checks.
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

# Check number of subjecta by condition.
flanker_data |>
  group_by(experiment, is_surprise_clip) |>
  summarize(
    n = n_distinct(subj_id)
  ) 
#   experiment is_surprise_clip     n
# 1 control    No                  81
# 2 surprise   No Surprise        120
# 3 surprise   Surprise           120


# ------------------------------------------------------------------------------
# Accuracy modeling
# ------------------------------------------------------------------------------

# TODO

# Get the accuracy split by congruency.
flanker_accuracy <- get_flanker_accuracy(flanker_data, overall = FALSE)

flanker_accuracy |> 
  group_by(experiment, is_congruent_trial) |> 
  summarize(
    acc = mean(accuracy)
  )
# experiment is_congruent_trial   acc
# <fct>      <chr>              <dbl>
# 1 control    Congruent          0.967
# 2 control    Incongruent        0.954
# 3 surprise   Congruent          0.963
# 4 surprise   Incongruent        0.961

surprise_acc_df <- flanker_accuracy |> 
  dplyr::filter(
    experiment == "surprise"
  )

mod1_acc <- brm(
  accuracy ~ is_congruent_trial +
    (1 + is_congruent_trial | subj_id),
  family = zero_one_inflated_beta(),
  data = surprise_acc_df, 
  backend = "cmdstanr"
)

pp_check(mod1_acc)

summary(mod1_acc)
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS
# Intercept                         3.54      0.09     3.36     3.73 1.01      501
# is_congruent_trialIncongruent    -0.11      0.11    -0.34     0.11 1.00      732

# get the congruency effect
flanker_congruency_accuracy <- get_flanker_congruency_accuracy(flanker_data)

# add the mean accuracy
mean_accuracy <- flanker_congruency_accuracy %>% 
  group_by(subj_id) %>% 
  summarise(mean_accuracy = mean(accuracy))


mean(mean_accuracy$mean_accuracy, na.rm = TRUE)

# accuracy 
thedat |> 
  group_by(experiment, is_surprise_clip, is_congruent_trial) |> 
  summarize(
    is_correct = mean(correct)
  )

# RTs
thedat |> 
  group_by(experiment, is_surprise_clip, block, is_congruent_trial) |> 
  summarize(
    rt = median(rt, na.rm = TRUE)
  ) |> 
  as.data.frame()


# ------------------------------------------------------------------------------
# Select correct trials only
# ------------------------------------------------------------------------------

dt_cor <- flanker_data |> 
  dplyr::filter(correct == 1)

nrow_total <- nrow(dt_cor)

# remove missing data on rt.
dt_cor <- dt_cor[!is.na(dt_cor$rt), ]

nrow_na_removed <- nrow(dt_cor)

# percent removed
(1 - nrow_na_removed / nrow_total) * 100
# 0.07513148

nrow_total - nrow_na_removed  
# [1] 45
nrow_total
# [1] 59895


# ------------------------------------------------------------------------------
# Select correct trials by experiment
# ------------------------------------------------------------------------------

# Select correct trials of the surprise experiment 
surprise_cor_df <- dt_cor[dt_cor$experiment == "surprise", ]

# Select correct trials of the control experiment 
control_cor_df <- dt_cor[dt_cor$experiment == "control", ]


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
       (1 | subj_id) + (1 | movie_id)
  ), 
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
       (1 + CT * SC * BF | subj_id) + (1 | movie_id)
  ), 
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
       (1 + CT * SC * BL | subj_id) + (1 | movie_id)
  ), 
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
       (1 + CT * BL | subj_id) + (1 | movie_id)
  ), 
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
       (1 + CT * SC + CT * BL + SC * BL | subj_id) + (1 | movie_id)
  ), 
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

mod <- m3

# Three-way interaction
conditions <- make_conditions(mod, "SC")
c_eff <- conditional_effects(mod, "BL:CT", conditions=conditions) 
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

# ------------------------------------------------------------------------------

# Compute contrasts.

mod <- m2  # Use block as factor

#get the adjusted means
em <- emmeans (mod, ~ BF:CT | SC)
em

#get all possible contrasts
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

#get the posterior draws from the contrasts
cont_posterior <- gather_emmeans_draws(cont)

#plot
ggplot(cont_posterior,
       aes(y = contrast, x = .value)) +
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
       (1 + CT * BL | subj_id) + (1 | movie_id)
  ), 
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
       (1 + CT + BL | subj_id) + (1 | movie_id)
  ), 
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




