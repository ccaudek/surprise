# Script name: 21_descr_stats_surpr.R
# Project: surprise with flanker task
# Script purpose: brms analysis of the congruenty effect
# @author: Corrado Caudek <corrado.caudek@unifi.it>
# Date Created: Tue Sep 12 08:39:19 2023
# Last Modified Date: Tue Sep 12 08:39:19 2023
#
# 👉 
# https://cran.r-project.org/web/packages/bayes4psy/vignettes/flanker.html


library("here")
suppressPackageStartupMessages(library("tidyverse")) 
library("brms")
library("cowplot")
library("reshape")
library("devtools")
library("mice")
library("tidybayes")
library("emmeans")
library("broom.mixed")
library("bayes4psy")
set.seed(123)

source(here("libraries", "fnct_surprise.R"))
source(here("libraries", "helpers.R"))


# The data have been created with the 01_, 10_, 11_ scripts in 
# the present directory.
mydata <- read_csv(here("data", "processed", "surprise_control_exps.csv"))

# Be sure that sequence of trial is in the right order
df <- mydata %>%
  dplyr::arrange(ID, date, time_of_day)

# check the number of blocks for each participant
n_subj <- length(unique(mydata$ID))
for (i in 1:n_subj) {
  one_subj <- mydata %>% 
    dplyr::filter(ID == i)
  cat("subject", i, ":", unique(one_subj$block), "\n")
}

# Remove blocks with less than 20 trials.
out <- df %>% 
  group_by(subj_name, block) %>% 
  summarise(
    ntrials_per_block = n()
  )
hist(out$ntrials_per_block)
df1 <- left_join(df, out)
df2 <- df1[df1$ntrials_per_block > 20, ]

# Remove subjects with less than 100 trials.
# add number of trials for each subject
df2 <- df2 %>% 
  group_by(subj_name) %>% 
  mutate(
    bysub_n = length(rtTukey)
  )
surprise_df <- df2 %>%
  dplyr::filter(bysub_n > 100)

# Remove block 5 and trial_number == 5
thedat <- surprise_df %>%
  dplyr::filter(
    trials_after_clip != "4" & block < 5
  )

# Number of subjects for each experiment.
thedat %>% 
  group_by(experiment) %>% 
  summarise(
    n = n_distinct(ID)
  )

# Check that all subjects completed all four blocks
n_subj <- length(unique(thedat$ID))
for (i in 1:n_subj) {
  one_subj <- thedat %>% 
    dplyr::filter(ID == i)
  cat("subject", i, ":", unique(thedat$block), "\n")
}

# data wrangling
thedat$experiment <- factor(thedat$experiment)
thedat$subj_name <- factor(thedat$ID)
thedat$resp <- factor(thedat$resp)
thedat$movie_id <- factor(thedat$movie_id)
thedat$is_surprise_clip <- factor(thedat$is_surprise_clip)
thedat$is_clip_trial <- factor(thedat$is_clip_trial)


# Better without imputation

# tobeimputed <- thedat %>% 
#   dplyr::select(
#     subj_name, experiment, ntrial, resp, is_congruent_trial, is_surprise_clip, 
#     correct, movie_id, trials_after_clip, rtTukey
#   )
# imp <- mice(tobeimputed, m = 1)
# newdat <- complete(imp)
# 
# thedat$rt <- newdat$rtTukey


# ------------------------------------------------------------------------------
# Descriptive stats

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
    rt = median(rtTukey, na.rm = TRUE)
  ) |> 
  as.data.frame()


#  -----------------------------------------------------------------------------
# Correct movie_id coding

part_to_remove <- "/Users/lab/Documents/surprise_project/surprise/stim/"
result <- sub(part_to_remove, "", thedat$movie_id)

part_to_remove <- "/Users/lab/Documents/magic3Definitivo/stim/"
result1 <- sub(part_to_remove, "", result)

part_to_remove <- "/Users/lab/Documents/magic_control_exp/stim/"
result2 <- sub(part_to_remove, "", result1)

part_to_remove <- "/Users/corrado/Dropbox/experiments/surprise_control_exp/stim/"
result3 <- sub(part_to_remove, "", result2)

part_to_remove <- "/Users/corrado/Documents/experiments/2019/surprise/stim/"
result4 <- sub(part_to_remove, "", result3)

part_to_remove <- "C:/Users/Miranda/Desktop/matlab_script/stim/"
result5 <- sub(part_to_remove, "", result4)

part_to_remove <- "/Users/lab/Documents/experiments/2019/surprise_5/stim/"
result6 <- sub(part_to_remove, "", result5)

part_to_remove <- "/Users/lab/Desktop/Surprise/stim/"
result7 <- sub(part_to_remove, "", result6)

thedat$movie_id <- result7
unique(thedat$movie_id)

# ------------------------------------------------------------------------------

# Select only correct trials.

dt_cor <- thedat |> 
  dplyr::filter(correct == 1)

nrow_total <- nrow(dt_cor)

# remove missing data on rt.
dt_cor <- dt_cor[!is.na(dt_cor$rt), ]

nrow_na_removed <- nrow(dt_cor)

# percent removed
(1 - nrow_na_removed / nrow_total) * 100
# [1] 0.07481794

nrow_total - nrow_na_removed  
# [1] 45
nrow_total
# [1] 60146

# ------------------------------------------------------------------------------

surprise_cor_df$RT <- surprise_cor_df$rt /1000

# Select trials of the surprise experiment 
surprise_cor_df <- dt_cor[dt_cor$experiment == "surprise", ]

# Select trials of the control experiment 
control_cor_df <- dt_cor[dt_cor$experiment == "control", ]

# ------------------------------------------------------------------------------

mydat <- surprise_cor_df |> 
  dplyr::filter(block == 1)

test_rt <- surprise_cor_df |> 
  dplyr::filter(
    block == 1 & 
      is_clip_trial == "Yes" & 
      is_congruent_trial == "Incongruent") |> 
  mutate(RT = rt / 1000) 
  
control_rt <- surprise_cor_df |> 
  dplyr::filter(
    block == 1 & 
      is_clip_trial == "Yes" & 
      is_congruent_trial == "Congruent") |> 
  mutate(RT = rt / 1000) 


# prior
uniform_prior <- b_prior(family="uniform", pars=c(0, 3))

# attach priors to relevant parameters
priors <- list(c("mu_m", uniform_prior))

control_rt$subject <- as.numeric(factor(as.character(control_rt$subj_id)))
test_rt$subject <- as.numeric(factor(as.character(test_rt$subj_id)))


# fit
rt_control_fit <- b_reaction_time(t=control_rt$RT,
                                  s=control_rt$subject,
                                  priors=priors,
                                  chains=4, iter=2000, warmup=1000)

rt_test_fit <- b_reaction_time(t=test_rt$RT,
                               s=test_rt$subject,
                               priors=priors,
                               chains=4, iter=2000, warmup=1000)


# plot trace
plot_trace(rt_control_fit)
plot_trace(rt_test_fit)

# check fits
plot(rt_control_fit)
plot(rt_test_fit)

# set rope (region of practical equivalence) interval to +/- 10ms
rope <- 0.01

# which mean is smaller/larger
rt_control_test <- compare_means(rt_control_fit, fit2=rt_test_fit, rope=rope)


# difference plot
plot_means_difference(rt_control_fit, fit2=rt_test_fit, rope=rope)

# visual comparsion of mean difference
plot_means(rt_control_fit, fit2=rt_test_fit)





# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# data wrangling
surprise_cor_df$blk <- factor(surprise_cor_df$block)
surprise_cor_df$zrt <- scale(surprise_cor_df$rt) |> as.numeric()
surprise_cor_df$CT <- surprise_cor_df$is_congruent_trial
surprise_cor_df$BF <- surprise_cor_df$blk
surprise_cor_df$BL <- surprise_cor_df$block
surprise_cor_df$SC <- surprise_cor_df$is_surprise_clip
surprise_cor_df$SC <- factor(surprise_cor_df$SC)

surprise_cor_df$CT <- factor(surprise_cor_df$CT)
surprise_cor_df$movie_id <- factor(surprise_cor_df$movie_id)

# I tried to remove first trial after the video-clip, but the results
# are worse.
d <- surprise_cor_df |>
  dplyr::select(rt, zrt, CT, SC, BL, BF, subj_id, movie_id) 

rio::export(d, "surprise_correct_data.csv")


# brms() analysis --------------------------------------------------------------


# Model without random slopes
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

# ------------------------------------------------------------------------------

# Target model
m3 <- brm(
  bf(zrt ~ CT * SC * BL +
       (1 + CT * SC * BL | subj_id) + (1 | movie_id)
  ), 
  algorithm = "meanfield",
  family = asym_laplace(),
  iter = 20000, # Increase the number of iterations
  data = d
)
loo_m3 <- loo(m3)

comp <- loo_compare(loo_m1, loo_m3)
print(comp, digits = 2)
#    elpd_diff se_diff
# m3    0.00      0.00 <- preferred model
# m1 -948.80     65.69
# The model with random slopes is better.

print(loo_m3)
# Computed from 1000 by 36032 log-likelihood matrix
# 
# Estimate    SE
# elpd_loo -34144.3 188.9
# p_loo      1237.9   9.2
# looic     68288.7 377.8
# ------
#   Monte Carlo SE of elpd_loo is 1.2.
# 
# Pareto k diagnostic values:
#   Count Pct.    Min. n_eff
# (-Inf, 0.5]   (good)     36012 99.9%   397       
# (0.5, 0.7]    (ok)          20  0.1%   698       
# (0.7, 1]      (bad)          0  0.0%   <NA>      
#   (1, Inf)    (very bad)     0  0.0%   <NA>      
#   
#   All Pareto k estimates are ok (k < 0.7).
# See help('pareto-k-diagnostic') for details.


# block as factor --------------------------------------------------------------

# best model
m3a <- brm(
  bf(zrt ~ CT * SC * BF +
       (1 + CT * SC * BF | subj_id) + (1 | movie_id)
  ), 
  algorithm = "meanfield",
  family = asym_laplace(),
  iter = 20000, # Increase the number of iterations
  data = d
)
loo_m3a <- loo(m3a)

comp <- loo_compare(loo_m3, loo_m3a)
print(comp, digits = 2)
#     elpd_diff se_diff
# m3a    0.00      0.00
# m3  -992.65     56.31
# The model with separate coefficients for each block is preferred.

# Interpretation: the random slopes are necessary for block.

pp_check(m3a)
# fit is ok.

# ------------------------------------------------------------------------------

# Test whether the kind of video is important within the surprise experiment.

m3b <- brm(
  bf(zrt ~ CT * BF +
       (1 + CT * BF | subj_id) + (1 | movie_id)
  ), 
  algorithm = "meanfield",
  family = asym_laplace(),
  iter = 20000, # Increase the number of iterations
  data = d
)
loo_m3b <- loo(m3b)


comp <- loo_compare(loo_m3a, loo_m3b)
print(comp, digits = 2)
#      elpd_diff se_diff 
# m3a     0.00      0.00 <- preferred model
# m3b -2381.65     84.03

# Interpretation: 
# The variable coding the kind of video within the surprise experiment 
# is necessary in the model.

# ------------------------------------------------------------------------------

# Remove the three-way interaction

m4 <- brm(
  bf(zrt ~ CT * SC + CT * BF + SC * BF + 
       (1 + CT * SC + CT * BF + SC * BF | subj_id) + (1 | movie_id)
  ), 
  algorithm = "meanfield",
  family = brms::asym_laplace(),
  iter = 20000, # Increase the number of iterations
  data = d
)
loo_m4 <- loo(m4)

# Test of the three-way interaction
comp <- loo_compare(loo_m3a, loo_m4)
print(comp, digits = 2)
#     elpd_diff se_diff
# m3a    0.00      0.00 <-- preferred model
# m4  -593.95     39.77

# ------------------------------------------------------------------------------

# Conditional plot.

mod <- m3a

conditional_effects(mod, "BF:CT")
conditional_effects(mod, "SC")

# Three-way interaction
conditions <- make_conditions(mod, "SC")
conditional_effects(mod, "BF:CT", conditions=conditions)


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

mod <- m3a

#get the adjusted means
em <- emmeans (mod, ~ BF:CT | SC)
em

#get all possible contrasts
cont <- contrast(em, "tukey")
cont
# SC = No Surprise:
#   contrast                          estimate lower.HPD upper.HPD
# BF1 Congruent - BF2 Congruent       0.2815   0.25319    0.3110
# BF1 Congruent - BF3 Congruent       0.3633   0.33370    0.3906
# BF1 Congruent - BF4 Congruent       0.4375   0.40349    0.4775
# BF1 Congruent - BF1 Incongruent     0.0933   0.07063    0.1215
# BF1 Congruent - BF2 Incongruent     0.1955   0.14066    0.2491
# BF1 Congruent - BF3 Incongruent     0.2661   0.21572    0.3221
# BF1 Congruent - BF4 Incongruent     0.3269   0.26240    0.3947
# BF2 Congruent - BF3 Congruent       0.0819   0.03849    0.1196
# BF2 Congruent - BF4 Congruent       0.1563   0.10692    0.2028
# BF2 Congruent - BF1 Incongruent    -0.1879  -0.22427   -0.1511
# BF2 Congruent - BF2 Incongruent    -0.0867  -0.12963   -0.0399
# BF2 Congruent - BF3 Incongruent    -0.0144  -0.07952    0.0454
# BF2 Congruent - BF4 Incongruent     0.0458  -0.03011    0.1142
# BF3 Congruent - BF4 Congruent       0.0749   0.03020    0.1193
# BF3 Congruent - BF1 Incongruent    -0.2700  -0.30638   -0.2310
# BF3 Congruent - BF2 Incongruent    -0.1675  -0.22396   -0.0997
# BF3 Congruent - BF3 Incongruent    -0.0953  -0.14053   -0.0466
# BF3 Congruent - BF4 Incongruent    -0.0363  -0.11133    0.0309
# BF4 Congruent - BF1 Incongruent    -0.3438  -0.39145   -0.3025
# BF4 Congruent - BF2 Incongruent    -0.2426  -0.30888   -0.1828
# BF4 Congruent - BF3 Incongruent    -0.1709  -0.23096   -0.1033
# BF4 Congruent - BF4 Incongruent    -0.1105  -0.16304   -0.0547
# BF1 Incongruent - BF2 Incongruent   0.1024   0.05715    0.1501
# BF1 Incongruent - BF3 Incongruent   0.1734   0.12903    0.2222
# BF1 Incongruent - BF4 Incongruent   0.2336   0.17757    0.2976
# BF2 Incongruent - BF3 Incongruent   0.0720   0.00413    0.1379
# BF2 Incongruent - BF4 Incongruent   0.1326   0.05796    0.2098
# BF3 Incongruent - BF4 Incongruent   0.0601  -0.01149    0.1383
# 
# SC = Surprise:
#   contrast                          estimate lower.HPD upper.HPD
# BF1 Congruent - BF2 Congruent       0.1883   0.13988    0.2314
# BF1 Congruent - BF3 Congruent       0.2663   0.22487    0.3148
# BF1 Congruent - BF4 Congruent       0.3503   0.28885    0.4083
# BF1 Congruent - BF1 Incongruent    -0.0650  -0.10301   -0.0233
# BF1 Congruent - BF2 Incongruent     0.1315   0.04469    0.2180
# BF1 Congruent - BF3 Incongruent     0.1691   0.08033    0.2470
# BF1 Congruent - BF4 Incongruent     0.2337   0.11573    0.3422
# BF2 Congruent - BF3 Congruent       0.0776   0.02145    0.1496
# BF2 Congruent - BF4 Congruent       0.1618   0.08601    0.2347
# BF2 Congruent - BF1 Incongruent    -0.2524  -0.31032   -0.1916
# BF2 Congruent - BF2 Incongruent    -0.0561  -0.13152    0.0180
# BF2 Congruent - BF3 Incongruent    -0.0188  -0.10845    0.0784
# BF2 Congruent - BF4 Incongruent     0.0472  -0.06911    0.1816
# BF3 Congruent - BF4 Congruent       0.0843   0.01326    0.1585
# BF3 Congruent - BF1 Incongruent    -0.3311  -0.39219   -0.2688
# BF3 Congruent - BF2 Incongruent    -0.1346  -0.22935   -0.0407
# BF3 Congruent - BF3 Incongruent    -0.0958  -0.16998   -0.0189
# BF3 Congruent - BF4 Incongruent    -0.0311  -0.14576    0.0840
# BF4 Congruent - BF1 Incongruent    -0.4153  -0.48230   -0.3389
# BF4 Congruent - BF2 Incongruent    -0.2201  -0.32004   -0.1017
# BF4 Congruent - BF3 Incongruent    -0.1831  -0.28677   -0.0805
# BF4 Congruent - BF4 Incongruent    -0.1140  -0.20561   -0.0268
# BF1 Incongruent - BF2 Incongruent   0.1953   0.12381    0.2785
# BF1 Incongruent - BF3 Incongruent   0.2336   0.15890    0.3043
# BF1 Incongruent - BF4 Incongruent   0.2983   0.20255    0.4027
# BF2 Incongruent - BF3 Incongruent   0.0393  -0.06391    0.1414
# BF2 Incongruent - BF4 Incongruent   0.1005  -0.02429    0.2261
# BF3 Incongruent - BF4 Incongruent   0.0666  -0.04402    0.1958
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

# control 

control_cor_df$blk <- factor(control_cor_df$block)
control_cor_df$zrt <- scale(control_cor_df$rt) |> as.numeric()
control_cor_df$CT <- control_cor_df$is_congruent_trial
control_cor_df$BF <- control_cor_df$blk
control_cor_df$BL <- control_cor_df$block

dc <- control_cor_df |> 
  dplyr::select(zrt, CT, BL, BF, subj_id, movie_id) 

dc$CT <- factor(dc$CT)
dc$movie_id <- factor(dc$movie_id)
dc$subj_id <- factor(dc$subj_id)


# ------------------------------------------------------------------------------

mc3a <- brm(
  bf(zrt ~ CT * BF +
       (1 + CT + BF | subj_id) + (1 | movie_id)
  ), 
  algorithm = "meanfield",
  family = brms::asym_laplace(),
  iter = 20000, # Increase the number of iterations
  data = dc
)
loo_mc3a <- loo(mc3a)

print(loo_mc3a)
# Computed from 1000 by 24069 log-likelihood matrix
# 
# Estimate    SE
# elpd_loo -25296.7 156.6
# p_loo       519.7   3.9
# looic     50593.5 313.2
# ------
#   Monte Carlo SE of elpd_loo is 0.8.
# 
# Pareto k diagnostic values:
#   Count Pct.    Min. n_eff
#   (-Inf, 0.5]   (good)     24059 100.0%  542       
#   (0.5, 0.7]    (ok)          10   0.0%  840       
#   (0.7, 1]      (bad)          0   0.0%  <NA>      
#   (1, Inf)      (very bad)     0   0.0%  <NA>      
#   
#   All Pareto k estimates are ok (k < 0.7).
# See help('pareto-k-diagnostic') for details.    

summary(mc3a)
# Population-Level Effects: 
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept            -0.50      0.02    -0.53    -0.46 1.00     1033      874
# CTIncongruent         0.14      0.01     0.12     0.17 1.00      723      907
# BF2                  -0.11      0.02    -0.15    -0.08 1.00      989      934
# BF3                  -0.18      0.02    -0.21    -0.15 1.00     1014      938
# BF4                  -0.24      0.02    -0.28    -0.21 1.00      978      884
# CTIncongruent:BF2     0.06      0.02     0.02     0.10 1.00      932      961
# CTIncongruent:BF3     0.06      0.02     0.02     0.10 1.00      870      909
# CTIncongruent:BF4     0.08      0.02     0.04     0.12 1.00      953      907

mod_c <- mc3a

conditional_effects(mod_c, "BF:CT")


# ------------------------------------------------------------------------------

# No interaction
mc3b <- brm(
  bf(zrt ~ CT + BF +
       (1 + CT + BF | subj_id) + (1 | movie_id)
  ), 
  algorithm = "meanfield",
  family = brms::asym_laplace(),
  iter = 20000, # Increase the number of iterations
  data = dc
)
loo_mc3b <- loo(mc3b)

# Test of the interaction
comp <- loo_compare(loo_mc3a, loo_mc3b)
print(comp, digits = 2)


# eof ----





# Only first block

blk1_dat <- surprise_cor_df |> 
  dplyr::filter(BL == 1) |> 
  dplyr::select(rt, zrt, rtTukey, CT, SC, subj_id, movie_id, 
                is_surprise_clip) 

blk1_dat$CT <- factor(blk1_dat$CT)
blk1_dat$movie_id <- factor(blk1_dat$movie_id)

m_blk1 <- brm(
  bf(rt ~ CT * SC +
       (1 + CT * SC | subj_id) + (1 | movie_id)
  ), 
  algorithm = "meanfield",
  # backend = "cmdstanr", 
  family = shifted_lognormal(),
  # iter = 20000, # Increase the number of iterations
  data = blk1_dat
)
pp_check(m_blk1)

conditional_effects(m_blk1, "SC:CT")

# Step 1: Get the conditional effects
ce <- conditional_effects(mod_c, "BF:CT")

# Step 2: Create a base plot
base_plot <- plot(ce)

plot_min_value <- -0.4
plot_max_value <-  0.3

# Step 3: Customize the plot
custom_plot <- base_plot +
  # labs(y = "Reaction Times (ms)", x = "Block") +
  theme(axis.title.y = element_text(size = 14, face = "bold")) + 
  scale_y_continuous(breaks = seq(plot_min_value, plot_max_value, by = 0.1)) # replace min_value, max_value and interval with appropriate values

# Display the plot
print(custom_plot)


















# Congruency effect ----

# Computing the size of the congruency effect for each subject 
bysub_df <- dt_cor %>% 
  group_by(
    subj_name, experiment, block, trials_after_clip, is_congruent_trial
  ) %>% 
  summarise(
    m_rt = mean(rt, na.rm = TRUE, trim = 0.1)
  )

# Spread is_congruent_trial on two columns.
bysub_df_wide <- bysub_df %>%
  pivot_wider(names_from = is_congruent_trial, values_from = m_rt)

# Impute missing values.
imp <- mice(bysub_df_wide, m = 1)
bysub_df_wide <- complete(imp)

bysub_df_wide <- bysub_df_wide %>% 
  mutate(
    congruency_eff = Incongruent - Congruent
  ) %>% 
  dplyr::select(-c(Congruent, Incongruent))

hist(bysub_df_wide$congruency_eff)

plot_df <- bysub_df_wide %>% 
  group_by(experiment, block, trials_after_clip) %>% 
  summarise(
    y = mean(congruency_eff, na.rm = TRUE),
    stderr = sqrt(var(congruency_eff) / n())
  )
plot_df$lower <- plot_df$y - plot_df$stderr
plot_df$upper <- plot_df$y + plot_df$stderr

pd <- position_dodge(0.2)

ggplot(plot_df, aes(x=trials_after_clip, y=y)) +
  geom_pointrange(aes(ymin=lower, ymax=upper), lwd=1.05, position=pd) +
  geom_line(position=pd, lwd=1.05) +
  geom_point(position=pd, size=5) +
  facet_grid(block ~experiment) +
  # ylim(-14, 65) +
  xlab("Trials after the video clip") +
  ylab("Congruency Effect (ms)") +
  # ggtitle(title)
  geom_hline(yintercept=0, lty=2) +
  theme(legend.position="bottom")



# Collapse trials after the clip: 0 vs. others

dt_cor$is_just_after_clip <- ifelse(
  dt_cor$trials_after_clip == 0, "yes", "no"
) %>% 
  as.factor()

bysub_df <- dt_cor %>% 
  group_by(
    subj_name, experiment, block, is_just_after_clip, is_congruent_trial
  ) %>% 
  summarise(
    m_rt = mean(rt, na.rm = TRUE, trim = 0.1)
  )

# Spread is_congruent_trial on two columns.
bysub_df_wide <- bysub_df %>%
  pivot_wider(names_from = is_congruent_trial, values_from = m_rt)

# Impute missing values.
imp <- mice(bysub_df_wide, m = 1)
bysub_df_wide <- complete(imp)

bysub_df_wide <- bysub_df_wide %>% 
  mutate(
    congruency_eff = Incongruent - Congruent
  ) %>% 
  dplyr::select(-c(Congruent, Incongruent))

hist(bysub_df_wide$congruency_eff)

plot_df <- bysub_df_wide %>% 
  group_by(experiment, block, is_just_after_clip) %>% 
  summarise(
    y = mean(congruency_eff, na.rm = TRUE),
    stderr = sqrt(var(congruency_eff) / n())
  )
plot_df$lower <- plot_df$y - plot_df$stderr
plot_df$upper <- plot_df$y + plot_df$stderr

pd <- position_dodge(0.2)

ggplot(plot_df, aes(x=block, y=y)) +
  geom_pointrange(aes(ymin=lower, ymax=upper), lwd=1.05, position=pd) +
  geom_line(position=pd, lwd=1.05) +
  geom_point(position=pd, size=5) +
  facet_grid(is_just_after_clip ~ experiment) +
  # ylim(-14, 65) +
  xlab("Block") +
  ylab("Congruency Effect (ms)") +
  # ggtitle(title)
  geom_hline(yintercept=0, lty=2) +
  theme(legend.position="bottom")


# Collapse all trials after the clip -------------------------------------------

bysub_df <- dt_cor %>% 
  group_by(
    subj_name, experiment, block, is_congruent_trial
  ) %>% 
  summarise(
    m_rt = mean(rt, na.rm = TRUE, trim = 0.1)
  ) |> 
  ungroup()

# Spread is_congruent_trial on two columns.
bysub_df_wide <- bysub_df %>%
  pivot_wider(names_from = is_congruent_trial, values_from = m_rt)

# Impute missing values.
imp <- mice(bysub_df_wide, m = 1)
bysub_df_wide <- complete(imp)

bysub_df_wide <- bysub_df_wide %>% 
  mutate(
    congruency_eff = Incongruent - Congruent
  ) %>% 
  dplyr::select(-c(Congruent, Incongruent))

hist(bysub_df_wide$congruency_eff)

plot_df <- bysub_df_wide %>% 
  group_by(experiment, block) %>% 
  summarise(
    y = mean(congruency_eff, na.rm = TRUE),
    stderr = sqrt(var(congruency_eff) / n())
  )
plot_df$lower <- plot_df$y - plot_df$stderr
plot_df$upper <- plot_df$y + plot_df$stderr

pd <- position_dodge(0.2)

plot_df |>
  ggplot(aes(x = block, y = y, color=experiment)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), lwd = 1.05, position = pd) +
  geom_line(position = pd, lwd = 1.05) +
  geom_point(position = pd, size = 5) +
  #facet_grid(~experiment) +
  # ylim(-14, 65) +
  xlab("Block") +
  ylab("Congruency Effect (ms)") +
  # ggtitle(title)
  geom_hline(yintercept = 0, lty = 2) +
  theme(legend.position = "bottom")




# Surprise effect on RTs, for each experiment
dt_cor %>% 
  group_by(experiment, is_surprise_clip) %>% 
  summarise(
    mrt = mean(rt, trim = 0.1)
  )

by_subj_cor <- dt_cor %>% 
  group_by(subj_name, experiment, is_surprise_clip) %>% 
  summarize(
    mrt = mean(rt, trim = 0.1, na.rm = TRUE)
  )



# brms analysis --------------------------------------------------------------

plot_df$user_id <- 1:nrow(plot_df)

m1 <- brm(
  bf(y ~ experiment * block  +
       (1 + block | user_id)
  ), 
  family = student(),
  # algorithm = "meanfield",
  data = plot_df, 
  cores = 4, 
  # iter = 20000,
  backend = "cmdstanr"
)

pp_check(m1)
print(m1, 4)

conditional_effects(
  m1, 
  "block",
  prob = 0.95
)

con_eff_surprise <- plot_df[plot_df$experiment == "surprise", ]
con_eff_surprise$zy <- scale(con_eff_surprise$y)

con_eff_surprise$blk <- factor(con_eff_surprise$block)

m2 <- brm(
  bf(zy ~ blk  +
       (1 + blk | user_id)
  ), 
  family = student(),
  # algorithm = "meanfield",
  data = con_eff_surprise, 
  cores = 4, 
  iter = 20000,
  backend = "cmdstanr"
)

pp_check(m2)
print(m2, 4)


by_subj_cor$experiment <- factor(by_subj_cor$experiment)
by_subj_cor$is_surprise_clip <- factor(by_subj_cor$is_surprise_clip)

mod_1 <- brm(
  bf(mrt ~ experiment  +
       (1 | subj_name)
  ), 
  family = gaussian(),
  # algorithm = "meanfield",
  data = by_subj_cor, 
  cores = 6, 
  iter = 20000,
  backend = "cmdstanr"
)

pp_check(mod_1)
print(mod_1, 4)
#                    Estimate Est.Error l-95% CI u-95% CI   Rhat Bulk_ESS Tail_ESS
# Intercept          612.2462   15.0592 582.6349 641.3495 1.0151      560     1115
# experimentsurprise  78.6705   19.4706  41.1843 116.1662 1.0108      434      899

conditional_effects(
  mod_1, 
  "experiment",
  prob = 0.95
)

# The presence of surprising video-clips (in half of the trials) made the overall
# performance slower of 78 ms with respect to the condition in which all video-clips
# presented non-surprising events.


bysub_df_wide$blk <- factor(bysub_df_wide$block)
fit_1 <- brm(
  bf(congruency_eff ~ 0 + blk), 
  family = gaussian(),
  data = bysub_df_wide %>% filter(experiment == "surprise" & is_just_after_clip == "no"), 
  cores = 6, 
  backend = "cmdstanr"
)
summary(fit_1)
# Population-Level Effects: 
#      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# blk1   -14.90      6.54   -27.70    -2.36 1.00     5281     2945
# blk2    -4.04      6.68   -16.61     9.13 1.00     4062     3002
# blk3     9.55      6.64    -3.49    22.58 1.00     5383     3000
# blk4    18.14      6.80     5.04    31.25 1.00     4700     3210

conditional_effects(
  fit_1, 
  "blk",
  prob = 0.95
)




fit_2 <- brm(
  bf(congruency_eff ~ 0 + blk), 
  family = gaussian(),
  data = bysub_df_wide %>% filter(experiment == "surprise" & is_just_after_clip == "yes"), 
  cores = 6, 
  backend = "cmdstanr"
)
summary(fit_2)
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# blk1   -18.38      7.72   -33.19    -3.03 1.00     4834     3244
# blk2    -9.07      8.03   -24.83     6.44 1.00     5342     3194
# blk3     4.23      7.95   -11.45    19.61 1.00     5609     2739
# blk4    10.11      8.08    -5.51    26.16 1.00     4927     3012

conditional_effects(
  fit_2, 
  "blk",
  prob = 0.95
)


fit_3 <- brm(
  bf(congruency_eff ~ 0 + blk), 
  family = gaussian(),
  data = bysub_df_wide %>% filter(experiment == "surprise"), 
  cores = 6, 
  backend = "cmdstanr"
)
summary(fit_3)
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# blk1   -17.07      5.26   -27.41    -6.95 1.00     5467     3070
# blk2    -6.50      5.20   -16.75     3.56 1.00     6010     3284
# blk3     7.03      5.18    -2.83    17.15 1.00     4564     3380
# blk4    13.71      5.29     3.29    24.15 1.00     5092     3249

conditional_effects(
  fit_3, 
  "blk",
  prob = 0.95
)


fit_5 <- brm(
  bf(congruency_eff ~ 0 + blk), 
  family = skew_normal(),
  data = bysub_df_wide %>% filter(experiment == "control" & is_just_after_clip == "no"), 
  cores = 6, 
  backend = "cmdstanr"
)
summary(fit_5)
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# blk1     7.76      7.79    -7.22    23.33 1.00     5793     2894
# blk2    34.30      7.67    19.93    49.62 1.00     5760     3268
# blk3    31.25      7.89    15.51    46.72 1.00     5544     2924
# blk4    33.24      7.91    17.76    48.60 1.00     5428     3080

conditional_effects(
  fit_5, 
  "blk",
  prob = 0.95
)


fit_6 <- brm(
  bf(congruency_eff ~ 0 + blk), 
  family = gaussian(),
  data = bysub_df_wide %>% filter(experiment == "control" & is_just_after_clip == "yes"), 
  cores = 6, 
  backend = "cmdstanr"
)
summary(fit_6)
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# blk1    -0.45      8.42   -16.40    16.27 1.00     5393     3165
# blk2    11.35      8.43    -5.69    27.86 1.00     4844     3043
# blk3     7.97      8.59    -8.86    25.04 1.00     4612     3456
# blk4    12.26      8.53    -4.32    28.61 1.00     4757     3085


marginal_effects(
  fit_6, 
  "blk"
)


fit_7 <- brm(
  bf(congruency_eff ~ 0 + blk), 
  family = gaussian(),
  data = bysub_df_wide %>% filter(experiment == "control"), 
  cores = 6, 
  backend = "cmdstanr"
)
summary(fit_7)
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# blk1     4.71      5.82    -6.68    16.39 1.00     5566     3052
# blk2    22.15      5.84    10.64    33.79 1.00     5678     2931
# blk3    19.21      5.94     7.33    30.77 1.00     5534     2888
# blk4    23.88      6.08    12.17    35.91 1.00     4478     3233




# ---------------------- end




# melsm: mixed effects location scale model -------------------------------


# me: mixed effects model
fit_1 <- brm(
  bf(
    RT ~ 1 + (1|c|ID),
    sigma ~ 1 + (1|c|ID)
  ),
  data = dt_cor %>% filter(is_congruent_trial == "Incongruent", experiment == "control"),
  inits = 0,
  cores = 4,
  chains = 2,
  iter = 2000,
  warmup = 1000,
  backend = "cmdstanr"
)


# remove scale model
fit_2 <- brm(
  bf(RT ~ 1 + (1|c|ID)), 
  data = mydata %>% filter(is_congruent_trial == "Incongruent", experiment == "control"),
  inits = 0,
  cores = 4,
  chains = 2,
  iter = 2000,
  warmup = 1000
)

# Compare with WAIC.
incong_waic <- waic(fit_1, fit_2)




form <- brmsformula(
  RT ~ is_congruent_trial + (is_congruent_trial |c| ID), 
  sigma  ~ is_congruent_trial + (is_congruent_trial |c| ID)
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
  set_prior("normal(0,0.25)",
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
  set_prior("normal(0,1)", 
            class = "sd", 
            coef = "Intercept", 
            group = "ID", 
            dpar = "sigma"), 
  # fixed effect priors
  set_prior("normal(0,5)", 
            class = "b"),
  set_prior("normal(0,5)", 
            class = "b", 
            dpar = "sigma"),
  set_prior("normal(0,5)", 
            class = "Intercept"),
  set_prior("normal(0,5)", 
            class = "Intercept", 
            dpar = "sigma"))




fit_flank <- brm(
  form, 
  data = dt_cor %>% filter(experiment == "control"), 
  inits = 0, 
  cores = 4, 
  chains = 2, 
  iter = 1500, 
  warmup = 1000,  
  prior = priors
)



dt_cor$trials_after <- factor(dt_cor$trials_after_clip)


form <- brmsformula(
  RT ~ is_congruent_trial * trials_after + (is_congruent_trial + trials_after |c| ID)
)


priors_1 <- c(
  set_prior("student_t(5, 0, 2.0)", class = "b"),
  set_prior("student_t(5, 0, 2.0)", class = "Intercept"),
  set_prior("cauchy(0, 2.0)", class = "sd")
)


fit3_cntr <- brm(
  form, 
  data = dt_cor %>% filter(experiment == "control", block == 1), 
  family = exgaussian(),
  inits = 0, 
  cores = 4, 
  chains = 2, 
  iter = 1500, 
  warmup = 1000,  
  prior = priors_1
)


pp_check(fit3_cntr)
print(fit3_cntr, 4)
pc <- marginal_effects(
  fit3_cntr, 
  "trials_after:is_congruent_trial")

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
