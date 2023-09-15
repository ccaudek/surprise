#' Surprise Experiment
#' Bayesian analysis of the control experiment.
#'
#' File: 35_bayes_control.R
#'
#' This file was last modified on "Mon May  7 10:30:51 2018"


library("lme4")
library("brms")
library("BEST")
library("effsize")

theme_set(theme_default())


# Remove NAs ---------------------------

clean_df <- final %>%
  dplyr::filter(!is.na(rtTukey) &
                  id != 11) # too slow with 100% correct performance!

clean_df$blockf <- factor(clean_df$block)

clean_df %>%
  group_by(is_congruent_trial) %>%
  summarise(
    mrt = mean(rtTukey, trim = 0.1, na.rm = TRUE),
    se = sqrt(var(rtTukey, na.rm = TRUE) / n())
  )
#   is_congruent_trial   mrt    se
# 1 Incongruent         591.  2.52
# 2 Congruent           545.  2.80


# Bayesian analysis with ex-gaussian distribution ---------------------------

fm1 <- lmer(interference_eff ~ block * is_clip_trial +
              (1 + is_clip_trial | subject_name),
            data = bysub_df
)

summary(fm1)
# Family: exgaussian 
# Links: mu = identity; sigma = identity; beta = identity 
# Formula: interference_eff ~ block * is_clip_trial + (1 + is_clip_trial | subject_name) 
# Data: bysub_df (Number of observations: 368) 
# Samples: 2 chains, each with iter = 2000; warmup = 1000; thin = 1; 
# total post-warmup samples = 2000
# ICs: LOO = NA; WAIC = NA; R2 = NA
# 
# Group-Level Effects: 
#   ~subject_name (Number of levels: 27) 
# Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
# sd(Intercept)                     53.64      8.01    39.89    71.45        711 1.00
# sd(is_clip_trialNo)               21.46      8.42     3.28    38.19        333 1.00
# cor(Intercept,is_clip_trialNo)    -0.21      0.28    -0.69     0.42       1713 1.00
# 
# Population-Level Effects: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
# Intercept                29.24     12.75     3.88    53.27        431 1.00
# block                     2.26      2.93    -3.29     8.21       1659 1.00
# is_clip_trialNo           0.55     11.47   -21.69    22.16       1403 1.00
# block:is_clip_trialNo    -3.60      4.05   -11.44     4.13       1468 1.00
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
# sigma    32.50      2.79    26.88    37.94        728 1.00
# beta     30.19      3.77    23.34    38.06       2000 1.00


mod1 <- brm(interference_eff ~ block * is_clip_trial +
              (1 + is_clip_trial | subject_name),
            family = exgaussian(),
            data = bysub_df,
            chains = 2,
            cores = 4
)





# One family that is especially suited to model reaction times is the exgaussian
# (‘exponentially modified Gaussian’) family, where β is the scale (inverse rate)
# of the exponential component, ξ is the mean of the Gaussian componenent, σ is
# the standard deviation of the Gaussian component. We parameterize μi = ξ + β so
# that the main predictor term equals the mean of the distribution.
mod1 <- brm(rtTukey ~ is_congruent_trial * block +
  (1 + is_congruent_trial | id) +
  (1 | video_id),
  family = exgaussian(),
  data = clean_df,
  chains = 2,
  cores = 4
)

summary(mod1)
# Family: exgaussian 
# Links: mu = identity; sigma = identity; beta = identity 
# Formula: rtTukey ~ is_congruent_trial * block + (1 + is_congruent_trial | id) + (1 | video_id) 
# Data: clean_df (Number of observations: 6614) 
# Samples: 2 chains, each with iter = 2000; warmup = 1000; thin = 1; 
# total post-warmup samples = 2000
# ICs: LOO = NA; WAIC = NA; R2 = NA
# 
# Group-Level Effects: 
#   ~id (Number of levels: 26) 
#                                            Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
# sd(Intercept)                                 59.92      9.02    45.03    79.78        560 1.00
# sd(is_congruent_trialCongruent)               32.03      5.70    22.74    44.78        629 1.00
# cor(Intercept,is_congruent_trialCongruent)    -0.00      0.22    -0.41     0.41        886 1.00
# 
# ~video_id (Number of levels: 9) 
#               Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
# sd(Intercept)     4.16      2.52     0.25     9.92        526 1.00
# 
# Population-Level Effects: 
#                                   Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
# Intercept                           632.65     12.29   609.08   656.99        444 1.00
# is_congruent_trialCongruent         -50.17      8.39   -66.36   -33.30        929 1.01
# block                               -11.21      1.49   -14.23    -8.36       2000 1.00
# is_congruent_trialCongruent:block    -2.58      2.17    -6.95     1.67       2000 1.00
# 
# Family Specific Parameters: 
#       Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
# sigma    56.73      1.34    54.14    59.37       2000 1.00
# beta    108.93      2.11   104.84   113.23       2000 1.00

marginal_effects(mod1, "block:is_congruent_trial")

# posterior_interval(mod1, prob = 0.95)

hypothesis(
  mod1, 
  "is_congruent_trialCongruent < 0", 
  class = "b", 
  group = "id"
)


# Leave-one-out cross-validation (LOO) ---------------------------

mod0 <- brm(rtTukey ~ block +
              (1 | id) +
              (1 | video_id),
            family = exgaussian(),
            data = clean_df,
            chains = 2,
            cores = 4
)

LOO(mod0, mod1)
#                LOOIC     SE
# mod0        81598.42 160.57
# mod1        80932.93 168.32
# mod0 - mod1   665.49  51.90

# the LOO information criterion for each model as well as the difference of the LOOs 
# each with its corresponding standard error are shown.
665.49 / 51.90
# [1] 12.82254



# Compute effect size on log-trasformed data ---------------------------

bysub <- clean_df %>%
  group_by(id, is_congruent_trial) %>%
  summarise(
    mrt = mean(log(rtTukey), trim = 0.0, na.rm = TRUE)
  )

bysub_con <- dplyr::filter(bysub, is_congruent_trial == "Congruent")
bysub_inc <- dplyr::filter(bysub, is_congruent_trial == "Incongruent")

t.test(bysub_inc$mrt, bysub_con$mrt, paired = TRUE)
# data:  bysub_inc$mrt and bysub_con$mrt
# t = 4.9498, df = 26, p-value = 3.837e-05
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.04529716 0.10963798
# sample estimates:
#   mean of the differences
# 0.07746757


# Compute effect size on log-trasformed RTs ---------------------------

# Eq. 9 from Lakens (2013)
mdiff <- mean(bysub_inc$mrt) - mean(bysub_con$mrt)
r <- cor(bysub_inc$mrt, bysub_con$mrt)

mdiff / sqrt(
  var(bysub_inc$mrt) + var(bysub_con$mrt) - 2 * r * sd(bysub_inc$mrt) * sd(bysub_con$mrt)
) * sqrt(2 * (1 - r))
# [1] 0.603886

mdiff / (
  (sd(bysub_inc$mrt) + sd(bysub_con$mrt)) / 2
  )
# [1] 0.6480844


# Compute effect size on raw data ---------------------------

bysubraw <- clean_df %>%
  group_by(id, is_congruent_trial) %>%
  summarise(
    mrt = mean(rtTukey, trim = 0.1, na.rm = TRUE)
  )

bysub_con <- dplyr::filter(bysubraw, is_congruent_trial == "Congruent")
bysub_inc <- dplyr::filter(bysubraw, is_congruent_trial == "Incongruent")

r <- cor(bysub_inc$mrt, bysub_con$mrt)

mdiff / sqrt(
  var(bysub_inc$mrt) + var(bysub_con$mrt) - 2 * r * sd(bysub_inc$mrt) * sd(bysub_con$mrt)
) * sqrt(2 * (1 - r))
# [1] 0.568167

mdiff <- mean(bysub_inc$mrt) - mean(bysub_con$mrt)
mdiff / (
  (sd(bysub_inc$mrt) + sd(bysub_con$mrt)) / 2
)
# [1] 0.5891324



# Lower and upper visual fields ---------------------------

clean_df$is_upper_visual_field <- ifelse(clean_df$delta_y > 0, 1, 0)


fm1 <- lmer(rtTukey ~ is_congruent_trial * delta_x * delta_y +
              (1 + is_congruent_trial | id) +
              (1 | video_id),
            data = clean_df
)

Anova(fm1)
#                                               Estimate Std. Error t value
# (Intercept)                                  6.993e+02  1.420e+01  49.252
# is_congruent_trialCongruent                 -2.493e+00  5.496e+00  -0.454
# delta_x                                     -1.302e-03  4.910e-03  -0.265
# delta_y                                      8.225e-03  4.918e-03   1.673
# is_congruent_trialCongruent:delta_x          2.160e-04  6.933e-03   0.031
# is_congruent_trialCongruent:delta_y         -1.484e-02  6.946e-03  -2.137
# delta_x:delta_y                             -2.440e-06  1.696e-05  -0.144
# is_congruent_trialCongruent:delta_x:delta_y  3.498e-05  2.398e-05   1.459


clean_df %>%
  group_by(is_upper_visual_field, is_surprise_clip, is_congruent_trial) %>%
  summarise(
    mrt = mean(rtTukey, na.rm = TRUE, trim = 0.1)
  )




##### procedure used to remove one outlying participant


bysub_df <- data.frame(id = bysub_con$id)

bysub_df$interference_eff <- 
  filter(bysubraw, is_congruent_trial == "Incongruent")$mrt - 
  filter(bysubraw, is_congruent_trial == "Congruent")$mrt

temp <- mydata %>%
  group_by(id) %>%
  summarise(
    acc = mean(correct, na.rm = TRUE),
    med_rt = median(rtTukey, na.rm = TRUE)
  )

bysub_df$acc <- temp$acc
bysub_df$med_rt <- temp$med_rt
bysub_df

mean(bysub_df$med_rt) + c(-1, 1) * 2.5 * sd(bysub_df$med_rt)


