#' Surprise Experiment
#' Bayesian analysis.
#'
#' File: 40_bayesian_lme.R
#'
#' This file was last modified on "Mon May  7 15:58:26 2018"


library("lme4")
library("brms")
library("BEST")


source("00_import_data.R")
source("10_tidying_data.R")


theme_set(theme_default())



# Effect of number of trials after the clip -------------------------------


final$lrt <- log(rt)

bysub_df <- final %>% 
  group_by(
    subject_name, trials_after_clip
    ) %>% 
  summarise(
    y = mean(lrt, na.rm = TRUE)
  )

# Randomly sample 10 id's
set.seed(111111111)
rand <- sample(unique(bysub_df$subject_name), 12)

# Plot those 10 id's individually over time. 
ggplot(
  data = bysub_df[bysub_df$subject_name %in% rand,], 
  aes(x = trials_after_clip, 
      y = y)
  ) +
  geom_point() +
  ylim(6, 6.8) +
  #geom_smooth(method = 'lm', formula = y ~ x) +
  facet_wrap("subject_name", nrow = 3)



# Remove NAs ---------------------------

clean_df <- final %>%
  dplyr::filter(!is.na(rtTukey) &
                  id != 11) # too slow with 100% correct performance!

clean_df$blockf <- factor(clean_df$block)

clean_df %>%
  group_by(is_surprise_clip, is_congruent_trial) %>%
  summarise(
    mrt = mean(rtTukey, trim = 0.1, na.rm = TRUE),
    se = sqrt(var(rtTukey, na.rm = TRUE) / n())
  )
#   is_surprise_clip is_congruent_trial   mrt    se
# 1 Surprise         Incongruent         673.  2.41
# 2 Surprise         Congruent           670.  2.71
# 3 No Surprise      Incongruent         671.  2.74
# 4 No Surprise      Congruent           668.  2.99


# Bayesian analysis with ex-gaussian distribution ---------------------------

# One family that is especially suited to model reaction times is the exgaussian
# (‘exponentially modified Gaussian’) family, where β is the scale (inverse rate)
# of the exponential component, ξ is the mean of the Gaussian componenent, σ is
# the standard deviation of the Gaussian component. We parameterize μi = ξ + β so
# that the main predictor term equals the mean of the distribution.

suprise_df <- clean_df %>%
  dplyr::filter(is_surprise_clip == "No Surprise")

model1 <- brm(rtTukey ~ is_congruent_trial * block +
              (1 + is_congruent_trial | id) +
              (1 | video_id),
            family = exgaussian(),
            data = suprise_df,
            chains = 2,
            cores = 4
)

summary(model1)
# Family: exgaussian 
# Links: mu = identity; sigma = identity; beta = identity 
# Formula: rtTukey ~ is_congruent_trial * block + (1 + is_congruent_trial | id) + (1 | video_id) 
# Data: suprise_df (Number of observations: 17503) 
# Samples: 2 chains, each with iter = 2000; warmup = 1000; thin = 1; 
# total post-warmup samples = 2000
# ICs: LOO = NA; WAIC = NA; R2 = NA
#   
# Group-Level Effects: 
#   ~id (Number of levels: 119) 
#                                            Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
# sd(Intercept)                                123.79      8.10   108.79   141.02        202 1.01
# sd(is_congruent_trialCongruent)               48.08      3.73    41.33    55.67        922 1.00
# cor(Intercept,is_congruent_trialCongruent)     0.26      0.09     0.07     0.44        819 1.00
# 
# ~video_id (Number of levels: 38) 
#               Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
# sd(Intercept)     7.54      1.78     4.36    11.37        757 1.00
# 
# Population-Level Effects: 
#                                   Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
# Intercept                           721.33     11.68   697.04   743.37         27 1.07
# is_congruent_trialCongruent         -12.92      6.72   -25.74     0.20       1011 1.00
# block                                -3.55      1.24    -6.00    -1.11       2000 1.00
# is_congruent_trialCongruent:block    -3.81      1.79    -7.17    -0.27       2000 1.00
# 
# Family Specific Parameters: 
#       Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
# sigma    71.93      1.32    69.44    74.57       2000 1.00
# beta    162.58      2.03   158.53   166.65       2000 1.00



marginal_effects(model1)









#######################################################################################

block1_congruent <- final %>% 
  dplyr::filter(# block == 1 & 
                  is_congruent_trial == "Congruent" &
                  rt > 200)

block1_congruent$blockf <- factor(block1_congruent$block)


# One family that is especially suited to model reaction times is the exgaussian 
# (‘exponentially modified Gaussian’) family, where β is the scale (inverse rate) 
# of the exponential component, ξ is the mean of the Gaussian componenent, σ is 
# the standard deviation of the Gaussian component. We parameterize μi = ξ + β so 
# that the main predictor term equals the mean of the distribution.
mod1 <- brm(rt ~ is_surprise_clip * blockf +
              (1 + is_surprise_clip + blockf | id) +
              (1 | video_id),
            family = exgaussian(),
            data = block1_congruent, 
            chains = 2, 
            cores = 2)

summary(mod1)
# Family: exgaussian 
# Links: mu = identity; sigma = identity; beta = identity 
# Formula: rt ~ is_surprise_clip * blockf + (1 + is_surprise_clip + blockf | id) + (1 | video_id) 
# Data: block1_congruent (Number of observations: 16728) 
# Samples: 2 chains, each with iter = 2000; warmup = 1000; thin = 1; 
# total post-warmup samples = 2000
# ICs: LOO = NA; WAIC = NA; R2 = NA
# 
# Group-Level Effects: 
#   ~id (Number of levels: 128) 
#                                           Estimate Est.Error l-95% CI u-95% CI
# sd(Intercept)                               164.19     10.71   143.86   187.54
# sd(is_surprise_clipNoSurprise)                4.55      3.16     0.16    11.39
# sd(blockf2)                                  60.46      5.40    50.44    71.68
# sd(blockf3)                                  83.26      6.58    71.78    97.17
# sd(blockf4)                                  84.09      6.80    71.35    97.81
# cor(Intercept,is_surprise_clipNoSurprise)    -0.22      0.34    -0.79     0.58
# cor(Intercept,blockf2)                       -0.20      0.11    -0.39     0.02
# cor(is_surprise_clipNoSurprise,blockf2)      -0.00      0.33    -0.65     0.64
# cor(Intercept,blockf3)                       -0.33      0.09    -0.50    -0.16
# cor(is_surprise_clipNoSurprise,blockf3)       0.16      0.34    -0.63     0.73
# cor(blockf2,blockf3)                          0.75      0.06     0.61     0.85
# cor(Intercept,blockf4)                       -0.41      0.08    -0.56    -0.24
# cor(is_surprise_clipNoSurprise,blockf4)       0.28      0.40    -0.63     0.86
# cor(blockf2,blockf4)                          0.40      0.10     0.18     0.58
# cor(blockf3,blockf4)                          0.79      0.05     0.67     0.87
#                                           Eff.Sample Rhat
# sd(Intercept)                                    276 1.00
# sd(is_surprise_clipNoSurprise)                   138 1.03
# sd(blockf2)                                      984 1.00
# sd(blockf3)                                      724 1.00
# sd(blockf4)                                      817 1.00
# cor(Intercept,is_surprise_clipNoSurprise)       2000 1.00
# cor(Intercept,blockf2)                           954 1.00
# cor(is_surprise_clipNoSurprise,blockf2)           71 1.01
# cor(Intercept,blockf3)                           987 1.00
# cor(is_surprise_clipNoSurprise,blockf3)           56 1.06
# cor(blockf2,blockf3)                             690 1.00
# cor(Intercept,blockf4)                          1028 1.00
# cor(is_surprise_clipNoSurprise,blockf4)           26 1.21
# cor(blockf2,blockf4)                             943 1.00
# cor(blockf3,blockf4)                             937 1.00
# 
# ~video_id (Number of levels: 38) 
#               Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
# sd(Intercept)     2.10      1.42     0.07     5.03        881 1.00
# 
# Population-Level Effects: 
#                                    Estimate Est.Error l-95% CI u-95% CI Eff.Sample
# Intercept                            725.75     14.18   695.77   751.68        159
# is_surprise_clipNoSurprise            -8.78      4.26   -17.18    -0.46       1326
# blockf2                              -12.64      6.94   -25.69     0.81        946
# blockf3                              -19.90      8.53   -36.57    -2.96        641
# blockf4                              -27.61      8.50   -43.87   -10.78        672
# is_surprise_clipNoSurprise:blockf2    11.13      5.91    -0.98    22.44       1521
# is_surprise_clipNoSurprise:blockf3    12.75      5.82     1.31    23.72       1538
# is_surprise_clipNoSurprise:blockf4    12.92      5.94     1.51    24.93       1521
#                                    Rhat
# Intercept                          1.01
# is_surprise_clipNoSurprise         1.00
# blockf2                            1.00
# blockf3                            1.00
# blockf4                            1.00
# is_surprise_clipNoSurprise:blockf2 1.00
# is_surprise_clipNoSurprise:blockf3 1.00
# is_surprise_clipNoSurprise:blockf4 1.00
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
# sigma    67.46      1.39    64.67    70.14       2000 1.00
# beta    191.96      2.21   187.60   196.27       2000 1.00


marginal_effects(mod1)





final$blockf <- factor(final$block)

mod1 <- lmer(lrt ~ blockf * is_surprise_clip * is_congruent_trial +
               (1 + is_surprise_clip + is_congruent_trial | id) +
               (1 | video_id),
             data = final)

Anova(mod1)


mod2 <- lmer(lrt ~ blockf * is_surprise_clip +
               (1 + is_surprise_clip + blockf | id) +
               (1 | video_id),
             data = final, 
             subset = is_congruent_trial == "Congruent")

Anova(mod2)





# -------------------------------------------------------------------
# rstanarm



standat_df <- correct_df %>%
  dplyr::filter(trials_after_clip != "4" & is_first_trial == "No")

# remove NAs
stan_clean <- standat_df[!is.na(standat_df$vel), ]

mod1 <- stan_lmer(lrt ~ block * is_surprise_clip * is_congruent_trial +
  (1 + block * is_surprise_clip * is_congruent_trial | id) +
  (1 | video_id),
  data = final)

summary(mod1)

posterior_interval(mod1, pars = "block:condition", prob = 0.95)

me <- marginal_effects(mod1, "block:condition")
plot(me, plot = FALSE)[[1]]

launch_shinystan(mod1)

me <- marginal_effects(mod1, "block:condition", re_formula = NULL)
plot(me, plot = FALSE)[[1]]


# -------------------------------------------------------------------

temp <- correct_df[!is.na(correct_df$trt), ]

fit <- brm(vel ~ block * condition * is_surprise_clip * is_clip_prev_1+
  (1 + block * condition * is_surprise_clip * is_clip_prev_1 | id) +
  (1 | video_id),
  data = correct_df,
  chains = 1)

summary(fit)

me <- marginal_effects(fit, "block:condition")
plot(me, plot = FALSE)[[1]]
