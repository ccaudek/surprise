# Script name: 201_congruency_effect.R
# Project: surprise
# Script purpose: congruency effect as a function of block, experiment
# @author: Corrado Caudek <corrado.caudek@unifi.it>
# Date Created: Sun Dec 15 11:07:14 2019
# Last Modified Date: Tue Jul 12 10:06:32 2022
# 
# Notes: 


# Load packages ----


library("here")
library("tidyverse")
library("brms")
# library("mgcv")
# library("hBayesDM")
library("emmeans")
library("tidybayes")
library("tidyr")
# library("tidylog", warn.conflicts = FALSE)
library("papaja")


# options(mc.cores = parallel::detectCores())
# theme_set(bayesplot::theme_default(base_family = "sans", base_size = 14))


# Source helper functions ----


source(here("libraries", "fnct_surprise.R"))
source(here("libraries", "helpers.R"))


# Import data ----

df1 <- read.csv(here("data", "processed", "control_2021b.csv"))
df1$experiment <- "control"

df2 <- read.csv(here("data", "processed", "surprise_2021b.csv"))
df2$experiment <- "surprise"

d1 <- rbind(df1, df2)

# d1 <- read.csv(here("data", "processed", "surprise_control_exps.csv"))

d1 %>% glimpse()
# Number of subjects per esperiment.
d1 %>%
  group_by(experiment) %>%
  summarise(
    n_distinct(subj_name)
  )



# Data wrangling ----


# Remove blocks with less than 10 trials
out <- d1 %>%
  group_by(subj_name, block) %>%
  summarise(
    ntrials_per_block = n()
  )
d2 <- left_join(d1, out)
d3 <- d2[d2$ntrials_per_block > 50 & d2$ntrials_per_block < 100, ]

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

# # Contrasts.
# contrasts(d4$experiment) <- contr.sum(2)
# contrasts(d4$experiment)
# 
# contr_mat <- t(cbind(
#   c(-1, -1),
#   c(0, 1),
#   c(1, 0)
# ))
# contrasts(d4$is_surprise_clip) <- contr_mat
# contrasts(d4$is_surprise_clip)
# 
# contrasts(d4$is_congruent_trial) <- contr.sum(2)
# contrasts(d4$is_congruent_trial)

# Overall proportion of errors
table(d4$correct)[1] / sum(table(d4$correct))

# model er ~ rt to eyeball where performance becomes better than chance
if (0) {
  foo <- data.frame(
    er = abs(1 - d4$correct),
    rt = d4$rt
  )

  foo <- foo[complete.cases(foo), ]

  y <- foo$er
  x <- foo$rt

  fit <- mgcv::bam(
    formula = y ~ s(x, bs = "ts", k = 10),
    family = binomial
  )

  preds <- data.frame(
    x = seq(min(x), max(x), length.out = 1e3)
  )

  preds$value <- plogis(predict(fit, newdata = preds))

  plot(preds, type = "l")
  abline(h = 0.5)
  abline(v = 240)
}

d4$RT <- d4$rt / 1000

d4$toss <- d4$RT < 0.240
d4$toss <- ifelse(is.na(d4$RT), TRUE, d4$toss)
d4$toss <- ifelse(is.na(d4$toss), TRUE, d4$toss)
d4$toss <- ifelse(d4$resp == "NORESP", TRUE, d4$toss)
d4$toss <- ifelse(d4$RT > 2.5, TRUE, d4$toss)

d5 <- d4[!d4$toss, ]
ds <- d5[!is.na(d5$subj_name), ]

rm(d1, d2, d3, d4, d5)

ds$error <- ifelse(ds$correct == 1, 0, 1)
ds$id <- as.numeric(factor(ds$subj_name))

ds$block <- factor(ds$block)

# remove subject with too few observations
foo <- ds %>%
  group_by(subj_name) %>%
  summarise(
    n_block = n_distinct(block)
  )

bad <- foo %>%
  dplyr::filter(n_block <= 2)
bad_subj_names <- bad$subj_name
rm(foo)

foo <- ds[!ds$subj_name %in% bad_subj_names, ]

ds_clean <- foo[foo$subj_name != "ca_ro_1996_12_13_249_f", ]
length(unique(ds_clean$id))

ds_clean$id <- as.numeric(factor(ds_clean$id))

# remove effect of block and trial number after the clip
ds_clean <- remove_eff_block_trial_all_subjs(ds_clean)
ds_clean <- ds_clean %>%
  ungroup()

ds_cor <- ds_clean[ds_clean$correct == 1, ]


dd <- ds_cor %>%
  dplyr::select(
    id, experiment, block, is_congruent_trial,
    RT
  )


foo <- dd %>%
  group_by(experiment, id, block, is_congruent_trial) %>%
  summarise(
    mrt = median(RT, na.rm = TRUE)
  )

bysubj <- reshape2::dcast(
  foo,
  id + experiment + block ~ is_congruent_trial,
  value.var = "mrt"
)

bysubj$ce <- 1000 * (bysubj$Incongruent - bysubj$Congruent)

bysubj <- bysubj %>%
  dplyr::select(-c("Congruent", "Incongruent")) %>%
  dplyr::rename(exp = experiment)

bysubj$exp <- relevel(bysubj$exp, ref = "surprise")
contrasts(bysubj$exp)


# Plots ----

plot_df <- bysubj %>%
  group_by(exp, block) %>%
  summarise(
    y = mean(ce, na.rm = TRUE),
    se = sqrt(var(ce) / n())
  )


# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge(0.1) # move them .05 to the left and right

ggplot(
  plot_df,
  aes(x = block, y = y, colour = exp)
) +
  geom_errorbar(aes(ymin = y - se, ymax = y + se), width = .1, position = pd) +
  geom_point(position = pd, size = 4) +
  geom_line(aes(group = exp), position = pd) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  # theme(text=element_text(family="LM Roman 10", size=16)) +
  labs(
    x = "Block",
    y = "Congruency effect (ms)",
    colour = "Experiment"
  ) +
  theme_apa()




# Congruency effect as a function of block and trial ----------------------


dd <- ds_cor %>%
  dplyr::select(
    id, experiment, block, is_first_trial, is_congruent_trial,
    RT
  )


foo <- dd %>%
  group_by(experiment, id, block, is_first_trial, is_congruent_trial) %>%
  summarise(
    mrt = median(RT, na.rm = TRUE)
  )

bysubj <- reshape2::dcast(
  foo,
  id + experiment + block + is_first_trial ~ is_congruent_trial,
  value.var = "mrt"
)

bysubj$ce <- 1000 * (bysubj$Incongruent - bysubj$Congruent)

bysubj <- bysubj %>%
  dplyr::select(-c("Congruent", "Incongruent")) %>%
  dplyr::rename(exp = experiment)

bysubj$exp <- relevel(bysubj$exp, ref = "surprise")
contrasts(bysubj$exp)


# Plots ----

plot_df <- bysubj %>%
  group_by(exp, block, is_first_trial) %>%
  summarise(
    y = mean(ce, na.rm = TRUE),
    se = sqrt(var(ce) / n())
  )


# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge(0.1) # move them .05 to the left and right

ggplot(
  plot_df,
  aes(x = block, y = y, colour = exp)
) +
  geom_errorbar(aes(ymin = y - se, ymax = y + se), width = .1, position = pd) +
  geom_point(position = pd, size = 4) +
  geom_line(aes(group = exp), position = pd) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  # theme(text=element_text(family="LM Roman 10", size=16)) +
  labs(
    x = "Block",
    y = "Congruency effect (ms)",
    colour = "Experiment"
  ) +
  facet_wrap(~ is_first_trial) +
  theme_apa()


# Congruency effect as a function of block and is_surprise_clip ----------------------


dd <- ds_cor %>%
  dplyr::select(
    id, experiment, block, is_surprise_clip, is_congruent_trial,
    RT
  )


foo <- dd %>%
  group_by(experiment, id, block, is_surprise_clip, is_congruent_trial) %>%
  summarise(
    mrt = median(RT, na.rm = TRUE)
  )

bysubj <- reshape2::dcast(
  foo,
  id + experiment + block + is_surprise_clip ~ is_congruent_trial,
  value.var = "mrt"
)

bysubj$ce <- 1000 * (bysubj$Incongruent - bysubj$Congruent)

bysubj <- bysubj %>%
  dplyr::select(-c("Congruent", "Incongruent")) %>%
  dplyr::rename(exp = experiment)

bysubj$exp <- relevel(bysubj$exp, ref = "surprise")
contrasts(bysubj$exp)


# Plots ----

plot_df <- bysubj %>%
  group_by(exp, block, is_surprise_clip) %>%
  summarise(
    y = mean(ce, na.rm = TRUE),
    se = sqrt(var(ce) / n())
  )


# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge(0.1) # move them .05 to the left and right

ggplot(
  plot_df,
  aes(x = block, y = y, colour = exp)
) +
  geom_errorbar(aes(ymin = y - se, ymax = y + se), width = .1, position = pd) +
  geom_point(position = pd, size = 4) +
  geom_line(aes(group = exp), position = pd) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  # theme(text=element_text(family="LM Roman 10", size=16)) +
  labs(
    x = "Block",
    y = "Congruency effect (ms)",
    colour = "Experiment"
  ) +
  facet_wrap(~ is_surprise_clip) +
  theme_apa()



# Congruency effect as a function of block and trials_after_clip ----------------------


dd <- ds_cor %>%
  dplyr::select(
    id, experiment, block, trials_after_clip, is_congruent_trial,
    RT
  )


foo <- dd %>%
  group_by(experiment, id, block, trials_after_clip, is_congruent_trial) %>%
  summarise(
    mrt = median(RT, na.rm = TRUE)
  )

bysubj <- reshape2::dcast(
  foo,
  id + experiment + block + trials_after_clip ~ is_congruent_trial,
  value.var = "mrt"
)

bysubj$ce <- 1000 * (bysubj$Incongruent - bysubj$Congruent)

bysubj <- bysubj %>%
  dplyr::select(-c("Congruent", "Incongruent")) %>%
  dplyr::rename(exp = experiment)

bysubj$exp <- relevel(bysubj$exp, ref = "surprise")
contrasts(bysubj$exp)


# Plots ----

plot_df <- bysubj %>%
  group_by(exp, block, trials_after_clip) %>%
  summarise(
    y = mean(ce, na.rm = TRUE),
    se = sqrt(var(ce) / n())
  )


# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge(0.1) # move them .05 to the left and right

ggplot(
  plot_df,
  aes(x = block, y = y, colour = exp)
) +
  geom_errorbar(aes(ymin = y - se, ymax = y + se), width = .1, position = pd) +
  geom_point(position = pd, size = 4) +
  geom_line(aes(group = exp), position = pd) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  # theme(text=element_text(family="LM Roman 10", size=16)) +
  labs(
    x = "Block",
    y = "Congruency effect (ms)",
    colour = "Experiment"
  ) +
  facet_wrap(~ trials_after_clip) +
  theme_apa()




# Congruency previous trial -----------------------------------------------

# In response-interference tasks, congruency effects are reduced in trials 
# that follow an incongruent trial. This congruence sequence effect (CSE) has 
# been taken to reflect top-down cognitive control processes that monitor for 
# and intervene in case of conflict. 

dd <- ds_cor %>%
  dplyr::select(
    id, experiment, block, trials_after_clip, is_congruent_trial,
    RT
  )

dd <- dd %>%
  group_by(id, block) %>%
  mutate(prev_congr = lag(is_congruent_trial))

foo <- dd %>%
  group_by(experiment, id, block, prev_congr, is_congruent_trial) %>%
  summarise(
    mrt = median(RT, na.rm = TRUE)
  )

bysubj <- reshape2::dcast(
  foo,
  id + experiment + block + prev_congr ~ is_congruent_trial,
  value.var = "mrt"
)

bysubj$ce <- 1000 * (bysubj$Incongruent - bysubj$Congruent)

bysubj <- bysubj %>%
  dplyr::select(-c("Congruent", "Incongruent")) %>%
  dplyr::rename(exp = experiment)

bysubj$exp <- relevel(bysubj$exp, ref = "surprise")
contrasts(bysubj$exp)


# Plots ----

plot_df <- bysubj %>%
  group_by(exp, block, prev_congr) %>%
  summarise(
    y = mean(ce, na.rm = TRUE),
    se = sqrt(var(ce) / n())
  )

plot_df <- na.omit(plot_df) 

# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge(0.1) # move them .05 to the left and right

ggplot(
  plot_df,
  aes(x = block, y = y, colour = exp)
) +
  geom_errorbar(aes(ymin = y - se, ymax = y + se), width = .1, position = pd) +
  geom_point(position = pd, size = 4) +
  geom_line(aes(group = exp), position = pd) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  # theme(text=element_text(family="LM Roman 10", size=16)) +
  labs(
    x = "Block",
    y = "Congruency effect (ms)",
    colour = "Experiment"
  ) +
  facet_wrap(~ prev_congr) +
  theme_apa()



plot_w <- plot_df %>%
  dplyr::select(-se) %>% 
  pivot_wider(names_from = prev_congr, values_from = y) %>% 
  mutate(cse = Incongruent - Congruent) %>% 
  dplyr::select(- c(Incongruent, Congruent))


ggplot(
  plot_w,
  aes(x = block, y = cse, colour = exp)
) +
  #geom_errorbar(aes(ymin = y - se, ymax = y + se), width = .1, position = pd) +
  geom_point(position = pd, size = 4) +
  geom_line(aes(group = exp), position = pd) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  # theme(text=element_text(family="LM Roman 10", size=16)) +
  labs(
    x = "Block",
    y = "Congruency effect (ms)",
    colour = "Experiment"
  ) +
  theme_apa()






























# Violin plot
ggplot(
  data = bysubj, aes(x = block, y = ce)
) +
  facet_wrap(~exp) +
  geom_violin(aes(fill = exp, color = exp)) +
  geom_boxplot(width = .4, outlier.shape = NA) +
  ylim(-300, 250) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "darkgray") +
  labs(
    x = "Block",
    y = "Congruency effect (ms)"
  ) +
  theme(legend.position = "none") 


# Prior Predictive Checks ----

bprior <- c(
  set_prior("student_t(3, 0, 100)", class = "Intercept"),
  set_prior("student_t(3, 0, 100)", class = "b"),
  set_prior("gamma(2, 0.1)", class = "nu"),
  set_prior("lkj(1)", class = "cor"),
  set_prior("student_t(3, 0, 50)", class = "sd")
)

bf <- bf(
  ce ~ exp * block + (block | id)
)

fit <- brm(
  bf,
  prior = bprior,
  data = bysubj,
  family = student(),
  chains = 2,
  cores = parallel::detectCores(),
  iter = 2000,
  sample_prior = TRUE
)

plot(hypothesis(fit, "expcontrol = 0")) 
plot(hypothesis(fit, "block2 = 0"))
plot(hypothesis(fit, "expcontrol:block2 = 0"))


# Model 1 ----

bf <- bf(
  ce ~ exp * block + (block | id)
)

mod1 <- brm(
  bf,
  data = bysubj,
  prior = bprior,
  family = student(),
  init_r = 0.05,
  chains = 4,
  cores = parallel::detectCores(),
  iter = 6000,
  warmup = 1000,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

mod1 <- add_criterion(mod1, "loo")

pp_check(mod1) + xlim(-400, 350)

loo1 <- loo(mod1)
plot(loo1)

print(mod1, 4)

conditional_effects(
  mod1,
  "block:exp"
)

plt_dat <- conditional_effects(mod1, effects = "block:exp")
plt_dat <- plt_dat$`block:exp`


ggplot(aes(x = block, y = `estimate__`, colour = exp), data = plt_dat) +
  geom_point(aes(y = `estimate__`), 
             position = position_dodge(width=0.2), size = 4) + 
  geom_pointrange(aes(ymax = `upper__`, ymin = `lower__`), 
                  position = position_dodge(width=0.2)) + 
  labs(
    x = "Block",
    y = "Congruency effect (ms)",
    colour = "Experiment"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") 



emmeans(mod1, specs = ~ block:exp, type = "ce")
# block exp      emmean lower.HPD upper.HPD
# 1     control   15.35     -4.44     35.49
# 2     control   33.32     12.30     54.70
# 3     control   30.42     10.43     50.18
# 4     control   35.55     17.14     54.73
# 1     surprise -28.64    -40.43    -16.75
# 2     surprise -10.07    -22.21      2.33
# 3     surprise   6.57     -4.72     18.60
# 4     surprise  12.00      1.25     22.27
# 
# Point estimate displayed: median 
# HPD interval probability: 0.95 

# mod1 %>%
#   emmeans("block", by = "exp") %>% 
#   pairs
# surprise_rg <- ref_grid(mod1)
# confint(surprise_rg, by = "block")

post <- posterior_samples(mod1) 

cond_r2 <- bayes_R2(mod1, re_formula = NULL, summary = TRUE)

cond_minus_slopes_r2 <- 
  bayes_R2(
    mod1, 
    re_formula = ~ (1 | subj_name), 
    summary = TRUE)

marginal_r2 <- bayes_R2(mod1, re_formula = NA, summary = TRUE)

cond_r2
cond_minus_slopes_r2
marginal_r2


hyp <- c(
  "Intercept < 0", 
  "block2 > 0",
  "block3 > 0",
  "block4 > 0"
)
hypothesis(mod1, hyp)


hyp <- c(
  "Intercept + expcontrol > 0", 
  "block2 + expcontrol + expcontrol:block2 > 0",
  "block3 + expcontrol + expcontrol:block3 > 0",
  "block4 + expcontrol + expcontrol:block4 > 0"
)
hypothesis(mod1, hyp)




# end --





# Model 2 -----------------------------------------------------------------



bf2 <- bf(
  ce ~ experiment + block + (block | id)
)


mod2 <- brm(
  bf2,
  data = wide,
  # family = shifted_lognormal(),
  init_r = 0.05,
  chains = 4,
  cores = parallel::detectCores(),
  iter = 2000,
  warmup = 1000,
  control = list(adapt_delta = 0.90, max_treedepth = 15)
)
mod2 <- add_criterion(mod2, "waic")

# compare both models
loo_compare(mod1, mod2, criterion = "waic")



bf3 <- bf(
  ce ~ experiment + (1 | id)
)

mod3 <- brm(
  bf3,
  data = wide,
  # family = shifted_lognormal(),
  init_r = 0.05,
  chains = 4,
  cores = parallel::detectCores(),
  iter = 2000,
  warmup = 1000,
  control = list(adapt_delta = 0.90, max_treedepth = 15)
)
mod3 <- add_criterion(mod3, "waic")


# compare both models
loo_compare(mod2, mod3, criterion = "waic")



bf4 <- bf(
  ce ~ block + (block | id)
)

mod4 <- brm(
  bf4,
  data = wide,
  # family = shifted_lognormal(),
  init_r = 0.05,
  chains = 4,
  cores = parallel::detectCores(),
  iter = 2000,
  warmup = 1000,
  control = list(adapt_delta = 0.90, max_treedepth = 15)
)
mod4 <- add_criterion(mod4, "waic")


# compare both models
loo_compare(mod2, mod4, criterion = "waic")



mod1 %>%
  emmeans("exp", by = "block") %>%
  pairs()

samp1 <- mod1 %>%
  emmeans("exp", by = "block") %>%
  pairs() %>%
  gather_emmeans_draws()

samp1 %>%
  median_hdi()


get_hdi <- function(x, level = 95) {
  tmp <- hdrcde::hdr(x, prob = level)
  list(data.frame(mode = tmp$mode[1], lower = tmp$hdr[1, 1], upper = tmp$hdr[1, 2]))
}

samp1 %>%
  summarise(hdi = get_hdi(.value)) %>%
  unnest()



loo1 <- loo(
  mod1,
  cores = parallel::detectCores()
)

plot(loo1)








# Other stuff -------------------------------------------------------------


exp_control <- thedat[thedat$experiment == "control", ]

bysubj_con <- exp_control %>%
  dplyr::filter(is_congruent_trial == "Congruent") %>%
  group_by(id) %>%
  summarise(
    rt_mean = mean(drt, na.rm = TRUE),
    rt_var = var(drt, na.rm = TRUE),
    prop_correct = mean(correct, na.rm = TRUE)
  )

bysubj_inc <- exp_control %>%
  dplyr::filter(is_congruent_trial == "Incongruent") %>%
  group_by(id) %>%
  summarise(
    rt_mean = mean(drt, na.rm = TRUE),
    rt_var = var(drt, na.rm = TRUE),
    prop_correct = mean(correct, na.rm = TRUE)
  )




ez_con <- ezDiffusion(
  data = bysubj_con,
  proportion_correct = "prop_correct",
  rt_variance = "rt_var",
  rt_mean = "rt_mean"
)
ez_con$subj_id <- 1:length(bysubj_con$id)
ez_con$condition <- "congruent"


ez_inc <- ezDiffusion(
  data = bysubj_inc,
  proportion_correct = "prop_correct",
  rt_variance = "rt_var",
  rt_mean = "rt_mean"
)
ez_inc$subj_id <- 1:length(bysubj_inc$id)
ez_inc$condition <- "incongruent"

# v: drift rate
# a: boundary separation
# t0: non decision time


tot <- rbind(ez_con, ez_inc)

hist(tot$v)



bf_v <- bf(
  v ~ condition + (condition | subj_id)
)


mod_v <- brm(
  bf_v,
  data = tot,
  family = shifted_lognormal(),
  init_r = 0.05,
  chains = 4,
  cores = parallel::detectCores(),
  iter = 2000,
  warmup = 1000,
  control = list(adapt_delta = 0.90, max_treedepth = 15)
)

pp_check(mod_v)



print(mod_v, 4)
