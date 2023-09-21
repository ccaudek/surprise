#' ---
#' title: "Analysis of the DSTP parameters"
#' output:
#'   pdf_document:
#'     keep_tex: true
#' ---

library("here")

suppressPackageStartupMessages({
  library("tidyverse")
  library("brms")
  library("cmdstanr")
  library("mice")
  library("tidybayes")
  library("emmeans")
  library("broom.mixed")
  library("patchwork")
})

theme_set(bayesplot::theme_default(base_family = "sans", base_size = 14))
set.seed(123)

params_cntl <- rio::import(
  here::here(
    "data", "processed", "params_data", "dstp", "control_DSTP_params.csv"
  )
)
params_cntl$exp <- "control"

params_surprise <- rio::import(
  here::here(
    "data", "processed", "params_data", "dstp", "surprise_DSTP_params.csv"
  )
)
params_surprise$exp <- "surprise"
# In both experiments, the subjects' id start from 1.
params_surprise$subject <- params_surprise$subject + 500

df <- bind_rows(params_cntl, params_surprise)
df$blk <- factor(df$block)

#' # Parameter A: 
#' Height of the boundary for the response selection diffusion process

m1 <- brm(
  A ~ exp * blk + (1 + blk | subject),
  family = shifted_lognormal(),
  backend = "cmdstanr",
  data = df
)
pp_check(m1)

summary(m1)

conditional_effects(m1, "blk:exp")

#' Parameter C: 
#' Height of the boundary for the stimulus selection diffusion process

hist(df$C)

m2 <- brm(
  C ~ exp * blk + (1 + blk | subject),
  family = shifted_lognormal(),
  backend = "cmdstanr",
  # algorithm = "meanfield",
  data = df
)
pp_check(m2)

summary(m2)

conditional_effects(m2, "blk:exp")

#' Parameter mu_ta: 
#' Drift rate for central target during response selection phase 1

plot(density(df$mu_ta))

m3 <- brm(
  mu_ta ~ exp * blk + (1 + blk | subject),
  family = student(),
  backend = "cmdstanr",
  # algorithm = "meanfield",
  data = df
)
pp_check(m3)

summary(m3)

conditional_effects(m3, "blk:exp")


#' Parameter mu_fl: 
#' Drift rate for the flankers during response selection phase 1

df |> 
  group_by(exp, blk) |> 
  summarize(
    mu_fl = mean(mu_fl, trim = 0.1),
    n = n(),
    stderr = sqrt(var(mu_fl) / n)
  )

# df$mu_fl <- 0.001 + df$mu_fl * 10
plot(density(df$mu_fl))

m4 <- brm(
  mu_fl ~ exp * blk + (1 + blk | subject),
  family = zero_inflated_beta(),
  backend = "cmdstanr",
  # control = list(max_treedepth = 15, adapt_delta = 0.95),
  # algorithm = "meanfield",
  # iter = 10000,
  data = df
)
pp_check(m4)

summary(m4)

conditional_effects(m4, "blk:exp")


#' Parameter mu_ss: 
#' Drift rate for stimulus selection

hist(df$mu_ss)

m5 <- brm(
  mu_ss ~ exp * blk + (1 + blk | subject),
  family = student(),
  backend = "cmdstanr",
  # algorithm = "meanfield",
  data = df
)
pp_check(m5)

summary(m5)

conditional_effects(m5, "blk:exp")

#' Parameter mu_rs2: 
#' Drift rate for phase 2 of response selection

plot(density(df$mu_rs2))

m6 <- brm(
  mu_rs2 ~ exp * blk + (1 + blk | subject),
  family = asym_laplace(),
  backend = "cmdstanr",
  # algorithm = "meanfield",
  data = df
)
pp_check(m6)

summary(m6)

conditional_effects(m6, "blk:exp")


#' Parameter ter: 
#' Drift rate for phase 2 of response selection

plot(density(df$ter))

m7 <- brm(
  ter ~ exp * blk + (1 + blk | subject),
  family = asym_laplace(),
  backend = "cmdstanr",
  # algorithm = "meanfield",
  data = df
)
pp_check(m7)

summary(m7)

conditional_effects(m7, "blk:exp")

