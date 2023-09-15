#' EZ Diffusion
#'
#' Calculate the parameters of the EZ diffusion model from mrt, vrt, pc...
#'
#' @param proportion_correct Character. The name of the variable that contains proportion of correct trials.
#' @param rt_variance Character. The name of the variable that contains the variance of reaction times for correct decisions
#' @param rt_mean Character. The name of the variable that contains the mean of reaction times for correct decisions.
#' @param data A \code{data.frame} that conitains the above-specified variables.
#' @param s Numeric. A scaling constant. Take a look at the literature, because different
#' authors use different scaling constants.
#'
#' @examples
#' data <- data.frame(
#'   prop_correct = runif(n = 100, min = 0.2, max = .8)
#'   , rt_mean = runif(n = 100, min = .2, max = .8)
#'   , rt_var = runif(n = 100, min = .01, max = .1)
#' )
#'
#' ezDiffusion(
#'   data = data
#'   , proportion_correct = "prop_correct"
#'   , rt_variance = "rt_var"
#'   , rt_mean = "rt_mean"
#'  )
#'
#' @export


ezDiffusion <- function(proportion_correct, rt_variance, rt_mean, s = 0.1, data){
  
  pc <- data[[proportion_correct]]
  vrt <- data[[rt_variance]]
  mrt <- data[[rt_mean]]
  
  parameters <- data.frame(
    "v" = rep(NA, nrow(data))
    , "a" = rep(NA, nrow(data))
    , "t0" = rep(NA, nrow(data))
  )
  # apply edge corrections if necessary
  if(any(pc==0)){
    pc[pc==0] <- 1e-9
  }
  
  if(any(pc==.5)){
    pc[pc==.5] <- .5 + 1e-9 * sample(rep(-1, 1), size = sum(pc==.5), replace = TRUE)
  }
  
  if(any(pc==1)){
    pc[pc==1] <- 1 - 1e-9
  }
  
  if(any(vrt==0)){
    warning("variances==0 supplied")
  }
  
  # The function qlogis calculates the logit.
  L <- qlogis(pc)
  
  # These give drift rate
  x <- L * (L * pc^2 - L * pc + pc -.5)/vrt
  v <- sign(pc - .5) * s * x^(1/4)
  
  # This gives boundary separation
  a <- s^2 * L / v
  
  # This gives non-decision time
  y <- -v*a/s^2
  mdt <- (a/(2*v)) * (1-exp(y))/(1+exp(y))
  t0 <- mrt - mdt
  
  parameters$v <- v
  parameters$a <- a
  parameters$t0 <- t0
  
  return(parameters)
}

#' The EZ-diffusion model
#'
#' This is a convenience wrapper function for the EZ-diffusion model.
#'
#' @param data A \code{data.frame} that contains the data.
#' @param id Character. The name of the variable that contains the participant identifier.
#' @param rt Character. The name of the variable that contains the reaction times.
#' @param correct Character. The name of the variable that codes correct vs. erroneous decision. Needs to be 0 = erroneous and 1 = correct.
#' @param within Character. Possibly a vector of names of variables that contain within-subjects factors.
#' @param s Numeric. An optional scaling constant.
#'
#' @examples
#' data <- data.frame(
#'   id = rep(1:150, each = 100)
#'   , rt = runif(150 * 100, min = .2, max = .8)
#'   , correct = sample(0:1, size = 150 * 100, replace = TRUE)
#' )
#' ez_diffusion_model(data = data, id = "id", rt = "rt", correct = "correct")

#' @export

ez_diffusion_model <- function(data, id, rt, correct, within = NULL, s = .1){
  correct_trials <- data[data[[correct]]==1, ]
  mrt <- papaja:::fast_aggregate(data = correct_trials, dv = rt, factors = c(id, within), fun = mean)
  vrt <- papaja:::fast_aggregate(data = correct_trials, dv = rt, factors = c(id, within), fun = var)
  pc <- papaja:::fast_aggregate(data = data, dv = correct, factors = c(id, within), fun = mean)
  
  colnames(mrt)[which(colnames(mrt)==rt)] <- "rt_mean"
  colnames(vrt)[which(colnames(vrt)==rt)] <- "rt_variance"
  colnames(pc)[which(colnames(pc)==correct)] <-"proportion_correct"
  
  dat_ <- merge(mrt, vrt)
  dat_ <- merge(dat_, pc)
  
  output <- ezDiffusion(
    data = dat_
    , proportion_correct = "proportion_correct"
    , rt_mean = "rt_mean"
    , rt_variance = "rt_variance"
    , s = s
  )
  cbind(dat_, output)
}

# ez2_moments <- function(data, rt, error, factors){
#
#   means <- papaja:::fast_aggregate(data = data, dv = rt, factors = c(factors, error), fun = mean)
#   variances <- papaja:::fast_aggregate(data = data, dv = rt, factors = c(factors, error), fun = var)
#   proportion_correct <- papaja:::fast_aggregate(data = data, dv = error, factors = c(factors), fun = mean)
#
#   data <- data.frame(
#     "pc" = proportion_correct[[error]]
#     , "mrt1" = means[[rt]][means[[error]]==1]
#     , "mrt0" = means[[rt]][means[[error]]==0]
#     , "vrt1" = variances[[rt]][variances[[error]]==1]
#     , "vrt0" = variances[[rt]][variances[[error]]==0]
#   )
#
#   data
# }

# ez2_moments <- function(data, rt, error, factors, within = NULL){
#
#   colnames(data) <- gsub(colnames(data), pattern = " ", replacement = "__")
#   factors <- gsub(factors, pattern = " ", replacement = "__")
#   within <- gsub(within, pattern = " ", replacement = "__")
#   for(i in factors){
#     data[[i]] <- droplevels(as.factor(data[[i]]))
#   }
#   # within <- gsub(within, pattern = " ", replacement = "__")
#
#   means <- papaja:::fast_aggregate(data = data, dv = rt, factors = c(factors), fun = mean)
#   colnames(means) <- gsub(colnames(means), pattern = rt, replacement = "mrt")
#   variances <- papaja:::fast_aggregate(data = data, dv = rt, factors = c(factors), fun = var)
#   colnames(variances) <- gsub(colnames(variances), pattern = rt, replacement = "vrt")
#   proportion_correct <- papaja:::fast_aggregate(data = data, dv = error, factors = c(factors), fun = mean)
#   colnames(proportion_correct) <- gsub(colnames(proportion_correct), pattern = error, replacement = "pc")
#
#   data <- merge(means, variances)
#   data <- merge(data, proportion_correct)
#
#   data <- data[, c(factors, "pc", "mrt", "vrt")]
#
#
#   if(!is.null(within)){
#     tmp1 <- data[data[[within]]==levels(data[[within]])[1], ]
#     tmp2 <- data[data[[within]]==levels(data[[within]])[2], ]
#     tmp1 <- tmp1[, c(setdiff(factors, within), "mrt", "vrt", "pc")]
#     tmp2 <- tmp2[, c(setdiff(factors, within), "mrt", "vrt", "pc")]
#
#     for(i in c("mrt", "vrt", "pc")){
#       colnames(tmp1)[which(colnames(tmp1)==i)] <- paste0(i, "_", levels(data[[within]])[1])
#       colnames(tmp2)[which(colnames(tmp2)==i)] <- paste0(i, "_", levels(data[[within]])[2])
#     }
#     data <- merge(tmp1, tmp2)
#   }
#
#   colnames(data) <- gsub(colnames(data), pattern = "__", replacement = " ")
#
#   return(data)
# }

# Load packages
library("here")
suppressPackageStartupMessages(library("tidyverse")) 
library("brms")
library("mgcv")
library("hBayesDM")

source(here("libraries", "fnct_surprise.R"))
source(here("libraries", "helpers.R"))

options(mc.cores = parallel::detectCores())

d1 <- read.csv(here("data", "processed", "surprise_control_exps.csv"))

d1 %>% glimpse()
# Number of subjects per esperiment.
d1 %>%
  group_by(experiment) %>% 
  summarise(
    n_distinct(subj_name)
  )

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

# Contrasts.
contrasts(d4$experiment) <- contr.sum(2)
contrasts(d4$experiment)

contr_mat <- t(cbind(
  c(-1, -1),
  c(0, 1),
  c(1, 0))
)
contrasts(d4$is_surprise_clip) <- contr_mat
contrasts(d4$is_surprise_clip)

contrasts(d4$is_congruent_trial) <- contr.sum(2)
contrasts(d4$is_congruent_trial)

# Overall proportion of errors
table(d4$correct)[1] / sum(table(d4$correct))

# model er ~ rt to eyeball where performance becomes better than chance
if (0) {
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
}

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

ds$block <- factor(ds$block)

# remove subject with too few observations
foo <- ds %>% 
  group_by(subj_name) %>% 
  summarise(
    n_block = n_distinct(block)
  )

bad <- foo %>% 
  dplyr::filter(n_block < 2)
bad_subj_names <- bad$subj_name
rm(foo)

foo <- ds[!ds$subj_name %in% bad_subj_names, ]

ds_clean <- foo[foo$subj_name != "ca_ro_1996_12_13_249_f", ]
length(unique(ds_clean$id))

ds_clean$id <- as.numeric(factor(ds_clean$id))


thedat <- remove_eff_block_trial_all_subjs(ds_clean)

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




fm <- lmer(
  v ~ condition +
    (1 | subj_id:condition),
  data = tot
)
summary(fm)



ezDiffusion(
  data = data, 
  proportion_correct = "prop_correct", 
  rt_variance = "rt_var", 
  rt_mean = "rt_mean"
  )

