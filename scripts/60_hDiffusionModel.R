#' ---
#' title: "DDM for flanker task in the surprise project"
#' author: "Corrado Caudek"
#' date: "`r date()`"
#' output:
#'  html_document:
#'    fig_caption: true
#'    toc: true
#'    highlight: tango
#' ---


#' Hierarchical Bayesian Modeling of choice/reaction time data with 
#' the following parameters:
#' "alpha" (boundary separation),
#' "beta" (bias),
#' "delta" (drift rate),
#' "tau" (non-decision time).



# ------------------------------------------------------------------------------
#' @title insure that each subject has at least one error
#'
#' @param input a data.frame and how to proceed
#' @return output a data.frame
#'
manage_errors_by_subj <- function(thedat, 
                                  PROCEDURE) {
  
  # find subject with all correct resonses (choice == 2)
  bysub_p <- thedat %>%
    group_by(subjID) %>%
    summarise(
      avg_p = mean(choice)
    ) %>%
    as.data.frame
  
  index <- bysub_p$avg_p == 2.00
  # list of subjID with perfect performance (no errors)
  bad_id <- bysub_p$subjID[index] 
  
  if (PROCEDURE == 1) { # code as error the first reponse of each subjects 
    # with perfect performance
    # select only subjects with 0 errors
    no_err_df <- thedat[thedat$subjID %in% bad_id, ]
    # find number of trials for those subjects
    no_err_n <- no_err_df %>%
      group_by(subjID) %>%
      summarise(
        n = n()
      ) 
    # add vector ii with trials from 1 to n for each subject
    info <- data.frame(start = 1, len = no_err_n$n)
    no_err_df$ii <- sequence(info$len) + rep(info$start-1, info$len)
    # code as error the first trial of each subject
    index <- no_err_df$ii == 1
    no_err_df$choice[index] <- 1 # 1 = error; 2 = correct 
    # remove ii column 
    no_err_df$ii <- NULL
    # select subjects with at least one error
    with_err_df <- thedat[!thedat$subjID %in% bad_id, ]
    # combine the two dataframes
    ddm_df <- rbind(no_err_df, with_err_df)
  } else { # remove subjects with perfect performance
    ddm_df <- thedat[!thedat$subjID %in% bad_id, ]
  }
  
  ddm_df
}



# Prelims -----------------------------------------------------------------

# rm(list = ls())

Sys.setenv(TZ = "Europe/Rome")

options(
  show.signif.stars = FALSE,
  encoding = "UTF-8"
)

knitr::opts_chunk$set(
  warning = FALSE,
  message = FALSE,
  tidy = FALSE
)

set.seed(998468235L, kind = "L'Ecuyer")
mcopts <- list(preschedule = FALSE, set.seed = TRUE)


if (!require("pacman")) install.packages("pacman")
library("pacman")
pacman::p_load(
  tidyverse, bayesplot, knitr, rmarkdown, forcats, rstan, lme4,
  car, effects, expss, gdata, ggthemes, RColorBrewer, 
  
)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
ggplot2::theme_set(theme_default(base_size = 12))



data_cntr <- read.csv("surprise_cntr.csv")
data_cntr$experiment <- "Control"
data_expr <- read.csv("surprise_expr.csv")
data_expr$experiment <- "Surprise"

surprise_df <- rbind(data_cntr, data_expr)



out <- data_cntr %>% 
  dplyr::filter(trials_after_clip > 0 & trials_after_clip < 4 &
                rtTukey > 200 & rtTukey < 1500 &
                  correct == 1) %>% 
  group_by(id, is_congruent_trial) %>% 
  summarise(
    avg = mean(rtTukey, na.rm = TRUE, trim = 0.1),
    med = median(rtTukey, na.rm = TRUE)
  )

data.frame(out)



#' ## Motivation
#'
#' Consider whether to define as outliers data points at 2.5 or 3.0 standard
#' deviations from the mean, for each participant.  Now the threshold is set
#' to 3.0.


selected_trials <- surprise_df %>% 
  dplyr::filter(
    experiment == "Control" & # "Control" "Surprise"
      is_congruent_trial == "Incongruent" &
      trials_after_clip > 0 & trials_after_clip < 4
  )

# Remove rows with NA on the rt variable
thedat <- selected_trials[!is.na(selected_trials$rtTukey), ]

# subject's id from 1 to n_subjects

thedat$id <- as.numeric(factor(as.character(thedat$id)))
sort(unique(thedat$id))

# select only one subject
temp <- thedat %>% 
  dplyr::filter(id < 3)

# Put data in proper format for choiceRT_ddm()
temp$subjID <- as.numeric(factor(temp$id))
temp$choice <- ifelse(temp$correct == 1, 2, 1)
temp$RT     <- temp$rtTukey / 1000

temp <- as.data.frame(temp)

# choiceRT_ddm() requires subjID to be sorted   
onesubj <- temp %>% 
  arrange(subjID) %>%
  dplyr::select(subjID, choice, RT) %>%
  ungroup()


# sample 15 observations per subject, to reduce the time for the MCMC run
set.seed(1)
res <- do.call(rbind, lapply(split(thedat, thedat$subjID), function(x){x[sample(nrow(x), 15), ]}))

dim(res)
length(unique(thedat$subjID))






# insure that each subject has at least one error response!
PROCEDURE <- 1 
# 1: code as error the first trial of each subjects with
# perfect performance; 0: remove subjects with no errors
# ddm_df <- manage_errors_by_subj(onesubj, PROCEDURE)
ddm_df <- onesubj

# sort again by subjID
ddm_df$subjID <- as.numeric(factor(ddm_df$subjID))
ddm_df <- ddm_df %>%
  dplyr::arrange(subjID)


out <- ddm_df %>% 
  group_by(subjID) %>% 
  summarise(
    n = n(),
    acc = mean(choice)
  )
summary(out)

# data must be provided in an external file, tab delimited
write.fwf(x = ddm_df, file = "ddm_data.txt", sep = "\t")



# Run the drift diffusion model -------------------------------------------

start <- proc.time()
res <- choiceRT_ddm(
  data = "ddm_data.txt",
  nchain = 4,
  ncore = 4
)
end <- proc.time()
end - start


# save workspace
# save.image(file = "partner_ufu_congr.RData")
save.image(file = "congr.RData")
# ---- cache=TRUE, warning=FALSE ----

#' ## Validating a fit in Stan
#' See <https://betanalpha.github.io/assets/case_studies/rstan_workflow.html>

#' ### Checking Split R̂ and Effective Sample Sizes


# Firstly we want to ensure that the split R̂ for each parameter is close to 1.
# Empirically we have found that R̂ > 1.1 is usually indicative of problems in
# the fit.

# Then we want to consider the effective sample size, or n_eff.
# When n_eff / n_transitions < 0.001 the estimators that we use are often biased
# and can significantly overestimate the true effective sample size.

#' ### Checking the Tree Depth
check_treedepth(res)


#' ### Checking the E-BFMI
check_energy(fit)


#' ### Checking Divergences
check_div(res)


# ---- save_results, cache=TRUE, warning=FALSE ----

#' save results with the appropriate name
# output_err_high_congr <- output
congr <- res

# save all the workspace
# save.image("self_fuf_congruent.RData")
# load("self_fuf_congruent.RData")


##  ............................................................................
##  View results                                                            ####
#+  chunk_view_results, cache=TRUE, warning=FALSE

res$fit








# control exp inc
#               mean se_mean    sd   2.5%    25%    50%    75%  97.5% n_eff Rhat
# mu_alpha      1.22    0.00  0.05   1.12   1.18   1.21   1.25   1.31   147 1.01
# mu_beta       0.31    0.00  0.02   0.28   0.30   0.31   0.33   0.35   145 1.01
# mu_delta      3.41    0.01  0.15   3.13   3.31   3.40   3.51   3.71   210 1.02
# mu_tau        0.42    0.00  0.01   0.40   0.42   0.42   0.43   0.43   101 1.01
# control exp con
#              mean se_mean    sd   2.5%    25%    50%    75%  97.5% n_eff Rhat
# mu_alpha      1.12    0.00  0.05   1.02   1.09   1.12   1.16   1.23   157 1.00
# mu_beta       0.44    0.00  0.02   0.40   0.42   0.43   0.45   0.48   300 1.00
# mu_delta      2.64    0.01  0.17   2.29   2.53   2.62   2.75   2.98   300 0.99
# mu_tau        0.36    0.00  0.01   0.34   0.35   0.36   0.36   0.37    59 1.02

# surp exp inc 
#                mean se_mean    sd   2.5%    25%    50%    75%  97.5% n_eff Rhat
# mu_alpha       1.47    0.01  0.04   1.38   1.44   1.47   1.49   1.55    60 1.01
# mu_beta        0.33    0.00  0.01   0.31   0.32   0.33   0.33   0.35   300 1.01
# mu_delta       2.85    0.01  0.10   2.68   2.79   2.85   2.91   3.05   300 1.00
# mu_tau         0.44    0.00  0.01   0.42   0.43   0.44   0.45   0.46    85 1.01
# surp exp con 
#                mean se_mean    sd   2.5%    25%    50%    75%  97.5% n_eff Rhat
# mu_alpha       1.48    0.00  0.05   1.37   1.44   1.48   1.51   1.57   133 1.02
# mu_beta        0.40    0.00  0.01   0.38   0.39   0.40   0.41   0.43   300 1.00
# mu_delta       2.36    0.00  0.08   2.22   2.30   2.35   2.41   2.50   300 1.00
# mu_tau         0.43    0.00  0.01   0.41   0.42   0.43   0.43   0.45   169 1.02






# Congruent ---------------------------------------------------------------

#                mean se_mean    sd    2.5%     25%     50%     75%   97.5% n_eff Rhat
# mu_alpha       1.55    0.00  0.08    1.40    1.50    1.55    1.61    1.72  1217    1
# mu_beta        0.39    0.00  0.02    0.36    0.38    0.39    0.40    0.42  2510    1
# mu_delta       3.47    0.01  0.24    3.03    3.31    3.46    3.62    3.97   877    1
# mu_tau         0.31    0.00  0.01    0.28    0.30    0.31    0.31    0.33  1004    1



# Incongruent -------------------------------------------------------------

#                mean se_mean    sd    2.5%     25%     50%     75%   97.5% n_eff Rhat
# mu_alpha       1.70    0.00  0.08    1.55    1.65    1.70    1.75    1.86  1159 1.00
# mu_beta        0.25    0.00  0.01    0.23    0.25    0.25    0.26    0.27  4000 1.00
# mu_delta       4.31    0.01  0.22    3.90    4.16    4.31    4.46    4.75  1109 1.00
# mu_tau         0.33    0.00  0.01    0.31    0.33    0.33    0.34    0.35   923 1.00



plot(res, type = 'trace', inc_warmup = FALSE)
plot(res)
plotInd(res)
res$allIndPars # congruent
#        alpha      beta    delta       tau subjID
# 1  2.0079758 0.3592145 4.128686 0.2601401      1
# 2  1.4371551 0.3919575 4.043810 0.2936203      2
# 3  1.4458564 0.3593240 4.982959 0.3081721      3
# 4  1.1071491 0.4527925 3.171517 0.3079722      4
# 5  1.4896395 0.4167392 3.406844 0.3154570      5
# 6  1.7671996 0.3225934 4.069348 0.3447223      6
# 7  1.6646303 0.3659685 4.265155 0.3243103      7
# 8  1.4850386 0.4156766 3.290241 0.2707442      8
# 9  1.5812907 0.3445877 5.066955 0.3422681      9
# 10 1.3539047 0.3955503 3.519424 0.3470880     10
# 11 2.7069176 0.3269610 3.538847 0.3302129     11
# 12 1.9870982 0.3847205 2.201355 0.1537250     12
# 13 1.8984548 0.3556841 4.468320 0.3221042     13
# 14 1.3136601 0.4125638 5.000140 0.2587330     14
# 15 0.9618005 0.4909229 2.774066 0.3272157     15
# 16 2.0642708 0.3886287 3.362932 0.3340651     16
# 17 1.7318303 0.4483814 3.093712 0.2741215     17
# 18 1.1924204 0.4254926 1.472097 0.2897705     18
# 19 1.5474193 0.3916401 1.860649 0.2888660     19
# 20 1.4165216 0.4288719 4.518214 0.2921745     20
# 21 2.1349404 0.3020515 2.554375 0.2324362     21
# 22 1.7327870 0.3672472 3.519438 0.3045147     22
# 23 1.3700582 0.3808770 3.061441 0.2832110     23
# 24 1.3005107 0.4025054 4.110249 0.3084359     24
# 25 1.6882720 0.3897076 4.121449 0.2520115     25
# 26 1.2973944 0.4476544 4.571329 0.3008428     26
# 27 1.4188297 0.4237771 3.610311 0.3079139     27

res$allIndPars # incongruent
#       alpha      beta    delta       tau subjID
# 1  1.606265 0.2321620 3.798744 0.3272841      1
# 2  1.440663 0.2636411 4.414022 0.3086217      2
# 3  1.425944 0.2484346 5.141230 0.3314108      3
# 4  1.248049 0.2293029 4.584840 0.3383537      4
# 5  1.672582 0.2650524 3.843497 0.2930422      5
# 6  2.018505 0.2470071 6.123203 0.3528161      6
# 7  1.835222 0.2369453 5.350896 0.3245167      7
# 8  1.850399 0.2369422 4.585816 0.2554665      8
# 9  1.745117 0.2556335 4.722073 0.3463231      9
# 10 1.998313 0.2501407 5.494355 0.3388960     10
# 11 2.340489 0.2549114 4.477693 0.3524307     11
# 12 2.047958 0.2868648 2.310239 0.1994133     12
# 13 2.037935 0.2178612 4.593785 0.3577244     13
# 14 1.066724 0.2629566 3.620897 0.2929443     14
# 15 1.398049 0.2376635 3.879317 0.2885679     15
# 16 2.334142 0.2600661 3.859580 0.3316270     16
# 17 1.870804 0.2357748 4.021689 0.3044071     17
# 18 1.362510 0.2801838 2.611434 0.2851824     18
# 19 1.669221 0.2740659 4.277256 0.3526604     19
# 20 1.629761 0.2542576 5.486544 0.3434289     20
# 21 2.267402 0.2584074 4.388757 0.3425261     21
# 22 2.238720 0.2334746 5.973152 0.3002868     22
# 23 1.269161 0.2884834 3.483017 0.3350857     23
# 24 1.662161 0.2597410 4.562455 0.2903716     24
# 25 1.707809 0.2632716 4.788533 0.2818723     25
# 26 1.440550 0.2444027 4.783592 0.3084092     26
# 27 1.645670 0.2549451 4.358368 0.3150213     27

printFit(res)
rhat(res)
plotInd(res, c("mu_alpha", "mu_beta", "mu_delta", "mu_tau"))




#' ### Original Computing Environment

devtools::session_info("rstan")



out <- df %>% 
  dplyr::filter(correct == 1) %>% 
  group_by(is_congruent_trial, trials_after_clip) %>% 
  summarise(
    p = mean(correct),
    m = mean(rtTukey, na.rm = TRUE), 
    n = n()
  )
data.frame(out)






# Simulated data ----------------------------------------------------------

#' The drift rate is higher in the incongruent than in the congruent condition,
#' also in the control condition. It could be that the video, even if they do
#' present any surprise, still generate some sort of expectation. I have to check
#' with the gabor reward data to see if this is true.

library("rtdists")



rt_sim <- rdiffusion(500, a = 1.5, v = 2.5, t0 = 0.25, z = 0.5)

ddm_df <- data.frame(
  subjID = rep(1:2, each = 250),
  choice = c(1, 1, 1, 1, 1, rep(2, 245), 1, 1, 1, 1, 1, rep(2, 245)),
  RT = rt_sim$rt
)


# data must be provided in an external file, tab delimited
write.fwf(x = ddm_df, file = "ddm_data.txt", sep = "\t")


# Run the drift diffusion model 
start <- proc.time()
res <- choiceRT_ddm(
  data = "ddm_data.txt",
  nchain = 4,
  ncore = 4
)
end <- proc.time()
end - start


# View results 
res$fit





#--------------------------------------------------------------------
# FAMILIAR
## ---- cache=TRUE, warning=FALSE ----

#'  ## Select condition
#'  bla bla

EXP               <- "Familiar"
ID_FACES          <- "UFU"
CONDITION         <- "incongruente"

thedat <- prepareData(self,
                      EXP,
                      ID_FACES,
                      CONDITION)

# mydata <- thedat %>%
#   dplyr::filter(subjID < 10)

thedat %>%
  group_by(subjID) %>%
  summarise(
    avg_p = mean(choice)
  ) %>%
  as.data.frame

bad_id <- c(6, 36)
cleaned_df <- thedat[!thedat$subjID %in% bad_id, ]

cleaned_df$subjID <- as.numeric(factor(cleaned_df$subjID))



# ---- cache=TRUE, warning=FALSE ----

# data must be provided in an external file, tab delimited
write.fwf(x = cleaned_df, file = "data_ddm.txt", sep = "\t")



# ---- cache=TRUE, warning=FALSE ----

ptm <- proc.time()
res1 <- choiceRT_ddm(data    = "data_ddm.txt",
                     niter   = 2000,
                     nwarmup = 1000,
                     nchain  = 4,
                     ncore   = 4)
proc.time() - ptm


# ---- save_results, cache=TRUE, warning=FALSE ----

#' save results with the appropriate name
self_ufu <- res # familiar
fam_ufu  <- res1



##  ............................................................................
##  View results                                                            ####
#+  chunk_view_results, cache=TRUE, warning=FALSE

res$fit

#                mean se_mean    sd    2.5%     25%     50%     75%   97.5% n_eff Rhat
# mu_alpha       1.22    0.00  0.04    1.15    1.20    1.22    1.25    1.30   419 1.01
# mu_beta        0.31    0.00  0.01    0.29    0.30    0.31    0.31    0.32  2833 1.00
# mu_delta       3.77    0.01  0.18    3.41    3.65    3.76    3.89    4.14   469 1.03
# mu_tau         0.30    0.00  0.00    0.29    0.30    0.30    0.30    0.30  1164 1.00


plot(res, type = 'trace', inc_warmup = FALSE)
plot(res)
plotInd(res)
res$allIndPars
printFit(res)
rhat(res)
plotInd(res, c("mu_alpha", "mu_beta", "mu_delta", "mu_tau"))























#-------------------------group comparisons-------------------------
#+ chunk_group_comparisons, cache=TRUE, warning=FALSE

#' After model fitting is complete for both groups, evaluate the group
#' differences on each paramter by examining the posterior distribution
#' of the group mean differences
diff_dist <-
  output_err_high_congr$parVals$mu_tau - output_err_low_congr$parVals$mu_tau
HDIofMCMC(diff_dist) # compute the 95% Highest Density Interval (HDI)
plotHDI(diff_dist)   # plot the group mean differences
# [1] "95% Highest Density Interval (HDI):"
# [1] "Lower bound=-0.0906, Upper bound=-0.0198"

diff_dist <-
  output_err_high_congr$parVals$mu_delta - output_err_low_congr$parVals$mu_delta
HDIofMCMC(diff_dist)
plotHDI(diff_dist)
# [1] "95% Highest Density Interval (HDI):"
# [1] "Lower bound=-0.3228, Upper bound=0.2364"

diff_dist <-
  output_err_high_congr$parVals$mu_alpha - output_err_low_congr$parVals$mu_alpha
HDIofMCMC(diff_dist)
plotHDI(diff_dist)
# [1] "95% Highest Density Interval (HDI):"
# [1] "Lower bound=-0.1817, Upper bound=0.1417"

diff_dist <-
  output_err_high_congr$parVals$mu_beta - output_err_low_congr$parVals$mu_beta
HDIofMCMC(diff_dist)
plotHDI(diff_dist)
# [1] "95% Highest Density Interval (HDI):"
# [1] "Lower bound=-0.0696, Upper bound=0.0248"


thedat <- self %>%
  dplyr::filter(Experiment == "Self" &
                  Identity   == "FUF" &
                  Condition  == "Incongruent")

thedat$response <- ifelse(thedat$correct == 1, "upper", "lower")
thedat$rt <- thedat$rt / 1000

thedat <- thedat[!is.na(thedat$rt), ]

recov <- nlminb(start, ll_diffusion, lower = 0,
                rt=thedat$rt, response=thedat$response)

round(recov$par, 3)



#----------------done!
