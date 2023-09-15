#' Surprise Experiment
#' IIV analysis for both the control and experimental conditions.
#'
#' File: 23_IIV.R
#'
#' This script was last modified on "Thu Jul  5 10:06:59 2018"



#' # R setup

rm(list = ls())
options(show.signif.stars = FALSE)
knitr::opts_chunk$set(warning = FALSE, message = FALSE, tidy = FALSE)

# Load the necessary libraries 
if (!require("pacman"))
  install.packages("pacman")
library("pacman")
pacman::p_load(
  tidyverse, knitr, lme4, car, effects, forcats,rstanarm, ggpubr, rstan, brms, 
  ggeffects, here, retimes, papaja, dplR, psych
)

source("../lib/helpers.R")


#' # Read data 

data_cntr <- read.csv("surprise_cntr.csv")
data_cntr$experiment <- "Control"
data_cntr$id <- as.numeric(factor(data_cntr$id))

data_expr <- read.csv("surprise_expr.csv")
data_expr$experiment <- "Surprise"
data_expr$id <- as.numeric(factor(data_expr$id))

surprise_df <- rbind(data_cntr, data_expr)

# create a unique identifier for each subject
surprise_df$id <- ifelse(surprise_df$experiment == "Surprise",
                         surprise_df$id + length(unique(data_cntr$id)),
                         surprise_df$id)


# Be sure that sequence of trial is in the right order --------------------

surprise_df <- surprise_df %>%
  dplyr::arrange(id, block, date, time_of_day)


# remove wrong trials
surprise_df$remove_trials <- factor(ifelse(
  surprise_df$subject_name == "corrado_caudek_maschio_24_05_1961" &
    surprise_df$date ==  "16-May-2018",
  "yes", "no"
))

surprise_df <- surprise_df[surprise_df$remove_trials == "no", ]





# surprise_df <- read.csv("surprise_data.csv")









# remove NAs on rtTukey

surprise_df <- surprise_df[!is.na(surprise_df$rtTukey), ]

# add number of trials for each subject
surprise_df <- surprise_df %>% 
  group_by(id) %>% 
  mutate(
    bysub_n = length(rtTukey)
  )

hist(surprise_df$bysub_n)


# Detrend by id and block -------------------------------------------------


drt_list <- detrend_id_block(surprise_df)
drt <- do.call("c", drt_list) 
surprise_df$drt <- drt



# Check sequence of trials
if (0) {
  onesub <- surprise_df[surprise_df$id == 92, ]
  # convert time to decimal hours
  onesub$hms <- format(onesub$time_of_day, format = "%H:%M:%S")
  onesub$hms <- as.POSIXct(onesub$hms, format = "%H:%M:%S")
  ggplot(onesub, aes(hms, rtTukey)) + geom_line()
  ggplot(onesub, aes(ntrial, rtTukey)) + geom_line()
  ggplot(onesub, aes(ntrial, drt)) + geom_line()
  mean(onesub$rtTukey)
  mean(onesub$drt)
  var(onesub$rtTukey)
  var(onesub$drt)
}




# Number of trials by subject
kable(
  surprise_df %>% 
    group_by(id) %>% 
    summarise(
      n = length(rtTukey)
    )
)

nrow(surprise_df) / 167
# 338.6587


surprise_df <- surprise_df %>%
  group_by(subject_name) %>%
  mutate(
    tot_trials = n()
  ) %>%
  ungroup()


kable(
  surprise_df %>% 
    group_by(experiment) %>% 
    summarise(
      n = length(unique(subject_name))
    )
)


#' # Are the RTs longer in the surprise experiment than in the control experiment?
#' Comparison between the average latencies of the two experiments.

kable(
  surprise_df %>% 
    dplyr::filter(correct == 1) %>% 
    group_by(experiment, is_surprise_clip) %>% 
    summarise(
      median_rt = median(rtTukey, na.rm = TRUE),
      avg_rt = mean(rtTukey, na.rm = TRUE, trim = 0.1),
      stderr = sqrt(var(rtTukey, na.rm = TRUE) / n())
    )
)


kable(
  surprise_df %>% 
    dplyr::filter(correct == 1 & trials_after_clip < 4) %>% 
    group_by(experiment, is_surprise_clip, trials_after_clip) %>% 
    summarise(
      median_rt = median(rtTukey, na.rm = TRUE),
      avg_rt = mean(rtTukey, na.rm = TRUE, trim = 0.1),
      stderr = sqrt(var(rtTukey, na.rm = TRUE) / n())
    )
)




# iSD ---------------------------------------------------------------------

# Intraindividual standard deviation (iSD) of log(RTs) for each individual's 
# time series.


onesub <- surprise_df %>% 
  dplyr::filter(experiment == "Surprise" & 
                  correct == 1 &
                  block == 1 &
                  #is_congruent_trial == "Incongruent" &
                   id == 50)

plot(onesub$ntrial, onesub$rtTukey, 
     type = 'l',
     ylim = c(300, 900))



bysub_sd <- surprise_df %>% 
  dplyr::filter(correct == 1 & trials_after_clip < 4) %>% 
  group_by(experiment, is_congruent_trial, is_surprise_clip, id, block) %>% 
  summarise(
    iSD = sd(log(drt), na.rm = TRUE)
  )


bysub_sd %>% 
  group_by(experiment, is_congruent_trial) %>% 
  summarise(
    avg_iSD = mean(iSD, na.rm = TRUE),
    sd_iSD = sd(iSD, na.rm = TRUE)
  )



bysub_mean <- surprise_df %>% 
  dplyr::filter(correct == 1 & trials_after_clip < 4) %>% 
  group_by(experiment, is_congruent_trial, id, block) %>% 
  summarise(
    m = mean(log(drt), na.rm = TRUE)
  )


bysub_sd$m <- bysub_mean$m
bysub_sd$iCV <- bysub_sd$iSD / bysub_mean$m


kable(
  bysub_sd %>% 
    group_by(experiment, is_congruent_trial, is_surprise_clip) %>% 
    summarise(
      avg_iSD = mean(iSD, na.rm = TRUE),
      sd_iSD = sd(iSD, na.rm = TRUE)
    )
)

kable(
  bysub_sd %>% 
    group_by(experiment, is_congruent_trial, is_surprise_clip) %>% 
    summarise(
      avg_CV = mean(iCV, na.rm = TRUE),
      sd_CV = sd(iCV, na.rm = TRUE)
    )
)






fit <- brm(iSD ~ experiment * is_congruent_trial +
             (1 + experiment * is_congruent_trial | id),
           data = bysub_sd)

print(fit, digits = 4)


fit1 <- brm(iCV ~ experiment * is_congruent_trial +
             (1 + experiment * is_congruent_trial | id),
           data = bysub_sd)

print(fit1, digits = 4)

marginal_effects(fit1, "is_congruent_trial")
marginal_effects(fit1, "experiment")


fit2 <- brm(iSD ~ experiment * is_congruent_trial + m + 
              (1 + experiment * is_congruent_trial | id),
            data = bysub_sd)

print(fit2, digits = 4)

fm <- lmer(iSD ~ experiment * is_congruent_trial +
             (1 | id),
           data = bysub_sd)
Anova(fm)


library("tseries")

dsub <- surprise_df %>% 
  dplyr::filter(
    id == 2 &
    correct == 1 &
    block == 2
  )

adf.test(dsub$rtTukey)



# Mean square successive difference (MSSD) --------------------------------


diff_rt <- diff(surprise_df$rtTukey, 1)
surprise_df$diff_rt <- c(NA, diff_rt)

surprise_df$diff_rt <- ifelse(surprise_df$trial_within_blk == 1, NA, surprise_df$diff_rt)

bysub_mssd <- surprise_df %>% 
  dplyr::filter(correct == 1 & 
                  trials_after_clip < 4 & 
                  trial_within_blk != 1) %>% 
  group_by(experiment, is_congruent_trial, id, block) %>% 
  summarise(
    mssd = sd(diff_rt, na.rm = TRUE)
  )


bysub_mssd %>% 
  group_by(experiment, is_congruent_trial) %>% 
  summarise(
    avg_mssd = mean(mssd, na.rm = TRUE),
    sd_mssd = sd(mssd, na.rm = TRUE)
  )

fit3 <- brm(mssd ~ experiment * is_congruent_trial + 
              (1 + experiment * is_congruent_trial | id),
            data = bysub_mssd)

print(fit3, digits = 4)

fit4 <- brm(mssd ~ is_congruent_trial + 
              (1 + experiment * is_congruent_trial | id),
            data = bysub_mssd)

LOO(fit3, fit4)


# head(data.frame(surprise_df$rtTukey, surprise_df$diff_rt))



temp <- surprise_df %>% 
  dplyr::select(id, rtTukey) %>% 
  ungroup()

out <- rmssd(temp, id)


onesub <- surprise_df %>% 
  dplyr::filter(id == 1)

df <- data.frame(
  gr = onesub$block, drt = onesub$drt / 1000
)

rmssd(df, gr)

onesub <- onesub %>% 
  mutate(gr = paste(id, block, is_congruent_trial)) 



out <- rmssd(temp, temp$shit)
rmssd(temp[, c(1, 3)])


bysub_df <- surprise_df %>% 
  dplyr::filter(correct == 1) %>% 
  group_by(experiment, id, is_surprise_clip) %>% 
  summarise(
    avg_rt = mean(rtTukey, na.rm = TRUE, trim = 0.1)
  )

fit <- brm(avg_rt ~ experiment * is_surprise_clip +
             (1 | id),
           data = bysub_df)
summary(fit)
marginal_effects(fit, "is_congruent_trial")
marginal_effects(fit, "experiment")




surprise_grp <- surprise_df %>% 
  dplyr::filter(experiment == "Surprise" & correct == 1)

fit <- brm(rtTukey ~ is_surprise_clip +
              (1 + is_surprise_clip | subject_name),
            data = surprise_grp, 
            family = exgaussian())

marginal_effects(fit, "is_surprise_clip")





bform <- bf(delta ~ experiment * q * trial_after + (1 + q + trial_after | subject),
            sigma ~ experiment * q * trial_after + (1 + q + trial_after | subject))

(prior <- get_prior(
  bform,
  data = data_long, 
  family = student())
)

prior$prior[2]  <- "normal(0, 20)" 
prior$prior[10] <- "normal(0, 20)" 


fit2 <- brm(bform,
            data = data_long, 
            prior = prior,
            family = student, 
            warmup = 1000, 
            iter = 2000, 
            chains = 4
)

summary(fit2)






#' # Data subsetting 



# Only correct trials
correct_df <- surprise_df %>%
  dplyr::filter(correct == 1)


# Exclude the last trial of each block (trials_after_clip is equal 
# to 4 only in the last trial of each block)
four_trials_df <- surprise_df %>%
  dplyr::filter(
    correct == 1,
    trials_after_clip != "4"
  )

final <- surprise_df %>%
  dplyr::filter(
    correct == 1 &
      trials_after_clip != "4" &
      block < 5
  ) # there are a few participants with more than 4 block, but
# this complicate the analysis


#' # Accuracy by subject 
#' 
#' Average accuracy for each participant. Red: control; blue = surprise condition.
create_plot_2(surprise_df)


#' # Plot congruency effect by experiment 

#' Boxplot of the congruency effect for each experiment.
create_plot_1(final)


#' # Sequential congruency effect
#' Congruency effect as a function of experiment and congruency in the previous trial.

create_plot_3(final)


#' ## Delta plot 
#' Delta plot for the congruency effect as a function of experiment, by considering all the four 
#' trials after the video clip.
create_plot_4(final, "ALL")
#' Delta plot as a function of experiment by considering only the first trial immediately after
#' the video clip.
create_plot_4(final, "FIRST_TRIAL")
#' Delta plot as a function of experiment by considering only the second trial immediately after
#' the video clip.
create_plot_4(final, "SECOND_TRIAL")
#' Delta plot as a function of experiment by considering only the third trial immediately after
#' the video clip.
create_plot_4(final, "THIRD_TRIAL")
#' Delta plot as a function of experiment by considering only the fourth trial immediately after
#' the video clip.
create_plot_4(final, "FOURTH_TRIAL")



data <- final

quantiles <- c(0.1, 0.3, 0.5, 0.7, 0.9)
exp <- c("Control", "Surprise")
df <- list()
congr_eff_mat <- list()
bytrialsubj_list <- list()

for (i_trial_after in 1:4) {
  
  data <- final[final$trials_after_clip == (i_trial_after - 1), ]
  
  for (i_exp in 1:2) {
    
    data_exp <- data[data$experiment == exp[i_exp], ]
    
    subs_names <- unique(data_exp$subject_name)
    n_subs <- length(subs_names)
    cdf_data <- matrix(0, nrow = length(quantiles), ncol = n_subs)
    
    condition <- c("Incongruent", "Congruent")
    avg_CDF <- list()
    
    for (i_cond in 1:2) {
      data_one_cond <-
        data_exp[data_exp$is_congruent_trial == condition[i_cond], ]
      
      for (i in 1:n_subs) {
        temp_data <- subset(
          data_one_cond,
          data_one_cond$subject_name == subs_names[i]
        )
        cdf_data[, i] <- quantile(
          temp_data$rtTukey,
          quantiles,
          na.rm = TRUE
        )
      }
      avg_CDF[[i_cond]] <- cdf_data
    }
    
    # Compute the congruency effect (inc - con) by subject
    congr_eff_mat[[i_exp]] <- avg_CDF[[1]] - avg_CDF[[2]]
    
  } # end of i_exp loop
  
  
  # number of control subjects
  n_control <- dim(congr_eff_mat[[1]])[2]
  
  # number of subjects in the experimental condition
  n_surprise <- dim(congr_eff_mat[[2]])[2]
  
  bytrialsubj_list[[i_trial_after]] <- rbind(t(congr_eff_mat[[1]]), t(congr_eff_mat[[2]]))
  
} # end of i_trial_after loop


delta_plot_df <- data.frame(
  rbind(
    bytrialsubj_list[[1]],
    bytrialsubj_list[[2]],
    bytrialsubj_list[[3]],
    bytrialsubj_list[[4]]
  )
)

delta_plot_df$trial_after <- rep(
  c(1:4), each = (n_control + n_surprise)
)

delta_plot_df$experiment <- 
  rep(
    c(rep("Control", n_control), rep("Surprise", n_surprise)), 4
  )

delta_plot_df$subject <- rep(c(1:(n_control + n_surprise)), 4)


data_long <- tidyr::gather(delta_plot_df, quantiles, delta, X1:X5, factor_key=TRUE)

data_long %>% 
  group_by(experiment, quantiles) %>% 
  summarise(
    avg_rt = mean(rt, na.rm = TRUE)
  )




data_long$trial_after <- factor(data_long$trial_after)
data_long$q <- as.numeric(data_long$quantiles)


bform <- bf(delta ~ experiment * q * trial_after + (1 + q + trial_after | subject),
            sigma ~ experiment * q * trial_after + (1 + q + trial_after | subject))

(prior <- get_prior(
  bform,
  data = data_long, 
  family = student())
)

prior$prior[2]  <- "normal(0, 20)" 
prior$prior[10] <- "normal(0, 20)" 


fit2 <- brm(bform,
            data = data_long, 
            prior = prior,
            family = student, 
            warmup = 1000, 
            iter = 2000, 
            chains = 4
)

summary(fit2)
# Family: student 
# Links: mu = identity; sigma = log; nu = identity 
# Formula: delta ~ experiment * q * trial_after + (1 + q + trial_after | subject) 
# sigma ~ experiment * q * trial_after + (1 + q + trial_after | subject)
# Data: data_long (Number of observations: 3320) 
# Samples: 1 chains, each with iter = 2000; warmup = 1000; thin = 1;
# total post-warmup samples = 1000
# 
# Group-Level Effects: 
#   ~subject (Number of levels: 166) 
# Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
# sd(Intercept)                                 43.23      2.97    37.53    49.25        431 1.00
# sd(q)                                         13.49      1.01    11.71    15.58        337 1.00
# sd(trial_after2)                              38.16      2.84    32.91    43.60        361 1.00
# sd(trial_after3)                              39.41      3.14    33.46    45.13        288 1.00
# sd(trial_after4)                              48.80      3.60    41.61    55.89        175 1.03
# sd(sigma_Intercept)                            0.45      0.08     0.30     0.60        167 1.00
# sd(sigma_q)                                    0.12      0.02     0.07     0.16         72 1.00
# sd(sigma_trial_after2)                         0.23      0.10     0.02     0.41         72 1.00
# sd(sigma_trial_after3)                         0.18      0.10     0.02     0.39         90 1.00
# sd(sigma_trial_after4)                         0.40      0.07     0.25     0.54        144 1.00
# cor(Intercept,q)                               0.04      0.10    -0.17     0.23        187 1.00
# cor(Intercept,trial_after2)                   -0.38      0.08    -0.52    -0.21        195 1.00
# cor(q,trial_after2)                            0.01      0.10    -0.17     0.21        235 1.01
# cor(Intercept,trial_after3)                   -0.39      0.08    -0.55    -0.22        223 1.00
# cor(q,trial_after3)                           -0.05      0.10    -0.25     0.14        231 1.00
# cor(trial_after2,trial_after3)                 0.53      0.08     0.37     0.67        154 1.00
# cor(Intercept,trial_after4)                   -0.34      0.08    -0.49    -0.16        292 1.00
# cor(q,trial_after4)                           -0.07      0.09    -0.25     0.12        286 1.00
# cor(trial_after2,trial_after4)                 0.44      0.08     0.28     0.58        262 1.00
# cor(trial_after3,trial_after4)                 0.56      0.07     0.42     0.68        228 1.00
# cor(sigma_Intercept,sigma_q)                  -0.43      0.21    -0.72     0.09         76 1.00
# cor(sigma_Intercept,sigma_trial_after2)       -0.28      0.30    -0.74     0.44        265 1.00
# cor(sigma_q,sigma_trial_after2)                0.05      0.32    -0.60     0.67        254 1.00
# cor(sigma_Intercept,sigma_trial_after3)        0.00      0.34    -0.63     0.67        371 1.01
# cor(sigma_q,sigma_trial_after3)               -0.10      0.33    -0.71     0.58        417 1.00
# cor(sigma_trial_after2,sigma_trial_after3)    -0.10      0.41    -0.79     0.70        101 1.00
# cor(sigma_Intercept,sigma_trial_after4)       -0.37      0.20    -0.72     0.05        156 1.03
# cor(sigma_q,sigma_trial_after4)                0.01      0.23    -0.45     0.46        183 1.02
# cor(sigma_trial_after2,sigma_trial_after4)     0.10      0.33    -0.60     0.68         91 1.02
# cor(sigma_trial_after3,sigma_trial_after4)     0.33      0.31    -0.39     0.83         56 1.00
# 
# Population-Level Effects: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
# Intercept                           29.83      6.77    17.09    43.18        207 1.00
# sigma_Intercept                      2.42      0.19     2.03     2.80        268 1.00
# experiment2                         -1.20      7.78   -17.47    13.21        231 1.01
# q                                   -9.45      2.56   -14.23    -4.29        162 1.00
# trial_after2                        14.30      7.19    -0.66    28.29        267 1.00
# trial_after3                        30.99      7.36    16.52    45.99        266 1.00
# trial_after4                        16.51      8.75    -0.81    34.41        319 1.00
# experiment2:q                        1.48      3.01    -4.71     7.44        156 1.00
# experiment2:trial_after2             3.34      8.73   -14.72    20.41        290 1.00
# experiment2:trial_after3           -16.95      8.89   -35.58    -0.80        277 1.00
# experiment2:trial_after4            -2.56     10.19   -23.44    16.35        304 1.00
# q:trial_after2                       3.97      1.74     0.64     7.49        461 1.00
# q:trial_after3                       2.52      1.96    -1.10     6.31        413 1.00
# q:trial_after4                       5.65      2.33     0.99    10.21        488 1.00
# experiment2:q:trial_after2         -11.97      2.29   -16.30    -7.69        512 1.00
# experiment2:q:trial_after3          -7.78      2.45   -12.35    -2.96        488 1.00
# experiment2:q:trial_after4         -10.03      2.69   -15.33    -4.79        534 1.00
# sigma_experiment2                    0.29      0.22    -0.14     0.73        272 1.00
# sigma_q                              0.27      0.05     0.16     0.37        256 1.00
# sigma_trial_after2                  -0.07      0.26    -0.58     0.44        343 1.00
# sigma_trial_after3                  -0.20      0.26    -0.71     0.31        285 1.00
# sigma_trial_after4                   0.17      0.26    -0.31     0.70        293 1.00
# sigma_experiment2:q                 -0.04      0.06    -0.16     0.08        273 1.00
# sigma_experiment2:trial_after2       0.32      0.29    -0.27     0.90        354 1.00
# sigma_experiment2:trial_after3       0.20      0.30    -0.40     0.79        283 1.00
# sigma_experiment2:trial_after4      -0.04      0.29    -0.62     0.51        306 1.00
# sigma_q:trial_after2                -0.04      0.07    -0.19     0.10        342 1.00
# sigma_q:trial_after3                 0.06      0.07    -0.09     0.21        291 1.00
# sigma_q:trial_after4                 0.03      0.07    -0.12     0.16        293 1.00
# sigma_experiment2:q:trial_after2     0.01      0.08    -0.16     0.18        356 1.00
# sigma_experiment2:q:trial_after3    -0.01      0.09    -0.19     0.16        292 1.00
# sigma_experiment2:q:trial_after4    -0.03      0.08    -0.18     0.14        341 1.00
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
# nu    27.02     13.90    10.67    63.82        692 1.00
# 
# Samples were drawn using sampling(NUTS). For each parameter, Eff.Sample 
# is a crude measure of effective sample size, and Rhat is the potential 
# scale reduction factor on split chains (at convergence, Rhat = 1).


plot(fit2, pars = c("experiment2:q:trial_after2", 
                    "experiment2:q:trial_after3", 
                    "experiment2:q:trial_after4")) 



bform <- bf(delta ~ experiment * q + (1 + q + trial_after | subject),
            sigma ~ experiment * q + (1 + q + trial_after | subject))

(prior <- get_prior(
  bform,
  data = data_long, 
  family = student())
)

prior$prior[2]  <- "normal(0, 20)" 


fit3 <- brm(bform,
            data = data_long, 
            prior = prior,
            family = student, 
            warmup = 1000, 
            iter = 2000, 
            chains = 1
)


# # Plot congr. effect by experiment, trial number, and subject -------------
# 
# 
# bysub_rt <- final %>% 
#   group_by(
#     experiment, subject_name, trials_after_clip, is_congruent_trial
#     ) %>% 
#   summarise(
#     mlrt = mean(log(rtTukey), na.rm = TRUE)
#   )
# 
# congr_df <- bysub_rt %>% 
#   dplyr::filter(is_congruent_trial == "Congruent")
# congr_df$is_congruent_trial <- NULL
# congr_df <- rename(congr_df, mlrt_con = mlrt)
# 
# incon_df <- bysub_rt %>% 
#   dplyr::filter(is_congruent_trial == "Incongruent")
# incon_df$is_congruent_trial <- NULL
# incon_df <- rename(incon_df, mlrt_inc = mlrt)
# 
# total <- merge(
#   congr_df, incon_df, 
#   by=c("experiment", "subject_name", "trials_after_clip")
# )
# 
# total$congr_eff <- total$mlrt_inc - total$mlrt_con
# 
# total$trials_after_clip <- factor(total$trials_after_clip)
# 
# total %>%
#   ggplot(aes(x=trials_after_clip, y=congr_eff)) +
#   geom_boxplot(fill="skyblue", notch=FALSE) +
#   geom_jitter(size=0.9, color="orange", width=0.1) +
#   geom_hline(yintercept=0, linetype="dashed", color = "black") +
#   facet_wrap(~ experiment) +
#   xlab("Trial Number After the Video Clip") +
#   ylab("Congruency Effect (log ms)") +
#   theme_apa()
# 
# 
# 
# # Only Surprise experiment ------------------------------------------------
# 
# 
# bysub_rt <- final %>% 
#   dplyr::filter(experiment == "Surprise") %>% 
#   group_by(
#     subject_name, is_surprise_clip, trials_after_clip, 
#     is_congruent_trial
#   ) %>% 
#   summarise(
#     mlrt = mean(log(rtTukey), na.rm = TRUE)
#   )
# 
# congr_df <- bysub_rt %>% 
#   dplyr::filter(is_congruent_trial == "Congruent")
# congr_df$is_congruent_trial <- NULL
# congr_df <- rename(congr_df, mlrt_con = mlrt)
# 
# incon_df <- bysub_rt %>% 
#   dplyr::filter(is_congruent_trial == "Incongruent")
# incon_df$is_congruent_trial <- NULL
# incon_df <- rename(incon_df, mlrt_inc = mlrt)
# 
# total <- merge(
#   congr_df, incon_df, 
#   by=c("subject_name", "is_surprise_clip", "trials_after_clip")
# )
# 
# total$congr_eff <- total$mlrt_inc - total$mlrt_con
# 
# total$trials_after_clip <- factor(total$trials_after_clip)
# 
# total %>%
#   ggplot(aes(x=trials_after_clip, y=congr_eff)) +
#   geom_boxplot(fill="skyblue", notch=TRUE) +
#   geom_jitter(size=0.9, color="orange", width=0.1) +
#   geom_hline(yintercept=0, linetype="dashed", color = "black") +
#   facet_wrap(~ is_surprise_clip) +
#   xlab("Trial Number After the Video Clip") +
#   ylab("Congruency Effect (log ms)") +
#   theme_apa()
# 
# # fine!!!
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# mydata <- mydata[!is.na(mydata$rtTukey), ]
# 
# mydata$rt_bins <- cut(
#   mydata$rtTukey,
#   c(seq(200, 1200, length.out = 10), 2000),
#   include.lowest = TRUE,
#   labels = c("q1", "q2", "q3", "q4", "q5", "q6", "q7", 
#              "q8", "q9", "q10")
# )
# 
# # mydata$rt_bins <- cut(
# #   mydata$rtTukey,
# #   quantile(
# #     mydata$rtTukey,
# #     seq(0, 1, length.out = 11)
# #   ),
# #   include.lowest = TRUE,
# #   labels = c(
# #     "q1", "q2", "q3", "q4", "q5", "q6", "q7",
# #     "q8", "q9", "q10"
# #   )
# # )
# 
# ntot <- nrow(mydata)
# 
# tapply(mydata$correct, mydata$is_congruent_trial, length)
# 
# acc_df <- mydata %>% 
#   group_by(is_congruent_trial, rt_bins) %>% 
#   summarise(
#     acc = mean(correct, na.rm = TRUE),
#     n = n() / ntot * 2
#   )
# 
# 
# ggplot(acc_df, aes(x=rt_bins, y=acc, group=is_congruent_trial)) +
#   geom_line(aes(color=is_congruent_trial)) +
#   geom_point(aes(color=is_congruent_trial))
# 
# 
# 
# 
# # Effect of the experiment on the RTs -------------------------------------
# 
# # Are RTs longer, on average, when the video clips signal a 
# # surprising event?
# 
# hist(final$rtTukey)
# rtdata <- complete.cases(final$rtTukey)
# retimes::mexgauss(rtdata)
# 
# # Parameter of the ex-gaussian distribution in the two conditions
# # control condition
# #         mu      sigma        tau 
# # 0.93562532 0.04419088 0.05892118 
# 
# # surprise condition
# #         mu      sigma        tau 
# # 0.89897785 0.06645159 0.08860212 
# # The difference is in the parameters sigma and tau.
# # Better do this analysis with a Bayesian hierarchic model!
# 
# plot_df <- final %>%
#   group_by(experiment, is_congruent_trial, trials_after_clip) %>%
#   summarise(
#     y = mean(rtTukey, na.rm = TRUE, trim = 0.05),
#     stderr = sqrt(var(rtTukey, na.rm = TRUE) / n())
#   )
# 
# plot_df$lower <- plot_df$y - plot_df$stderr
# plot_df$upper <- plot_df$y + plot_df$stderr
# 
# pd <- position_dodge(0.2)
# 
# ggplot(plot_df, aes(x = trials_after_clip, y = y, color = is_congruent_trial)) +
#   geom_pointrange(aes(ymin = lower, ymax = upper), lwd = 1.05, position = pd) +
#   geom_line(position = pd, lwd = 1.05) +
#   geom_point(position = pd) +
#   ylim(500, 750) +
#   facet_wrap(~ experiment) +
#   xlab("Trials after the video clip") +
#   ylab("Response Latency (ms)") +
#   # ggtitle(title)
#   geom_hline(yintercept = 0, lty = 2) +
#   theme_apa() +
#   theme(legend.position = "bottom") 
# 
# 
# 
# 
# 
# fit1 <- brm(formula = rtTukey ~ is_surprise_trial + (1 + is_surprise_trial | subject_name) ,
#             data = data, family = lognormal(),
#             warmup = 10000, iter = 20000, chains = 4, cores = 8,
#             control = list(adapt_delta = 0.9999))
# 
# 
# 
# 
# 
# 
# # Congruency effect for each subject and condition ------------------------
# 
# bysubraw <- final %>%
#   group_by(
#     subject_name, block, trials_after_clip, 
#     is_congruent_trial
#   ) %>%
#   summarise(
#     mrt = mean(rtTukey, trim = 0.1, na.rm = TRUE)
#   )
# 
# bysub_con <- bysubraw %>% dplyr::filter(
#   is_congruent_trial == "Congruent"
# )
# bysub_inc <- bysubraw %>% dplyr::filter(
#   is_congruent_trial == "Incongruent"
# )
# 
# bysub_df <- merge(
#   bysub_con,
#   bysub_inc,
#   by = c("subject_name", "block", 
#          "trials_after_clip")
# )
# 
# bysub_df <- bysub_df %>% 
#   mutate(
#     interference_eff = mrt.y - mrt.x
#   )
# 
# 
# 
# plot_df <- bysub_df %>%
#     group_by(block, trials_after_clip) %>%
#     summarise(
#       y = mean(interference_eff, trim = .1, na.rm = TRUE),
#       stderr = sqrt(var(interference_eff, na.rm = TRUE) / n()),
#       n()
#     )
# plot_df$lower <- plot_df$y - plot_df$stderr
# plot_df$upper <- plot_df$y + plot_df$stderr
# 
# pd <- position_dodge(0.2)
# 
# ggplot(plot_df, aes(x=trials_after_clip, y=y)) +
#   geom_pointrange(aes(ymin=lower, ymax=upper), lwd=1.05, position=pd) +
#   geom_line(position=pd, lwd=1.05) +
#   geom_point(position=pd, size=5) +
#   facet_wrap(~ block) +
#   xlab("Trials after the video clip") +
#   ylab("Congruency Effect (ms)") +
#   # ggtitle(title)
#   geom_hline(yintercept=0, lty=2) +
#   theme(legend.position="bottom")
# 
# 
# plot_df <- bysub_df %>%
#   group_by(trials_after_clip) %>%
#   summarise(
#     y = mean(interference_eff, trim = .1, na.rm = TRUE),
#     stderr = sqrt(var(interference_eff, na.rm = TRUE) / n()),
#     n()
#   )
# plot_df$lower <- plot_df$y - plot_df$stderr
# plot_df$upper <- plot_df$y + plot_df$stderr
# 
# pd <- position_dodge(0.2)
# 
# ggplot(plot_df, aes(x=trials_after_clip, y=y)) +
#   geom_pointrange(aes(ymin=lower, ymax=upper), lwd=1.05, position=pd) +
#   geom_line(position=pd, lwd=1.05) +
#   geom_point(position=pd, size=5) +
#   xlab("Trials after the video clip") +
#   ylab("Congruency Effect (ms)") +
#   # ggtitle(title)
#   geom_hline(yintercept=0, lty=2) +
#   theme(legend.position="bottom")
# 
# 
# 
# mydf <- bysub_df %>%
#   group_by(subject_name, trials_after_clip) %>%
#   summarise(
#     y = mean(interference_eff, trim = .1, na.rm = TRUE)
#   )
# 
# 
# m1 <- brm(formula = y ~ trials_after_clip + 
#             (1 + trials_after_clip | subject_name),
#           data = mydf, 
#           #family = lognormal(),
#           warmup = 1000, 
#           iter = 2000, 
#           chains = 4)
# 
# summary(m1)
# 
# plot(m1, pars = c("trials_after_clip")) 
# 
# plot(marginal_effects(m1, effects = "trials_after_clip"))
# 
# 
# 
# 
# 
# bysub_df %>%
#   group_by(block, is_clip_trial) %>%
#   summarise(
#     y = mean(interference_eff, na.rm = TRUE)
#   )
# 
# length(unique(final$subject_name))
# 
# 
# dd <- merge(bysub_df, temp, by = "subject_name")
# 
# bysub_df$acc <- temp$acc
# bysub_df$med_rt <- temp$med_rt
# bysub_df
# 
# mean(dd$med_rt) + c(-1, 1) * 2.5 * sd(dd$med_rt)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # Stimulus sequence
# with(
#   df,
#   data.frame(
#     trial_within_blk, is_surprise_clip, is_clip_trial, trials_after_clip
#   )[1:50, ]
# )
# 
# 
# temp <- mydata %>%
#   dplyr::filter(
#     block < 5
#   )
# 
# ggerrorplot(temp,
#             x = "block", 
#             y = "correct",
#             desc_stat = "mean_se",
#             color = "is_congruent_trial", 
#             palette = "jco",
#             facet.by = "is_surprise_clip",
#             ylab = "Accuracy",
#             xlab = "Block",
#             position = position_dodge(0.3) # Adjust the space between bars
# )
# 
# 
# 
# ggerrorplot(temp,
#             x = "block", 
#             y = "rtTukey",
#             desc_stat = "mean_se",
#             color = "is_congruent_trial", 
#             palette = "jco",
#             facet.by = "is_surprise_clip",
#             ylab = "Reaction Times (ms)",
#             xlab = "Block",
#             position = position_dodge(0.3) # Adjust the space between bars
# )
# 
# 
# ggerrorplot(temp,
#             x = "is_clip_trial", 
#             y = "rtTukey",
#             desc_stat = "mean_se",
#             color = "is_congruent_trial", 
#             palette = "jco",
#             facet.by = "is_surprise_clip",
#             ylab = "Reaction Times (ms)",
#             xlab = "Block",
#             position = position_dodge(0.3) # Adjust the space between bars
# )
# 
# 
# 
# 
# 
# # Histogram of accuracy
# out <- mydata %>%
#   group_by(subject_name, is_surprise_clip) %>%
#   summarise(
#     acc = mean(correct, na.rm = TRUE)
#   )
# # data.frame(out)
# 
# p <- ggdensity(out,
#   x = "acc",
#   add = "mean",
#   rug = FALSE,
#   color = "is_surprise_clip",
#   fill = "is_surprise_clip",
#   palette = c("#00AFBB", "#E7B800", "#FC4E07")
# )
# 
# p2 <- ggpar(p,
#   title = "Density Plot for Accuracy",
#   subtitle = "Flanker task after video clip",
#   caption = "Surprise project",
#   xlab = "Accuracy",
#   ylab = "Density",
#   legend.title = "Surprise Video"
# )
# p2
# 
# 
# # Histogram of median RTs as a function of surprise
# out <- mydata %>%
#   dplyr::filter(correct == 1) %>%
#   group_by(subject_name, is_surprise_clip) %>%
#   summarise(
#     median_rt = median(rt)
#   )
# # data.frame(out)
# 
# p <- ggdensity(out,
#   x = "median_rt",
#   add = "mean", rug = FALSE,
#   color = "is_surprise_clip", fill = "is_surprise_clip",
#   palette = c("#00AFBB", "#E7B800", "#FC4E07")
# )
# 
# p2 <- ggpar(p,
#   title = "Density Plot for Response Latencies",
#   subtitle = "Flanker task after video clip",
#   caption = "Surprise project",
#   xlab = "Median Reaction Times (ms)",
#   ylab = "Density",
#   legend.title = "Surprise Video"
# )
# p2
# 
# 
# # Histogram of median RTs as a function of number of trials after the clip
# 
# temp <- mydata
# temp$f_n_after <- factor(temp$trials_after_clip)
# temp <- temp[temp$trials_after_clip != 4, ]
# 
# ggplot() + 
#   geom_density(data=temp, aes(x=lrt, group=f_n_after, fill=f_n_after),alpha=0.5, adjust=2) + 
#   xlab("Log reaction time") +
#   ylab("Density")
# 
# 
# # -------------------------------------------------------------------
# # Log transformation of RTs works best
# par(mfrow = c(1, 3))
# qqPlot(mydata$trt, ylab = "Response Latencies (s)", main = "Reaction Times")
# qqPlot(log(mydata$trt), ylab = "Log Response Latencies (s)", main = "Log Reaction Times")
# qqPlot(1 / mydata$trt, ylab = "1 / RT", main = "Reciprocal of Reaction Times")
# par(mfrow = c(1, 1))
# 
# 
# ggerrorplot(final,
#   x = "trials_after_clip", 
#   y = "lrt",
#   desc_stat = "mean_se",
#   color = "is_congruent_trial", 
#   palette = "jco",
#   # facet.by = "block",
#   ylab = "Log RTs",
#   xlab = "Trial after the Video Clip",
#   position = position_dodge(0.3) # Adjust the space between bars
# )
# 
# 
# 
# ggerrorplot(four_trials_df,
#   x = "trials_after_clip", 
#   y = "lrt",
#   desc_stat = "mean_se",
#   color = "is_congruent_trial", 
#   palette = "jco",
#   facet.by = "is_surprise_clip",
#   ylab = "Log RTs",
#   xlab = "Trial after the Video Clip",
#   position = position_dodge(0.3) # Adjust the space between bars
# )
# 
# 
# b1 <- final %>%
#   dplyr::filter(block == 1)
# 
# ggerrorplot(b1,
#   x = "trials_after_clip",
#   y = "lrt",
#   desc_stat = "mean_se",
#   color = "is_congruent_trial", 
#   palette = "jco",
#   facet.by = "is_surprise_clip",
#   ylab = "Log RTs",
#   xlab = "Trials after the clip",
#   position = position_dodge(0.3) # Adjust the space between bars
# )
# 
# 
# # -------------------------------------------------------------------
# #  table of congruency effect and surprise
# out <- final %>%
#   dplyr::filter(correct == 1) %>%
#   # group_by(block, is_surprise_clip, is_congruent_trial) %>%
#   group_by(is_first_trial, is_congruent_trial) %>%
#   summarise(
#     m = mean(rt, na.rm = TRUE, trim = 0.1),
#     se = sqrt(var(rt, na.rm = TRUE) / n()),
#     n = n()
#   )
# data.frame(out)
# 
# 
# # -------------------------------------------------------------------
# # Plot reaction times
# temp <- final %>%
#   filter(trials_after_clip == 1)
# 
# ggerrorplot(final,
#   x = "block",
#   y = "lrt",
#   desc_stat = "mean_se",
#   color = "is_congruent_trial",
#   palette = "jco",
#   facet.by = "is_surprise_clip",
#   ylab = "Log Reaction Times",
#   xlab = "Block of Trials",
#   position = position_dodge(0.3) # Adjust the space between bars
# )
# 
# 
# ggerrorplot(final,
#             x = "block",
#             y = "lrt",
#             desc_stat = "mean_se",
#             color = "is_surprise_clip",
#             palette = "jco",
#             facet.by = "is_congruent_trial",
#             ylab = "Log Reaction Times",
#             xlab = "Block of Trials",
#             position = position_dodge(0.3) # Adjust the space between bars
# )
# 
# 
# bysub_df <- final %>%
#   group_by(id, block, is_congruent_trial, is_surprise_clip) %>%
#   summarise(
#     mrt = mean(lrt, na.rm = TRUE), 
#     n = n()
#   )
# 
# bysub_df$mrt <- ifelse(bysub_df$n < 5, NA, bysub_df$mrt)
# 
# con_df <- dplyr::filter(bysub_df, is_congruent_trial == "Congruent")
# con_df$is_congruent_trial <- NULL
# inc_df <- dplyr::filter(bysub_df, is_congruent_trial == "Incongruent")
# inc_df$is_congruent_trial <- NULL
# 
# con_df$n <- NULL
# inc_df$n <- NULL
# 
# int_eff_df <- merge(con_df, inc_df, 
#                     by = c("id", "block", "is_surprise_clip"))
# 
# int_eff_df$interf_eff <- int_eff_df$mrt.y - int_eff_df$mrt.x
# 
# 
# ggerrorplot(int_eff_df,
#             x = "block", 
#             y = "interf_eff",
#             desc_stat = "mean_se",
#             color = "is_surprise_clip", 
#             palette = "jco",
#             ylab = "Log Reaction Times",
#             xlab = "Block",
#             position = position_dodge(0.3) # Adjust the space between bars
# )
# 
# 
# 
# 
# 
# 
# b1 <- final %>%
#   filter(trials_after_clip == 1)
# 
# ggerrorplot(b1,
#   x = "block", 
#   y = "lrt",
#   desc_stat = "mean_se",
#   color = "is_congruent_trial", 
#   palette = "jco",
#   facet.by = c("is_surprise_clip"),
#   ylab = "Log Reaction Times",
#   xlab = "Block",
#   position = position_dodge(0.3) # Adjust the space between bars
# )
# 
# 
# ggerrorplot(final,
#   x = "is_surprise_clip", 
#   y = "rt",
#   desc_stat = "mean_se",
#   palette = "jco",
#   # facet.by = "is_surprise_clip",
#   ylab = "Reaction Time (ms)",
#   xlab = "Condition",
#   position = position_dodge(0.3) # Adjust the space between bars
# )
# 
# 
# ggerrorplot(mydata,
#   x = "is_surprise_clip", 
#   y = "correct",
#   desc_stat = "mean_se",
#   palette = "jco",
#   # facet.by = "is_surprise_clip",
#   ylab = "Accuracy",
#   xlab = "Surprise",
#   position = position_dodge(0.3) # Adjust the space between bars
# )
# 
# 
# # -------------------------------------------------------------------
# # Plot accuracy
# all_trials_df <- df
# 
# all_trials_df <- all_trials_df %>%
#   mutate(is_surprise_clip = fct_recode(is_surprise_clip,
#     "No Surprise" = "No",
#     "Surprise" = "Yes"
#   ))
# 
# all_trials_df <- rename(all_trials_df, Condition = is_congruent_trial)
# all_trials_df <- all_trials_df %>%
#   mutate(Condition = fct_recode(Condition,
#     "Incongruent" = "No",
#     "Congruent" = "Yes"
#   ))
# 
# ggerrorplot(all_trials_df,
#   x = "block", y = "correct",
#   desc_stat = "mean_se",
#   color = "Condition", palette = "jco",
#   facet.by = "is_surprise_clip",
#   ylab = "Accuracy",
#   xlab = "Block of Trials",
#   ylim = c(0.75, 1.0),
#   position = position_dodge(0.3) # Adjust the space between bars
# )
# 
# ggerrorplot(all_trials_df,
#   x = "is_surprise_clip", y = "correct",
#   desc_stat = "mean_se",
#   color = "condition", palette = "jco",
#   facet.by = "block",
#   ylab = "Accuracy",
#   xlab = "Video Clip",
#   ylim = c(0.75, 1.0),
#   position = position_dodge(0.3) # Adjust the space between bars
# )
# 
# all_trials_df <- all_trials_df[!is.na(all_trials_df$is_clip_prev_1) &
#   all_trials_df$block < 5, ]
# 
# ggerrorplot(all_trials_df,
#   x = "is_clip_prev_1", y = "correct",
#   desc_stat = "mean_se",
#   color = "condition", palette = "jco",
#   facet.by = c("block"),
#   ylab = "Accuracy",
#   xlab = "First Trial After the Video Clip",
#   ylim = c(0.75, 1.0),
#   position = position_dodge(0.3) # Adjust the space between bars
# )
# 
# # -------------------------------------------------------------------
# 
# names_df <- c(
#   "id", "vel", "block", "i_trial", "video_id", "is_surprise_clip",
#   "is_congruent_trial", "trials_after_clip"
# )
# 
# dd <- dplyr::select(correct_df, names_df)
# 
# block_1_df <- dd %>%
#   dplyr::filter(block == 4)
# 
# subject_means <- block_1_df %>%
#   group_by(id, is_congruent_trial) %>%
#   summarize(vel = mean(vel, na.rm = TRUE))
# subject_means
# 
# # barplot <- ggplot(subject_means, aes(x = is_congruent_trial, y = vel)) +
# #     stat_summary(
# #     geom = "bar",
# #     fun.y = "mean",
# #     col = "black",
# #     fill = "gray70"
# #     )
# # barplot
# 
# subject_means_wide <- spread(subject_means,
#   key = is_congruent_trial,
#   value = vel,
#   sep = "_"
# )
# subject_means_wide
# 
# lims <- c(min(correct_df$vel, na.rm = TRUE), max(correct_df$vel, na.rm = TRUE))
# wsplot <-
#   ggplot(subject_means_wide, aes(
#     x = is_congruent_trial_No,
#     y = is_congruent_trial_Yes
#   )) +
#   geom_point() +
#   geom_abline() +
#   scale_x_continuous("Incongruent", limits = lims) +
#   scale_y_continuous("Congruent", limits = lims) +
#   theme(aspect.ratio = 1)
# wsplot
# 
# # -------------------------------------------------------------------
# 
# message("\nDescriptive step: done!")





out <- surprise_df %>% 
  group_by(id) %>% 
  summarise(
    n = unique(bysub_n)
  )

l <- list()
a <- c(3, 5, 2)
for (i in 1:length(out$n)) {
  l[[i]] = 1:out$n[i]
}

surprise_df$ntrial <- do.call("c", l) 



drt_list <- detrend(surprise_df)
drt <- do.call("c", drt_list) 
surprise_df$drt <- drt



# Check sequence of trials
if (0) {
  onesub <- surprise_df[surprise_df$id == 92, ]
  # convert time to decimal hours
  onesub$hms <- format(onesub$time_of_day, format = "%H:%M:%S")
  onesub$hms <- as.POSIXct(onesub$hms, format = "%H:%M:%S")
  ggplot(onesub, aes(hms, rtTukey)) + geom_line()
  ggplot(onesub, aes(ntrial, rtTukey)) + geom_line()
  ggplot(onesub, aes(ntrial, drt)) + geom_line()
  mean(onesub$rtTukey)
  mean(onesub$drt)
  var(onesub$rtTukey)
  var(onesub$drt)
}




# Number of trials by subject
kable(
  surprise_df %>% 
    group_by(id) %>% 
    summarise(
      n = length(rtTukey)
    )
)

nrow(surprise_df) / 167
# 338.6587


surprise_df <- surprise_df %>%
  group_by(subject_name) %>%
  mutate(
    tot_trials = n()
  ) %>%
  ungroup()


kable(
  surprise_df %>% 
    group_by(experiment) %>% 
    summarise(
      n = length(unique(subject_name))
    )
)


#' # Are the RTs longer in the surprise experiment than in the control experiment?
#' Comparison between the average latencies of the two experiments.

kable(
  surprise_df %>% 
    dplyr::filter(correct == 1) %>% 
    group_by(experiment, is_surprise_clip) %>% 
    summarise(
      median_rt = median(rtTukey, na.rm = TRUE),
      avg_rt = mean(rtTukey, na.rm = TRUE, trim = 0.1),
      stderr = sqrt(var(rtTukey, na.rm = TRUE) / n())
    )
)


kable(
  surprise_df %>% 
    dplyr::filter(correct == 1 & trials_after_clip < 4) %>% 
    group_by(experiment, is_surprise_clip, trials_after_clip) %>% 
    summarise(
      median_rt = median(rtTukey, na.rm = TRUE),
      avg_rt = mean(rtTukey, na.rm = TRUE, trim = 0.1),
      stderr = sqrt(var(rtTukey, na.rm = TRUE) / n())
    )
)




# iSD ---------------------------------------------------------------------

# Intraindividual standard deviation (iSD) of log(RTs) for each individual's 
# time series.


onesub <- surprise_df %>% 
  dplyr::filter(experiment == "Surprise" & 
                  correct == 1 &
                  block == 1 &
                  #is_congruent_trial == "Incongruent" &
                  id == 50)

plot(onesub$ntrial, onesub$rtTukey, 
     type = 'l',
     ylim = c(300, 900))



bysub_sd <- surprise_df %>% 
  dplyr::filter(correct == 1 & trials_after_clip < 4) %>% 
  group_by(experiment, is_congruent_trial, id, block) %>% 
  summarise(
    iSD = sd(log(drt), na.rm = TRUE)
  )


bysub_sd %>% 
  group_by(experiment, is_congruent_trial) %>% 
  summarise(
    avg_iSD = mean(iSD, na.rm = TRUE),
    sd_iSD = sd(iSD, na.rm = TRUE)
  )



bysub_mean <- surprise_df %>% 
  dplyr::filter(correct == 1 & trials_after_clip < 4) %>% 
  group_by(experiment, is_congruent_trial, id, block) %>% 
  summarise(
    m = mean(log(drt), na.rm = TRUE)
  )


bysub_sd$m <- bysub_mean$m
bysub_sd$iCV <- bysub_sd$iSD / bysub_mean$m


kable(
  bysub_sd %>% 
    group_by(experiment, is_congruent_trial) %>% 
    summarise(
      avg_iSD = mean(iSD, na.rm = TRUE),
      sd_iSD = sd(iSD, na.rm = TRUE)
    )
)

kable(
  bysub_sd %>% 
    group_by(experiment, is_congruent_trial) %>% 
    summarise(
      avg_CV = mean(cv, na.rm = TRUE),
      sd_CV = sd(cv, na.rm = TRUE)
    )
)






fit <- brm(iSD ~ experiment * is_congruent_trial +
             (1 + experiment * is_congruent_trial | id),
           data = bysub_sd)

print(fit, digits = 4)


fit1 <- brm(iCV ~ experiment * is_congruent_trial +
              (1 + experiment * is_congruent_trial | id),
            data = bysub_sd)

print(fit1, digits = 4)

marginal_effects(fit1, "is_congruent_trial")
marginal_effects(fit1, "experiment")


fit2 <- brm(iSD ~ experiment * is_congruent_trial + m + 
              (1 + experiment * is_congruent_trial | id),
            data = bysub_sd)

print(fit2, digits = 4)

fm <- lmer(iSD ~ experiment * is_congruent_trial +
             (1 | id),
           data = bysub_sd)
Anova(fm)


library("tseries")

dsub <- surprise_df %>% 
  dplyr::filter(
    id == 2 &
      correct == 1 &
      block == 2
  )

adf.test(dsub$rtTukey)



# Mean square successive difference (MSSD) --------------------------------


diff_rt <- diff(surprise_df$rtTukey, 1)
surprise_df$diff_rt <- c(NA, diff_rt)

surprise_df$diff_rt <- ifelse(surprise_df$trial_within_blk == 1, NA, surprise_df$diff_rt)

bysub_mssd <- surprise_df %>% 
  dplyr::filter(correct == 1 & 
                  trials_after_clip < 4 & 
                  trial_within_blk != 1) %>% 
  group_by(experiment, is_congruent_trial, id, block) %>% 
  summarise(
    mssd = sd(diff_rt, na.rm = TRUE)
  )


bysub_mssd %>% 
  group_by(experiment, is_congruent_trial) %>% 
  summarise(
    avg_mssd = mean(mssd, na.rm = TRUE),
    sd_mssd = sd(mssd, na.rm = TRUE)
  )

fit3 <- brm(mssd ~ experiment * is_congruent_trial + 
              (1 + experiment * is_congruent_trial | id),
            data = bysub_mssd)

print(fit3, digits = 4)

fit4 <- brm(mssd ~ is_congruent_trial + 
              (1 + experiment * is_congruent_trial | id),
            data = bysub_mssd)

LOO(fit3, fit4)


# head(data.frame(surprise_df$rtTukey, surprise_df$diff_rt))



temp <- surprise_df %>% 
  dplyr::select(id, rtTukey) %>% 
  ungroup()

out <- rmssd(temp, id)


onesub <- surprise_df %>% 
  dplyr::filter(id == 1)

df <- data.frame(
  gr = onesub$block, drt = onesub$drt / 1000
)

rmssd(df, gr)

onesub <- onesub %>% 
  mutate(gr = paste(id, block, is_congruent_trial)) 



out <- rmssd(temp, temp$shit)
rmssd(temp[, c(1, 3)])


bysub_df <- surprise_df %>% 
  dplyr::filter(correct == 1) %>% 
  group_by(experiment, id, is_surprise_clip) %>% 
  summarise(
    avg_rt = mean(rtTukey, na.rm = TRUE, trim = 0.1)
  )

fit <- brm(avg_rt ~ experiment * is_surprise_clip +
             (1 | id),
           data = bysub_df)
summary(fit)
marginal_effects(fit, "is_congruent_trial")
marginal_effects(fit, "experiment")




surprise_grp <- surprise_df %>% 
  dplyr::filter(experiment == "Surprise" & correct == 1)

fit <- brm(rtTukey ~ is_surprise_clip +
             (1 + is_surprise_clip | subject_name),
           data = surprise_grp, 
           family = exgaussian())

marginal_effects(fit, "is_surprise_clip")





bform <- bf(delta ~ experiment * q * trial_after + (1 + q + trial_after | subject),
            sigma ~ experiment * q * trial_after + (1 + q + trial_after | subject))

(prior <- get_prior(
  bform,
  data = data_long, 
  family = student())
)

prior$prior[2]  <- "normal(0, 20)" 
prior$prior[10] <- "normal(0, 20)" 


fit2 <- brm(bform,
            data = data_long, 
            prior = prior,
            family = student, 
            warmup = 1000, 
            iter = 2000, 
            chains = 4
)

summary(fit2)




