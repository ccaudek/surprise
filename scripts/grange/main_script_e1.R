#==============================================================================
### depression & flanker - experiment 1
# if you have any questions about the scripts, please email 
# grange.jim@gmail.com
#==============================================================================



#==============================================================================
# load required libraries
#==============================================================================
library(tidyverse)
theme_set(theme_light()) # ggplot theme change
library(gridExtra)
library(brms)
library(here)
# devtools::install_github("JimGrange/flankr", ref = "development")
library(flankr)
library(parallel)
library(foreach)
library(doParallel)
library(corrplot)

source("functions.R")
source("model_fit.R")


#==============================================================================
# import the data 
#==============================================================================
setwd(here())

# surprise and control experiments, with both correct and error trials.
data <- get_data()
length(unique(data$subj_id))


#==============================================================================
# perform some participant flanker checks
#==============================================================================


#====================================
# participant flanker accuracy checks
#======================================

# first, tidy the data frame
flanker_data <- tidy_flanker(data)

# get the overall accuracy
flanker_accuracy_overall <- get_flanker_accuracy(flanker_data, overall = TRUE)

# get a list of participants who scored below 80% accuracy
accuracy_removal <- flanker_accuracy_overall %>% 
  filter(accuracy < 0.80) %>% 
  pull(subj_id)

# on prolific, you can "reject" participants who do not perform the task
# properly. So, if some of these poorly-performing participants were low
# accuracy because they were just pressing a single response button, I 
# want to be able to remove these. The below function returns the responses
# for the participants identified in the list of <80% accuracy
# participants. I then manually check each one to see whether they are 
# producing low-effort responses. This is then fedback to prolifics.


# remove the <80% accuracy participants from the flanker data
flanker_data <- flanker_data %>% 
  filter(!subj_id %in% accuracy_removal)

flanker_data |> 
  group_by(experiment) |> 
  summarise(
    n = n_distinct(subj_id)
  )
#     experiment      n
# 1   control         81
# 2   surprise        120

rio::export(
  flanker_data, 
  "flanker_data_clean.csv"
)

#==============================================================================
# prepare flanker data for analysis
#==============================================================================

#======================================
# accuracy
#======================================

# get the accuracy split by congruency
flanker_accuracy <- get_flanker_accuracy(flanker_data, overall = FALSE)

# get the congruency effect
flanker_congruency_accuracy <- get_flanker_congruency_accuracy(flanker_data)

# add the mean accuracy
mean_accuracy <- flanker_congruency_accuracy %>% 
  group_by(subj_id) %>% 
  summarise(mean_accuracy = mean(accuracy))



#======================================
# response time
#======================================

# tidy the data frame: rename the columns, round rts 
flanker_data <- flanker_data %>% 
  rename(id = Participant.Public.ID, rt = Reaction.Time, accuracy = Correct)
flanker_data$rt <- round(flanker_data$rt, 0)

# trim the data >250 < 1500 (as in Hubner & Tobel:
# https://www.frontiersin.org/articles/10.3389/fpsyg.2012.00434/full)
# this removes about 1% of data

flanker_data <- flanker_data %>% 
  filter(rt > 250) %>% 
  filter(rt < 1500)

# how many trials does this leave?
n_trials <- flanker_data %>% 
  group_by(id, congruency) %>% 
  summarise(n = length(rt)) 

# minimum number of trials per condition is 145
# mean is 178
# se is 0.16
mean(n_trials$n)
sd(n_trials$n) / sqrt(length(n_trials$n))
min(n_trials$n)

# get the congruency effect per participant 
# (this function removes error trials)
flanker_congruency_rt <- get_flanker_congruency_rt(flanker_data)

# add to the main data file
all_data <- merge(flanker_congruency_rt, all_data, by = "id")

mean_rt <- flanker_data %>% 
  filter(accuracy == 1) %>% 
  group_by(id) %>% 
  summarise(mean_rt = mean(rt))
all_data <- merge(mean_rt, all_data, by = "id")


flanker_data <- flanker_data %>% 
  mutate(experiment = "experiment_1")
write.csv(flanker_data, "flanker_data_exp_1.csv", row.names = FALSE)

# tidy the demographics data
demographics_data <- tidy_demographics(data$demographics)

# remove the <80% accuracy participants from the demographics data
demographics_data <- demographics_data %>% 
  filter(!Participant.Public.ID %in% accuracy_removal)

# age
age <- demographics_data %>% 
  filter(Question.Key == "dem_dob-year") %>% 
  mutate(Response = as.numeric(as.character(Response))) %>% 
  group_by(Participant.Public.ID) %>% 
  summarise(age = 2019 - Response) %>% 
  summarise(mean_age = mean(age), 
            se_age = sd(age) / sqrt(length(age)))

# gender distribution
gender <- demographics_data %>% 
  filter(Question.Key == "dem_gender") %>% 
  count(Response)

# education
education <- demographics_data %>% 
  filter(Question.Key == "dem_education") %>% 
  mutate(Response = as.numeric(as.character(Response))) %>% 
  summarise(mean = mean(Response), 
            se = sd(Response) / sqrt(length(Response)))

# depression diagnosis
depression_diagnosis <- demographics_data %>% 
  filter(Question.Key == "dem_depression") %>% 
  count(Response)

# depression medication
depression_medication <- demographics_data %>% 
  filter(Question.Key == "dem_medication") %>% 
  count(Response)

# other medication
other_medication <- demographics_data %>% 
  filter(Question.Key == "dem_other_med") %>% 
  count(Response)

# depression_episodes
depression_episodes <- demographics_data %>% 
  filter(Question.Key == "response-2") %>% 
  count(Response)


#==============================================================================
# prepare flanker data for analysis
#==============================================================================


#======================================
# accuracy
#======================================

# get the accuracy split by congruency
flanker_accuracy <- get_flanker_accuracy(flanker_data, overall = FALSE)

# get the congruency effect
flanker_congruency_accuracy <- get_flanker_congruency_accuracy(flanker_data)

# update the main data file
all_data <- merge(flanker_congruency_accuracy, all_data, by = "id")
all_data$accuracy <- round(all_data$accuracy, 3)

# add the mean accuracy
mean_accuracy <- flanker_accuracy %>% 
  group_by(id) %>% 
  summarise(mean_accuracy = mean(accuracy))
all_data <- merge(mean_accuracy, all_data, by = "id")

#======================================
# response time
#======================================

# tidy the data frame: rename the columns, round rts 
flanker_data <- flanker_data %>% 
  rename(id = Participant.Public.ID, rt = Reaction.Time, accuracy = Correct)
flanker_data$rt <- round(flanker_data$rt, 0)

# trim the data >250 < 1500 (as in Hubner & Tobel:
# https://www.frontiersin.org/articles/10.3389/fpsyg.2012.00434/full)
# this removes about 1% of data

flanker_data <- flanker_data %>% 
  filter(rt > 250) %>% 
  filter(rt < 1500)

# how many trials does this leave?
n_trials <- flanker_data %>% 
  group_by(id, congruency) %>% 
  summarise(n = length(rt)) 

# minimum number of trials per condition is 145
# mean is 178
# se is 0.16
mean(n_trials$n)
sd(n_trials$n) / sqrt(length(n_trials$n))
min(n_trials$n)

# get the congruency effect per participant 
# (this function removes error trials)
flanker_congruency_rt <- get_flanker_congruency_rt(flanker_data)

# add to the main data file
all_data <- merge(flanker_congruency_rt, all_data, by = "id")

mean_rt <- flanker_data %>% 
  filter(accuracy == 1) %>% 
  group_by(id) %>% 
  summarise(mean_rt = mean(rt))
all_data <- merge(mean_rt, all_data, by = "id")


flanker_data <- flanker_data %>% 
  mutate(experiment = "experiment_1")
write.csv(flanker_data, "flanker_data_exp_1.csv", row.names = FALSE)




#==============================================================================
# behavioural analysis of the flanker effect
#==============================================================================


#======================================
# response time
#======================================

### model the data

# bayesian model of the response time data
rt_model_data <- flanker_data %>% 
  group_by(id, congruency) %>% 
  summarise(mean_rt = mean(rt))

rt_model <- brm(mean_rt ~ congruency, data = rt_model_data, 
                cores = 4, family = exgaussian())

# Population-Level Effects: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
# Intercept               497.53      4.62   488.91   506.60       2602 1.00
# congruencyincongruent    47.84      5.77    36.67    59.15       3017 1.00
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
# sigma    44.94      3.03    39.18    51.03       3168 1.00
# beta     74.59      4.85    65.48    84.23       2587 1.00


### plot the data
# plot the response time data per congruency level
rt_data <- flanker_data %>% 
  filter(accuracy == 1) %>% 
  group_by(id, congruency) %>% 
  summarise(mean_rt = mean(rt))

rt_congruency_plot <- ggplot(rt_data, aes(x = mean_rt, group = congruency)) + 
  geom_density(aes(colour = congruency,
                   fill = congruency), 
               alpha = 0.8) + 
  scale_colour_manual(values = c("#3182bd", "#9ecae1")) +
  scale_fill_manual(values = c("#3182bd", "#9ecae1")) +
  labs(x = "Response Time (ms)") + 
  theme(legend.position = c(0.75, 0.8)) + 
  ggtitle("A")

rt_flanker_plot <- ggplot(all_data, aes(x = rt)) + 
  geom_density(colour = "black", 
               fill = "#67a9cf") + 
  labs(x = "Flanker Effect (ms)") + 
  ggtitle("B")

setwd(here("paper_figures"))
# pdf("rt_data.pdf", width = 8, height = 4)
# grid.arrange(rt_congruency_plot, rt_flanker_plot, ncol = 2)
# dev.off()
ggsave(
  "rt_data.png",
  grid.arrange(rt_congruency_plot, rt_flanker_plot, ncol = 2),
  width = 8,
  height = 4,
  dpi = 1200
)
setwd(here())



#======================================
# accuracy
#======================================

### model the data

# bayesian model of the response time data
acc_model_data <- flanker_data %>% 
  group_by(id, congruency) %>% 
  summarise(mean_acc = mean(accuracy))

acc_model <- brm(mean_acc ~ congruency, data = acc_model_data, 
                cores = 4, family = skew_normal())

print(summary(acc_model), digits = 4)

# Population-Level Effects: 
#                        Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# Intercept               0.9596    0.0019   0.9559   0.9633       1907 1.0026
# congruencyincongruent  -0.0100    0.0019  -0.0141  -0.0064       2367 1.0017
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# sigma   0.0367    0.0011   0.0347   0.0390       1645 1.0024
# alpha -14.6938    2.1588 -19.3560 -10.8037       1710 1.0031



# plot the data
accuracy_data <- flanker_data %>% 
  group_by(id, congruency) %>% 
  summarise(mean_accuracy = mean(accuracy))

acc_congruency_plot <- ggplot(accuracy_data, aes(x = mean_accuracy, 
                                                 group = congruency)) + 
  geom_density(aes(colour = congruency, 
                   fill = congruency), 
               alpha = 0.8) + 
  scale_colour_manual(values = c("#3182bd", "#9ecae1")) +
  scale_fill_manual(values = c("#3182bd", "#9ecae1")) +
  labs(x = "Accuracy (%)") + 
  theme(legend.position = c(0.25, 0.8)) + 
  ggtitle("C")

accuracy_flanker_plot <- ggplot(all_data, aes(x = accuracy)) + 
  geom_density(colour = "black", 
               fill = "#67a9cf") + 
  labs(x = "Flanker Effect (%)") + 
  ggtitle("D")


setwd(here("paper_figures"))
# pdf("all_behavioural_data.pdf", width = 6, height = 6)
# grid.arrange(rt_congruency_plot, rt_flanker_plot, 
#              acc_congruency_plot, accuracy_flanker_plot, ncol = 2)
# dev.off()
ggsave(
  "all_behavioural_data.png",
  grid.arrange(rt_congruency_plot, rt_flanker_plot,
               acc_congruency_plot, accuracy_flanker_plot, ncol = 2),
  width = 6,
  height = 6,
  dpi = 1200
)
setwd(here())




setwd(here("paper_figures"))
pdf("all_behavioural_data.pdf", width = 6, height = 6)
grid.arrange(rt_congruency_plot, rt_flanker_plot,
             acc_congruency_plot, accuracy_flanker_plot, ncol = 2)
dev.off()
setwd(here())



#==============================================================================
# predicting flanker effects from questionnaire responses
#==============================================================================


#======================================
# QIDS response time
#======================================

# mean response time
qids_mean_rt <- brm(mean_rt ~ qids, data = all_data, cores = 4, 
                    family = exgaussian())

summary(qids_mean_rt)

# Population-Level Effects: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
# Intercept   524.27      8.93   506.95   541.97       3585 1.00
# qids         -0.11      0.66    -1.40     1.20       3732 1.00
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
# sigma    45.12      4.14    37.19    53.37       3512 1.00
# beta     70.86      6.46    58.94    83.92       3032 1.00



# resonse time flanker effect
qids_rt <- brm(rt ~ qids, data = all_data, cores = 4)

summary(qids_rt)


# Population-Level Effects: 
#           Estimate  Est.Error l-95% CI u-95%CI    Eff.Sample Rhat
# Intercept    61.13      3.64    53.77    68.36       4189 1.00
# qids          0.02      0.29    -0.53     0.58       4149 1.00
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
# sigma    29.50      1.20    27.28    31.87       4215 1.00


# get the posterior samples & plot them (but not included in the paper)
qids_rt_samples <-  posterior_samples(qids_rt, pars = c("b_Intercept", 
                                                        "b_qids"))
qids_rt_samples %>% 
  gather(key = parameter, value = value, b_Intercept:b_qids, 
         factor_key = TRUE) %>% 
  ggplot(aes(x = value)) +
  geom_density(fill = "skyblue", colour = "darkblue") +
  facet_wrap(~parameter, scales = "free") +
  theme_bw()

# plot the relationship
qids_rt_plot <- all_data %>% 
  ggplot(aes(x = qids, y = rt)) + 
  geom_abline(intercept = qids_rt_samples$b_Intercept[1:200], 
              slope     = qids_rt_samples$b_qids[1:200],
              size = 1/3, alpha = .3, colour = "#ef8a62") + 
  geom_point(alpha = 0.6, colour = "#67a9cf", size = 2) + 
  labs(x = "QIDS Score", y = "RT Flanker Effect (ms)") + 
  ggtitle("A")


#======================================
# SHPS response time
#======================================

# mean response time
shps_mean_rt <- brm(mean_rt ~ shps, data = all_data, cores = 4, 
                    family = exgaussian())

summary(shps_mean_rt)

# Population-Level Effects: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
# Intercept   552.71     14.69   523.62   580.63       3427 1.00
# shps         -1.07      0.50    -2.04    -0.09       3480 1.00
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
# sigma    43.89      4.08    36.07    51.99       2858 1.00
# beta     71.37      6.37    59.50    84.27       2567 1.00




# resonse time flanker effect
shps_rt <- brm(rt ~ shps, data = all_data, cores = 4)               

summary(shps_rt)

# Population-Level Effects: 
#             Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
# Intercept    63.44      2.44    58.70    68.17       4539 1.00
# shps         -0.65      0.52    -1.67     0.40       4634 1.00
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
# sigma    29.46      1.19    27.26    31.95       3647 1.00

# get the posterior samples & plot them (but not included in the paper)
shps_rt_samples <-  posterior_samples(shps_rt, pars = c("b_Intercept", 
                                                        "b_shps"))


# plot the relationship
shps_rt_plot <- all_data %>% 
  ggplot(aes(x = shps, y = rt)) + 
  geom_abline(intercept = shps_rt_samples$b_Intercept[1:200], 
              slope     = shps_rt_samples$b_shps[1:200],
              size = 1/3, alpha = .3, colour = "#ef8a62") + 
  geom_point(alpha = 0.6, colour = "#67a9cf", size = 2) + 
  labs(x = "SHPS Score", y = "RT Flanker Effect (ms)") + 
  ggtitle("C")


#======================================
# QIDS accuracy
#======================================

# mean accuracy
qids_mean_acc <- brm(mean_accuracy ~ qids, family = skew_normal(), 
                     data = all_data, cores = 4)
print(summary(qids_mean_acc), digits = 4)


# Population-Level Effects: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# Intercept   0.9617    0.0024   0.9569   0.9663       2487 1.0004
# qids       -0.0002    0.0001  -0.0005   0.0000       4089 0.9992
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# sigma   0.0310    0.0013   0.0286   0.0335       1129 1.0034
# alpha -12.4548    2.2476 -17.4954  -8.4460       1469 1.0000


# accuracy flanker effect
qids_acc <- brm(accuracy ~ qids, family = skew_normal(), 
                data = all_data, cores = 4)

print(summary(qids_acc), digits = 4)

# Population-Level Effects: 
#              Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# Intercept  -0.0472    0.0050  -0.0571  -0.0373       3418 0.9998
# qids       -0.0003    0.0004  -0.0010   0.0005       4015 0.9995
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# sigma   0.0465    0.0021   0.0426   0.0505       1723 1.0005
# alpha  -3.8108    0.4795  -4.8030  -2.9187       1716 1.0019




# get the posterior samples & plot them (but not included in the paper)
qids_acc_samples <-  posterior_samples(qids_acc, pars = c("b_Intercept", 
                                                          "b_qids"))

# plot the relationship
qids_acc_plot <- all_data %>% 
  ggplot(aes(x = qids, y = accuracy)) + 
  geom_abline(intercept = qids_acc_samples$b_Intercept[1:200], 
              slope     = qids_acc_samples$b_qids[1:200],
              size = 1/3, alpha = .3, colour = "#ef8a62") + 
  geom_point(alpha = 0.6, colour = "#67a9cf", size = 2) + 
  labs(x = "QIDS Score", y = "Accuracy Flanker Effect (%)") + 
  ggtitle("B")


#======================================
# SHPS accuracy
#======================================

# mean accuracy
shps_mean_acc <- brm(mean_accuracy ~ shps, family = skew_normal(), 
                     data = all_data, cores = 4)               

print(summary(shps_mean_acc), digits = 4)

# Population-Level Effects: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# Intercept   0.9591    0.0043   0.9506   0.9675       3241 1.0009
# shps       -0.0000    0.0001  -0.0003   0.0003       4158 1.0001
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# sigma   0.0310    0.0013   0.0285   0.0336       1163 1.0028
# alpha -12.2550    2.1817 -16.7382  -8.1670       1648 1.0015

# accuracy flanker effect
shps_acc <- brm(accuracy ~ shps, family = skew_normal(), 
                data = all_data, cores = 4)               

print(summary(shps_acc), digits = 4)

# Population-Level Effects: 
#            Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# Intercept  -0.0493    0.0035  -0.0562  -0.0427       2785 1.0007
# shps       -0.0003    0.0006  -0.0015   0.0009       4751 0.9993
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# sigma   0.0464    0.0022   0.0424   0.0510       1756 1.0040
# alpha  -3.7760    0.5070  -4.8364  -2.8586       1715 1.0031

# get the posterior samples & plot them (but not included in the paper)
shps_acc_samples <-  posterior_samples(shps_acc, pars = c("b_Intercept", 
                                                          "b_shps"))

# plot the relationship
shps_acc_plot <- all_data %>% 
  ggplot(aes(x = shps, y = accuracy)) + 
  geom_abline(intercept = shps_acc_samples$b_Intercept[1:200], 
              slope     = shps_acc_samples$b_shps[1:200],
              size = 1/3, alpha = .3, colour = "#ef8a62") + 
  geom_point(alpha = 0.6, colour = "#67a9cf", size = 2) + 
  labs(x = "SHPS Score", y = "RT Flanker Effect (ms)") + 
  ggtitle("D")



setwd(here("paper_figures"))
# pdf("all_behavioural_questionnaire_models.pdf", width = 8, height = 8)
# grid.arrange(qids_rt_plot, qids_acc_plot, 
#              shps_rt_plot, shps_acc_plot, ncol = 2)
# dev.off()
ggsave("all_behavioural_questionnaire_models.png", 
       grid.arrange(qids_rt_plot, qids_acc_plot, 
                    shps_rt_plot, shps_acc_plot, ncol = 2), 
       width = 8, 
       height = 8, 
       dpi = 1200)
setwd(here())


setwd(here("paper_figures"))
pdf("all_behavioural_questionnaire_models.pdf", width = 8, height = 8)
grid.arrange(qids_rt_plot, qids_acc_plot, 
             shps_rt_plot, shps_acc_plot, ncol = 2)
dev.off()
setwd(here())

#==============================================================================
# model fits (this takes a LONG time!)
#==============================================================================

# save the flanker data to the model data directory so we can fit the models
# in a different location (it generates more files, so better for tidyness)
setwd(here("model_data"))

# transform the response times to seconds (which is used by the flankr package)
flanker_data$rt <- flanker_data$rt / 1000
write.csv(flanker_data, "flanker_data.csv", row.names = FALSE)


#======================================
# set some parameters for the model 
# fit routines
#======================================

# during the fit, how many sets of starting parameters should be explored?
n_start_parms <- 50

# what should the variance across starting parameters be?
var_start_parms <- 20

# how many trials to simulate during each iteration of the fit routine whilst
# exploring multiple starting parameters?
n_first_pass <- 1000

# how many trials to simulate during the final fit routine?
n_final_pass <- 50000



# pass the data to the model_fit functions. This is the funcion that fits the
# SSP and DSTP model. Note there are two functions being called:
# 1. fit_[model]_group: fit each model to group-level data. This is used to 
#    assess on a group level which model fits the data better. [TODO]
# 2. fit_[model]_individual: fit each model to individual participants' data.
#    This is used to get parameters per participant, used for the correlations
#    to depression & pleasure scales.


#---
## 1. fit each model to group-level data
# [TODO]
# # fit the dstp model using parallel processing
# fit_dstp_group <- fit_dstp_group_par(flanker_data, n_start_parms,
#                                      var_start_parms, n_first_pass, 
#                                      n_final_pass)
# 
# plotFitDSTP(fit, data)
# 
# plotFitDSTP(first_fit, data)


#----
## 2. fit each model to individual subjects

# fit the dstp model using parallel processing
dstp_individual <- fit_dstp_individual_par(flanker_data, n_start_parms,
                                           var_start_parms, n_first_pass, 
                                           n_final_pass)

# fit the ssp model using parallel processing
ssp_individual <- fit_ssp_individual_par(flanker_data, n_start_parms,
                                         var_start_parms, n_first_pass, 
                                         n_final_pass)

# add dstp to the main data file
all_data_dstp <- read.csv("dstp_fit.csv", header = TRUE)
all_data_dstp <- merge(all_data_dstp, all_data, by = "id")
write.csv(all_data_dstp, "all_data_dstp.csv", row.names = FALSE)

# add ssp to the main data file
all_data_ssp <- read.csv("ssp_fit.csv", header = TRUE)
all_data_ssp <- merge(all_data_ssp, all_data, by = "id")
write.csv(all_data_ssp, "all_data_ssp.csv", row.names = FALSE)

setwd(here())


#======================================
# assess model fit via qq plots
#======================================

# this function does all the work and saves the plot to the plot directory
qq_plot(human_data = flanker_data, dstp_model_data = all_data_dstp, 
        ssp_model_data = all_data_ssp, n_trials = 50000)

#======================================
# compare the DSTP & SSP models
# which fit subject data better?
# (see paper Appendices)
#======================================

bic_values <- tibble(
  dstp = all_data_dstp$bBIC, 
  ssp = all_data_ssp$bBIC
)

delta_bic <- tibble(delta = all_data_dstp$bBIC - all_data_ssp$bBIC) 

# what proportion of participants had better SSP fits than DSTP?
better_ssp <- round(sum(delta_bic$delta > 0) / nrow(delta_bic) * 100, 1)

# what proportion of participants had better DSTP fits than DSTP?
better_dstp <- round(sum(delta_bic$delta <0) / nrow(delta_bic) * 100, 1)

bic_scatter <- ggplot(bic_values, aes(x = dstp, y = ssp)) + 
  geom_point(size = 2.3, alpha = 0.6, colour = '#67a9cf') + 
  geom_abline(intercept = 1,
              slope = 1,
              colour = '#ef8a62', 
              size = 1.2) +
  labs(x = "bBIC DSTP", y = "bBIC SSP") + 
  geom_text(x = 1900, y = 2200, size = 5, label = "DSTP Better") + 
  geom_text(x = 2200, y = 1900, size = 5, label = "SSP Better")

# do a histogram of the difference between SSP and DSTP fit
bic_delta <- ggplot(delta_bic, aes(x = delta)) + 
  geom_histogram(colour = 'black', fill = '#67a9cf', binwidth = 5) + 
  geom_vline(xintercept = 0, size = 2, colour = '#ef8a62') + 
  labs(x = "Delta bBIC", y = "Count") +
  geom_text(x = -50, y = 25, size = 5, 
            label = paste("DSTP Better", 
                          paste("(", better_dstp, "%", ")", sep = ""), 
                          sep = "\n")) + 
  geom_text(x = 50, y = 25, size = 5, 
            label = paste("SSP Better", 
                          paste("(", better_ssp, "%", ")", sep = ""), 
                          sep = "\n"))

# what proportion of participants had better SSP fits than DSTP?
better_ssp <- round(sum(delta_bic$delta > 0) / nrow(delta_bic) * 100, 1)

# draw both plots
setwd(here("paper_figures"))
# pdf("model_bic_comparison.pdf", width = 6, height = 7)
# grid.arrange(bic_scatter, bic_delta, ncol = 1)
# dev.off()

ggsave(
  "model_bic_comparison.png", 
  grid.arrange(bic_scatter, bic_delta, ncol = 1), 
  width = 6, 
  height = 7, 
  dpi = 1200
)
setwd(here())


setwd(here("paper_figures"))
pdf("model_bic_comparison.pdf", width = 6, height = 7)
grid.arrange(bic_scatter, bic_delta, ncol = 1)
dev.off()
setwd(here())


#==============================================================================
# model & questionnaire regressions
#==============================================================================


#======================================
# DSTP QIDS
#======================================

# first, look at the distribution of all 7 parameters
par(mfrow = c(3, 3))
for(i in 1:7){
  plot(density(all_data_dstp[, i + 1]), 
       main = colnames(all_data_dstp[i + 1]))
}


# parameter A
dstp_qids_a <- brm(A ~ qids, family = skew_normal(), 
                   data = all_data_dstp, cores = 4)               

print(summary(dstp_qids_a), digits = 4)

# Population-Level Effects: 
#           Estimate Est.Error l-95% CI u-95% CI     Eff.Sample   Rhat
# Intercept   0.1045    0.0037   0.0972   0.1119       2740 1.0015
# qids       -0.0001    0.0003  -0.0006   0.0004       4192 0.9992
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# sigma   0.0403    0.0018   0.0371   0.0438       1332 1.0055
# alpha   8.0406    1.4649   5.5545  11.2240       1312 1.0027

# get the posterior samples & plot them (but not included in the paper)
dstp_qids_a_samples <-  posterior_samples(dstp_qids_a, pars = c("b_Intercept", 
                                                                "b_qids"))

# plot the relationship
dstp_qids_a_plot <- all_data_dstp %>% 
  ggplot(aes(x = qids, y = A)) + 
  geom_abline(intercept = dstp_qids_a_samples$b_Intercept[1:200], 
              slope     = dstp_qids_a_samples$b_qids[1:200],
              size = 1/3, alpha = .3, colour = "#ef8a62") + 
  geom_point(alpha = 0.6, colour = "#67a9cf", size = 2) + 
  labs(x = "QIDS Score", y = "Parameter  Value") + 
  ggtitle("Parameter: A")



# parameter C
dstp_qids_c <- brm(C ~ qids, family = skew_normal(), 
                   data = all_data_dstp, cores = 4)               

print(summary(dstp_qids_c), digits = 4)

# Population-Level Effects: 
#          Estimate Est.Error l-95% CI u-95% CI    Eff.Sample   Rhat
# Intercept   0.0883    0.0038   0.0808   0.0960       3538 0.9996
# qids        0.0004    0.0003  -0.0002   0.0009       3875 0.9996
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# sigma   0.0337    0.0016   0.0308   0.0369       1258 1.0020
# alpha   3.9279    0.6702   2.8100   5.4010       1450 1.0011

# get the posterior samples & plot them (but not included in the paper)
dstp_qids_c_samples <-  posterior_samples(dstp_qids_c, pars = c("b_Intercept", 
                                                                "b_qids"))

# plot the relationship
dstp_qids_c_plot <- all_data_dstp %>% 
  ggplot(aes(x = qids, y = C)) + 
  geom_abline(intercept = dstp_qids_c_samples$b_Intercept[1:200], 
              slope     = dstp_qids_c_samples$b_qids[1:200],
              size = 1/3, alpha = .3, colour = "#ef8a62") + 
  geom_point(alpha = 0.6, colour = "#67a9cf", size = 2) + 
  labs(x = "QIDS Score", y = "Parameter Value") + 
  ggtitle("Parameter: C")



# parameter mu target
dstp_qids_muT <- brm(muT ~ qids, family = skew_normal(), 
                   data = all_data_dstp, cores = 4)               

print(summary(dstp_qids_muT), digits = 4)

# Population-Level Effects: 
#          Estimate Est.Error l-95% CI u-95% CI    Eff.Sample   Rhat
# Intercept   0.1587    0.0077   0.1435   0.1735       3848 0.9993
# qids       -0.0013    0.0006  -0.0025  -0.0002       4750 0.9991
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# sigma   0.0670    0.0030   0.0615   0.0734       1841 1.0001
# alpha   2.5157    0.5433   1.5276   3.7059       1928 1.0003

# get the posterior samples & plot them (but not included in the paper)
dstp_qids_muT_samples <-  posterior_samples(dstp_qids_muT, 
                                            pars = c("b_Intercept", 
                                                     "b_qids"))

# plot the relationship
dstp_qids_muT_plot <- all_data_dstp %>% 
  ggplot(aes(x = qids, y = muT)) + 
  geom_abline(intercept = dstp_qids_muT_samples$b_Intercept[1:200], 
              slope     = dstp_qids_muT_samples$b_qids[1:200],
              size = 1/3, alpha = .3, colour = "#ef8a62") + 
  geom_point(alpha = 0.6, colour = "#67a9cf", size = 2) + 
  labs(x = "QIDS Score", y = "Parameter Value") + 
  ggtitle("Parameter: mu target")


# parameter mu flanker
dstp_qids_muFl <- brm(muFl ~ qids, family = skew_normal(), 
                     data = all_data_dstp, cores = 4)               

print(summary(dstp_qids_muFl), digits = 4)

# Population-Level Effects: 
#         Estimate Est.Error l-95% CI u-95% CI     Eff.Sample   Rhat
# Intercept   0.1999    0.0109   0.1790   0.2222       3011 1.0003
# qids       -0.0007    0.0008  -0.0022   0.0008       4314 0.9998
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# sigma   0.1082    0.0050   0.0992   0.1189       1677 1.0000
# alpha   4.6073    0.8864   3.0962   6.5762       1558 1.0002

# get the posterior samples & plot them (but not included in the paper)
dstp_qids_muFl_samples <-  posterior_samples(dstp_qids_muFl, 
                                            pars = c("b_Intercept", 
                                                     "b_qids"))

# plot the relationship
dstp_qids_muFl_plot <- all_data_dstp %>% 
  ggplot(aes(x = qids, y = muFl)) + 
  geom_abline(intercept = dstp_qids_muFl_samples$b_Intercept[1:200], 
              slope     = dstp_qids_muFl_samples$b_qids[1:200],
              size = 1/3, alpha = .3, colour = "#ef8a62") + 
  geom_point(alpha = 0.6, colour = "#67a9cf", size = 2) + 
  labs(x = "QIDS Score", y = "Parameter Value") + 
  ggtitle("Parameter: mu flanker")



# parameter mu stimulus selection
dstp_qids_muSS <- brm(muSS ~ qids, family = skew_normal(), 
                      data = all_data_dstp, cores = 4)               

print(summary(dstp_qids_muSS), digits = 4)

# Population-Level Effects: 
#        Estimate Est.Error l-95% CI u-95% CI      Eff.Sample   Rhat
# Intercept   0.5390    0.0159   0.5079   0.5703       4553 0.9998
# qids       -0.0010    0.0012  -0.0034   0.0013       5315 0.9999
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# sigma   0.1375    0.0059   0.1269   0.1497       2429 1.0014
# alpha   1.5358    0.4596   0.5400   2.3174       1407 1.0016

# get the posterior samples & plot them (but not included in the paper)
dstp_qids_muSS_samples <-  posterior_samples(dstp_qids_muSS, 
                                             pars = c("b_Intercept", 
                                                      "b_qids"))

# plot the relationship
dstp_qids_muSS_plot <- all_data_dstp %>% 
  ggplot(aes(x = qids, y = muSS)) + 
  geom_abline(intercept = dstp_qids_muSS_samples$b_Intercept[1:200], 
              slope     = dstp_qids_muSS_samples$b_qids[1:200],
              size = 1/3, alpha = .3, colour = "#ef8a62") + 
  geom_point(alpha = 0.6, colour = "#67a9cf", size = 2) + 
  labs(x = "QIDS Score", y = "Parameter Value") + 
  ggtitle("Parameter: mu SS")



# parameter mu response selection, stage 2
dstp_qids_muRS2 <- brm(muRS2 ~ qids, family = skew_normal(), 
                      data = all_data_dstp, cores = 4)               

print(summary(dstp_qids_muRS2), digits = 4)

# Population-Level Effects: 
#         Estimate Est.Error l-95% CI u-95% CI       Eff.Sample   Rhat
# Intercept   1.2398    0.0309   1.1781   1.3010       4637 0.9994
# qids        0.0008    0.0024  -0.0039   0.0054       5813 0.9992
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# sigma   0.2544    0.0111   0.2340   0.2778       2738 1.0000
# alpha   1.5563    0.7318  -0.2938   2.7914       1630 1.0009

# get the posterior samples & plot them (but not included in the paper)
dstp_qids_muRS2_samples <-  posterior_samples(dstp_qids_muRS2, 
                                             pars = c("b_Intercept", 
                                                      "b_qids"))

# plot the relationship
dstp_qids_muRS2_plot <- all_data_dstp %>% 
  ggplot(aes(x = qids, y = muRS2)) + 
  geom_abline(intercept = dstp_qids_muRS2_samples$b_Intercept[1:200], 
              slope     = dstp_qids_muRS2_samples$b_qids[1:200],
              size = 1/3, alpha = .3, colour = "#ef8a62") + 
  geom_point(alpha = 0.6, colour = "#67a9cf", size = 2) + 
  labs(x = "QIDS Score", y = "Parameter Value") + 
  ggtitle("Parameter: mu RS2")


# parameter non-decision time
dstp_qids_ter<- brm(ter ~ qids, family = skew_normal(), 
                       data = all_data_dstp, cores = 4)               

print(summary(dstp_qids_ter), digits = 4)

# Population-Level Effects: 
#        Estimate Est.Error l-95% CI u-95% CI      Eff.Sample   Rhat
# Intercept   0.2881    0.0056   0.2771   0.2993       3637 0.9999
# qids        0.0000    0.0004  -0.0008   0.0009       3805 1.0001
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# sigma   0.0449    0.0019   0.0416   0.0490       1548 1.0012
# alpha  -0.3652    0.7087  -1.5221   0.9737       1283 1.0003

# get the posterior samples & plot them (but not included in the paper)
dstp_qids_ter_samples <-  posterior_samples(dstp_qids_ter, 
                                              pars = c("b_Intercept", 
                                                       "b_qids"))

# plot the relationship
dstp_qids_ter_plot <- all_data_dstp %>% 
  ggplot(aes(x = qids, y = ter)) + 
  geom_abline(intercept = dstp_qids_ter_samples$b_Intercept[1:200], 
              slope     = dstp_qids_ter_samples$b_qids[1:200],
              size = 1/3, alpha = .3, colour = "#ef8a62") + 
  geom_point(alpha = 0.6, colour = "#67a9cf", size = 2) + 
  labs(x = "QIDS Score", y = "Parameter Value") + 
  ggtitle("Parameter: ter")


# draw all plots
setwd(here("paper_figures"))
ggsave("dstp_qids.png", 
       grid.arrange(dstp_qids_a_plot, dstp_qids_c_plot, dstp_qids_muT_plot, 
                    dstp_qids_muFl_plot, dstp_qids_muSS_plot, dstp_qids_muRS2_plot, 
                    dstp_qids_ter_plot, ncol = 3), 
       width = 10, 
       height = 10, 
       dpi = 1200)
setwd(here())


setwd(here("paper_figures"))
pdf("dstp_qids.pdf", width = 10, height = 10)
grid.arrange(dstp_qids_a_plot, dstp_qids_c_plot, dstp_qids_muT_plot, 
             dstp_qids_muFl_plot, dstp_qids_muSS_plot, dstp_qids_muRS2_plot, 
             dstp_qids_ter_plot, ncol = 3)
dev.off()
setwd(here())



#======================================
# DSTP SHPS
#======================================


# parameter A
dstp_shps_a <- brm(A ~ shps, family = skew_normal(), 
                   data = all_data_dstp, cores = 4)               

print(summary(dstp_shps_a), digits = 4)

# Population-Level Effects: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# Intercept   0.1056    0.0058   0.0940   0.1168       2948 1.0005
# shps       -0.0001    0.0002  -0.0005   0.0003       3340 0.9998
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# sigma   0.0402    0.0017   0.0370   0.0439       1379 1.0018
# alpha   8.0047    1.5170   5.4833  11.3615       1006 1.0031

# get the posterior samples & plot them (but not included in the paper)
dstp_shps_a_samples <-  posterior_samples(dstp_shps_a, pars = c("b_Intercept", 
                                                                 "b_shps"))

# plot the relationship
dstp_shps_a_plot <- all_data_dstp %>% 
  ggplot(aes(x = shps, y = A)) + 
  geom_abline(intercept = dstp_shps_a_samples$b_Intercept[1:200], 
              slope     = dstp_shps_a_samples$b_shps[1:200],
              size = 1/3, alpha = .3, colour = "#ef8a62") + 
  geom_point(alpha = 0.6, colour = "#67a9cf", size = 2) + 
  labs(x = "SHPS Score", y = "Parameter  Value") + 
  ggtitle("Parameter: A")



# parameter C
dstp_shps_c <- brm(C ~ shps, family = skew_normal(), 
                   data = all_data_dstp, cores = 4)               

print(summary(dstp_shps_c), digits = 4)

# Population-Level Effects: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# Intercept   0.0872    0.0066   0.0743   0.0998       3412 1.0010
# shps        0.0002    0.0002  -0.0002   0.0006       3587 1.0012
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# sigma   0.0337    0.0016   0.0308   0.0369       1277 1.0019
# alpha   3.8494    0.6403   2.7775   5.2774       1384 1.0021

# get the posterior samples & plot them (but not included in the paper)
dstp_shps_c_samples <-  posterior_samples(dstp_shps_c, pars = c("b_Intercept", 
                                                                "b_shps"))

# plot the relationship
dstp_shps_c_plot <- all_data_dstp %>% 
  ggplot(aes(x = shps, y = C)) + 
  geom_abline(intercept = dstp_shps_c_samples$b_Intercept[1:200], 
              slope     = dstp_shps_c_samples$b_shps[1:200],
              size = 1/3, alpha = .3, colour = "#ef8a62") + 
  geom_point(alpha = 0.6, colour = "#67a9cf", size = 2) + 
  labs(x = "SHPS Score", y = "Parameter  Value") + 
  ggtitle("Parameter: C")



# parameter muT
dstp_shps_muT <- brm(muT ~ shps, family = skew_normal(), 
                   data = all_data_dstp, cores = 4)               

print(summary(dstp_shps_muT), digits = 4)

# 
# Population-Level Effects: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# Intercept   0.1547    0.0134   0.1281   0.1806       4237 0.9995
# shps       -0.0004    0.0005  -0.0013   0.0005       4510 0.9993
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# sigma   0.0677    0.0031   0.0619   0.0741       1546 1.0005
# alpha   2.5639    0.5348   1.6016   3.6820       1543 1.0010

# get the posterior samples & plot them (but not included in the paper)
dstp_shps_muT_samples <-  posterior_samples(dstp_shps_muT, 
                                            pars = c("b_Intercept", 
                                                     "b_shps"))

# plot the relationship
dstp_shps_muT_plot <- all_data_dstp %>% 
  ggplot(aes(x = shps, y = muT)) + 
  geom_abline(intercept = dstp_shps_muT_samples$b_Intercept[1:200], 
              slope     = dstp_shps_muT_samples$b_shps[1:200],
              size = 1/3, alpha = .3, colour = "#ef8a62") + 
  geom_point(alpha = 0.6, colour = "#67a9cf", size = 2) + 
  labs(x = "SHPS Score", y = "Parameter  Value") + 
  ggtitle("Parameter: muT")



# parameter muFl
dstp_shps_muFl <- brm(muFl ~ shps, family = skew_normal(), 
                     data = all_data_dstp, cores = 4)               

print(summary(dstp_shps_muFl), digits = 4)

# Population-Level Effects: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# Intercept   0.2124    0.0184   0.1768   0.2476       4069 0.9996
# shps       -0.0007    0.0006  -0.0019   0.0005       4786 0.9992
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# sigma   0.1083    0.0049   0.0992   0.1186       1558 1.0009
# alpha   4.6948    0.9050   3.1651   6.6417       1774 1.0020

# get the posterior samples & plot them (but not included in the paper)
dstp_shps_muFl_samples <-  posterior_samples(dstp_shps_muFl, 
                                            pars = c("b_Intercept", 
                                                     "b_shps"))

# plot the relationship
dstp_shps_muFl_plot <- all_data_dstp %>% 
  ggplot(aes(x = shps, y = muFl)) + 
  geom_abline(intercept = dstp_shps_muFl_samples$b_Intercept[1:200], 
              slope     = dstp_shps_muFl_samples$b_shps[1:200],
              size = 1/3, alpha = .3, colour = "#ef8a62") + 
  geom_point(alpha = 0.6, colour = "#67a9cf", size = 2) + 
  labs(x = "SHPS Score", y = "Parameter  Value") + 
  ggtitle("Parameter: muFl")



# parameter muSS
dstp_shps_muSS <- brm(muSS ~ shps, family = skew_normal(), 
                      data = all_data_dstp, cores = 4)               

print(summary(dstp_shps_muSS), digits = 4)

# Population-Level Effects: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# Intercept   0.5286    0.0294   0.4717   0.5857       4299 1.0005
# shps       -0.0000    0.0010  -0.0020   0.0019       4663 1.0002
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# sigma   0.1373    0.0058   0.1267   0.1491       2230 1.0003
# alpha   1.5418    0.4696   0.4443   2.3634       2001 1.0038


# get the posterior samples & plot them (but not included in the paper)
dstp_shps_muSS_samples <-  posterior_samples(dstp_shps_muSS, 
                                             pars = c("b_Intercept", 
                                                      "b_shps"))

# plot the relationship
dstp_shps_muSS_plot <- all_data_dstp %>% 
  ggplot(aes(x = shps, y = muSS)) + 
  geom_abline(intercept = dstp_shps_muSS_samples$b_Intercept[1:200], 
              slope     = dstp_shps_muSS_samples$b_shps[1:200],
              size = 1/3, alpha = .3, colour = "#ef8a62") + 
  geom_point(alpha = 0.6, colour = "#67a9cf", size = 2) + 
  labs(x = "SHPS Score", y = "Parameter  Value") + 
  ggtitle("Parameter: muSS")



# parameter muRS2
dstp_shps_muRS2 <- brm(muRS2 ~ shps, family = skew_normal(), 
                      data = all_data_dstp, cores = 4)               

print(summary(dstp_shps_muRS2), digits = 4)

# Population-Level Effects: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# Intercept   1.2096    0.0546   1.1020   1.3209       5166 0.9995
# shps        0.0014    0.0019  -0.0025   0.0051       5601 0.9995
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# sigma   0.2540    0.0113   0.2328   0.2776       2761 1.0014
# alpha   1.5690    0.6899  -0.2565   2.7453       1705 1.0013

# get the posterior samples & plot them (but not included in the paper)
dstp_shps_muRS2_samples <-  posterior_samples(dstp_shps_muRS2, 
                                             pars = c("b_Intercept", 
                                                      "b_shps"))

# plot the relationship
dstp_shps_muRS2_plot <- all_data_dstp %>% 
  ggplot(aes(x = shps, y = muRS2)) + 
  geom_abline(intercept = dstp_shps_muRS2_samples$b_Intercept[1:200], 
              slope     = dstp_shps_muRS2_samples$b_shps[1:200],
              size = 1/3, alpha = .3, colour = "#ef8a62") + 
  geom_point(alpha = 0.6, colour = "#67a9cf", size = 2) + 
  labs(x = "SHPS Score", y = "Parameter  Value") + 
  ggtitle("Parameter: muRS2")


# parameter ter
dstp_shps_ter <- brm(ter ~ shps, family = skew_normal(), 
                       data = all_data_dstp, cores = 4)               

print(summary(dstp_shps_ter), digits = 4)

# Population-Level Effects: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# Intercept   0.2971    0.0098   0.2778   0.3165       4111 1.0013
# shps       -0.0003    0.0003  -0.0010   0.0004       4275 1.0007
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# sigma   0.0448    0.0018   0.0415   0.0487       1859 1.0000
# alpha  -0.2952    0.7302  -1.5216   1.0502       1438 1.0048


# get the posterior samples & plot them (but not included in the paper)
dstp_shps_ter_samples <-  posterior_samples(dstp_shps_ter, 
                                              pars = c("b_Intercept", 
                                                       "b_shps"))

# plot the relationship
dstp_shps_ter_plot <- all_data_dstp %>% 
  ggplot(aes(x = shps, y = ter)) + 
  geom_abline(intercept = dstp_shps_ter_samples$b_Intercept[1:200], 
              slope     = dstp_shps_ter_samples$b_shps[1:200],
              size = 1/3, alpha = .3, colour = "#ef8a62") + 
  geom_point(alpha = 0.6, colour = "#67a9cf", size = 2) + 
  labs(x = "SHPS Score", y = "Parameter  Value") + 
  ggtitle("Parameter: ter")




# draw all plots
setwd(here("paper_figures"))
ggsave("dstp_shps.png", 
       grid.arrange(dstp_shps_a_plot, dstp_shps_c_plot, dstp_shps_muT_plot, 
                    dstp_shps_muFl_plot, dstp_shps_muSS_plot, dstp_shps_muRS2_plot, 
                    dstp_shps_ter_plot, ncol = 3), 
       width = 10, 
       height = 10, 
       dpi = 1200)
setwd(here())


setwd(here("paper_figures"))
pdf("dstp_shps.pdf", width = 10, height = 10)
grid.arrange(dstp_shps_a_plot, dstp_shps_c_plot, dstp_shps_muT_plot, 
             dstp_shps_muFl_plot, dstp_shps_muSS_plot, dstp_shps_muRS2_plot, 
             dstp_shps_ter_plot, ncol = 3)
dev.off()
setwd(here())

#======================================
# SSP QIDS
#======================================

# first, look at the distribution of all 5 parameters
par(mfrow = c(3, 2))
for(i in 1:5){
  plot(density(all_data_ssp[, i + 1]), 
       main = colnames(all_data_ssp[i + 1]))
}


# parameter A
ssp_qids_a <- brm(a ~ qids, family = skew_normal(), 
                   data = all_data_ssp, cores = 4)               

print(summary(ssp_qids_a), digits = 4)

# Population-Level Effects: 
#         Estimate Est.Error l-95% CI u-95% CI    Eff.Sample   Rhat
# Intercept   0.0711    0.0021   0.0670   0.0754       2819 1.0007
# qids       -0.0001    0.0001  -0.0004   0.0002       3790 1.0009
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# sigma   0.0237    0.0010   0.0218   0.0256       1152 1.0027
# alpha   8.5469    1.6252   5.7808  12.2065       1002 1.0030

# get the posterior samples & plot them (but not included in the paper)
ssp_qids_a_samples <-  posterior_samples(ssp_qids_a, pars = c("b_Intercept", 
                                                                "b_qids"))

# plot the relationship
ssp_qids_a_plot <- all_data_ssp %>% 
  ggplot(aes(x = qids, y = a)) + 
  geom_abline(intercept = ssp_qids_a_samples$b_Intercept[1:200], 
              slope     = ssp_qids_a_samples$b_qids[1:200],
              size = 1/3, alpha = .3, colour = "#ef8a62") + 
  geom_point(alpha = 0.6, colour = "#67a9cf", size = 2) + 
  labs(x = "QIDS Score", y = "Parameter  Value") + 
  ggtitle("Parameter: A")


# parameter ter
ssp_qids_ter <- brm(ter ~ qids, family = skew_normal(), 
                  data = all_data_ssp, cores = 4)               

print(summary(ssp_qids_ter), digits = 4)

# Population-Level Effects: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# Intercept   0.3449    0.0063   0.3325   0.3569       3872 0.9992
# qids        0.0001    0.0005  -0.0009   0.0010       4404 0.9993
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# sigma   0.0505    0.0021   0.0467   0.0547       1711 1.0007
# alpha   0.9103    0.5367  -0.5484   1.6186        582 1.0019

# get the posterior samples & plot them (but not included in the paper)
ssp_qids_ter_samples <-  posterior_samples(ssp_qids_ter, pars = c("b_Intercept", 
                                                              "b_qids"))

# plot the relationship
ssp_qids_ter_plot <- all_data_ssp %>% 
  ggplot(aes(x = qids, y = ter)) + 
  geom_abline(intercept = ssp_qids_ter_samples$b_Intercept[1:200], 
              slope     = ssp_qids_ter_samples$b_qids[1:200],
              size = 1/3, alpha = .3, colour = "#ef8a62") + 
  geom_point(alpha = 0.6, colour = "#67a9cf", size = 2) + 
  labs(x = "QIDS Score", y = "Parameter  Value") + 
  ggtitle("Parameter: ter")


# parameter p
ssp_qids_p <- brm(p ~ qids, family = skew_normal(), 
                    data = all_data_ssp, cores = 4)               

print(summary(ssp_qids_p), digits = 4)

# Population-Level Effects: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# Intercept   0.5021    0.0142   0.4740   0.5303       4554 1.0000
# qids       -0.0021    0.0011  -0.0043   0.0001       5095 1.0000
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# sigma   0.1206    0.0053   0.1108   0.1312       2260 1.0002
# alpha   1.7679    0.4273   0.9484   2.6275       2358 1.0007

# get the posterior samples & plot them (but not included in the paper)
ssp_qids_p_samples <-  posterior_samples(ssp_qids_p, pars = c("b_Intercept", 
                                                                  "b_qids"))

# plot the relationship
ssp_qids_p_plot <- all_data_ssp %>% 
  ggplot(aes(x = qids, y = p)) + 
  geom_abline(intercept = ssp_qids_p_samples$b_Intercept[1:200], 
              slope     = ssp_qids_p_samples$b_qids[1:200],
              size = 1/3, alpha = .3, colour = "#ef8a62") + 
  geom_point(alpha = 0.6, colour = "#67a9cf", size = 2) + 
  labs(x = "QIDS Score", y = "Parameter  Value") + 
  ggtitle("Parameter: p")


# parameter rd
ssp_qids_rd <- brm(rd ~ qids, family = skew_normal(), 
                  data = all_data_ssp, cores = 4)               

print(summary(ssp_qids_rd), digits = 4)

# Population-Level Effects: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# Intercept   0.0364    0.0017   0.0332   0.0397       2319 0.9996
# qids       -0.0002    0.0001  -0.0004   0.0000       3921 0.9999
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# sigma   0.0204    0.0009   0.0187   0.0222       1015 1.0040
# alpha  10.8643    1.9555   7.4440  15.3489       1230 1.0029


# get the posterior samples & plot them (but not included in the paper)
ssp_qids_rd_samples <-  posterior_samples(ssp_qids_rd, pars = c("b_Intercept", 
                                                              "b_qids"))

# plot the relationship
ssp_qids_rd_plot <- all_data_ssp %>% 
  ggplot(aes(x = qids, y = rd)) + 
  geom_abline(intercept = ssp_qids_rd_samples$b_Intercept[1:200], 
              slope     = ssp_qids_rd_samples$b_qids[1:200],
              size = 1/3, alpha = .3, colour = "#ef8a62") + 
  geom_point(alpha = 0.6, colour = "#67a9cf", size = 2) + 
  labs(x = "QIDS Score", y = "Parameter  Value") + 
  ggtitle("Parameter: rd")


# parameter sda
ssp_qids_sda <- brm(sda ~ qids, family = skew_normal(), 
                   data = all_data_ssp, cores = 4)               

print(summary(ssp_qids_sda), digits = 4)

# Population-Level Effects: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# Intercept   1.7708    0.0397   1.6921   1.8491       4408 0.9996
# qids       -0.0023    0.0031  -0.0083   0.0038       4830 0.9994
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# sigma   0.3354    0.0149   0.3070   0.3654       2976 0.9999
# alpha   1.2106    0.7111  -0.5077   2.3413       2032 0.9996


# get the posterior samples & plot them (but not included in the paper)
ssp_qids_sda_samples <-  posterior_samples(ssp_qids_sda, pars = c("b_Intercept", 
                                                                "b_qids"))

# plot the relationship
ssp_qids_sda_plot <- all_data_ssp %>% 
  ggplot(aes(x = qids, y = sda)) + 
  geom_abline(intercept = ssp_qids_sda_samples$b_Intercept[1:200], 
              slope     = ssp_qids_sda_samples$b_qids[1:200],
              size = 1/3, alpha = .3, colour = "#ef8a62") + 
  geom_point(alpha = 0.6, colour = "#67a9cf", size = 2) + 
  labs(x = "QIDS Score", y = "Parameter  Value") + 
  ggtitle("Parameter: sda")

# draw all plots
setwd(here("paper_figures"))
ggsave("ssp_qids.png", 
       grid.arrange(ssp_qids_a_plot, ssp_qids_ter_plot, ssp_qids_p_plot, 
                    ssp_qids_rd_plot, ssp_qids_sda_plot, ncol = 3), 
       width = 10, 
       height = 7, 
       dpi = 1200)
setwd(here())


setwd(here("paper_figures"))
pdf("ssp_qids.pdf", width = 10, height = 7)
grid.arrange(ssp_qids_a_plot, ssp_qids_ter_plot, ssp_qids_p_plot, 
             ssp_qids_rd_plot, ssp_qids_sda_plot, ncol = 3)
dev.off()
setwd(here())


#======================================
# SSP SHPS
#======================================

# first, look at the distribution of all 5 parameters
par(mfrow = c(3, 2))
for(i in 1:5){
  plot(density(all_data_ssp[, i + 1]), 
       main = colnames(all_data_ssp[i + 1]))
}


# parameter A
ssp_shps_a <- brm(a ~ shps, family = skew_normal(), 
                  data = all_data_ssp, cores = 4)               

print(summary(ssp_shps_a), digits = 4)

# Population-Level Effects: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# Intercept   0.0716    0.0032   0.0653   0.0783       3282 1.0005
# shps       -0.0001    0.0001  -0.0003   0.0002       3716 1.0001
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# sigma   0.0236    0.0010   0.0217   0.0258       1354 1.0014
# alpha   8.5436    1.5365   5.9261  11.8015       1356 1.0006

# get the posterior samples & plot them (but not included in the paper)
ssp_shps_a_samples <-  posterior_samples(ssp_shps_a, pars = c("b_Intercept", 
                                                              "b_shps"))

# plot the relationship
ssp_shps_a_plot <- all_data_ssp %>% 
  ggplot(aes(x = shps, y = a)) + 
  geom_abline(intercept = ssp_shps_a_samples$b_Intercept[1:200], 
              slope     = ssp_shps_a_samples$b_shps[1:200],
              size = 1/3, alpha = .3, colour = "#ef8a62") + 
  geom_point(alpha = 0.6, colour = "#67a9cf", size = 2) + 
  labs(x = "SHPS Score", y = "Parameter  Value") + 
  ggtitle("Parameter: A")



# parameter ter
ssp_shps_ter <- brm(ter ~ shps, family = skew_normal(), 
                  data = all_data_ssp, cores = 4)               

print(summary(ssp_shps_ter), digits = 4)

# Population-Level Effects: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# Intercept   0.3559    0.0113   0.3342   0.3777       3920 1.0021
# shps       -0.0004    0.0004  -0.0011   0.0004       4402 1.0017
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# sigma   0.0505    0.0022   0.0465   0.0549       1082 1.0016
# alpha   1.0051    0.4305  -0.3072   1.6058        963 1.0025


# get the posterior samples & plot them (but not included in the paper)
ssp_shps_ter_samples <-  posterior_samples(ssp_shps_ter, pars = c("b_Intercept", 
                                                              "b_shps"))

# plot the relationship
ssp_shps_ter_plot <- all_data_ssp %>% 
  ggplot(aes(x = shps, y = ter)) + 
  geom_abline(intercept = ssp_shps_ter_samples$b_Intercept[1:200], 
              slope     = ssp_shps_ter_samples$b_shps[1:200],
              size = 1/3, alpha = .3, colour = "#ef8a62") + 
  geom_point(alpha = 0.6, colour = "#67a9cf", size = 2) + 
  labs(x = "SHPS Score", y = "Parameter  Value") + 
  ggtitle("Parameter: ter")


# parameter p
ssp_shps_p <- brm(p ~ shps, family = skew_normal(), 
                    data = all_data_ssp, cores = 4)               

print(summary(ssp_shps_p), digits = 4)

# Population-Level Effects: 
# Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# Intercept   0.4893    0.0255   0.4391   0.5389       4417 0.9997
# shps       -0.0004    0.0009  -0.0021   0.0013       4553 0.9996
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# sigma   0.1214    0.0053   0.1117   0.1321       2415 0.9997
# alpha   1.8113    0.4410   1.0201   2.6842       1851 1.0010


# get the posterior samples & plot them (but not included in the paper)
ssp_shps_p_samples <-  posterior_samples(ssp_shps_p, pars = c("b_Intercept", 
                                                                  "b_shps"))

# plot the relationship
ssp_shps_p_plot <- all_data_ssp %>% 
  ggplot(aes(x = shps, y = p)) + 
  geom_abline(intercept = ssp_shps_p_samples$b_Intercept[1:200], 
              slope     = ssp_shps_p_samples$b_shps[1:200],
              size = 1/3, alpha = .3, colour = "#ef8a62") + 
  geom_point(alpha = 0.6, colour = "#67a9cf", size = 2) + 
  labs(x = "SHPS Score", y = "Parameter  Value") + 
  ggtitle("Parameter: p")



# parameter rd
ssp_shps_rd <- brm(rd ~ shps, family = skew_normal(), 
                  data = all_data_ssp, cores = 4)               

print(summary(ssp_shps_rd), digits = 4)

# Population-Level Effects: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# Intercept   0.0368    0.0029   0.0312   0.0423       3211 0.9995
# shps       -0.0001    0.0001  -0.0003   0.0001       4241 0.9995
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# sigma   0.0204    0.0009   0.0188   0.0222       1276 1.0003
# alpha  10.7363    1.8954   7.4406  14.8807       1250 1.0012


# get the posterior samples & plot them (but not included in the paper)
ssp_shps_rd_samples <-  posterior_samples(ssp_shps_rd, pars = c("b_Intercept", 
                                                              "b_shps"))

# plot the relationship
ssp_shps_rd_plot <- all_data_ssp %>% 
  ggplot(aes(x = shps, y = rd)) + 
  geom_abline(intercept = ssp_shps_rd_samples$b_Intercept[1:200], 
              slope     = ssp_shps_rd_samples$b_shps[1:200],
              size = 1/3, alpha = .3, colour = "#ef8a62") + 
  geom_point(alpha = 0.6, colour = "#67a9cf", size = 2) + 
  labs(x = "SHPS Score", y = "Parameter  Value") + 
  ggtitle("Parameter: rd")


# parameter sda
ssp_shps_sda <- brm(sda ~ shps, family = skew_normal(), 
                   data = all_data_ssp, cores = 4)               

print(summary(ssp_shps_sda), digits = 4)

# Population-Level Effects: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# Intercept   1.8068    0.0718   1.6655   1.9467       4975 1.0006
# shps       -0.0022    0.0025  -0.0069   0.0026       5181 1.0003
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# sigma   0.3359    0.0145   0.3094   0.3663       2884 0.9997
# alpha   1.2938    0.7314  -0.4652   2.4812       2021 1.0010




# get the posterior samples & plot them (but not included in the paper)
ssp_shps_sda_samples <-  posterior_samples(ssp_shps_sda, pars = c("b_Intercept", 
                                                                "b_shps"))

# plot the relationship
ssp_shps_sda_plot <- all_data_ssp %>% 
  ggplot(aes(x = shps, y = sda)) + 
  geom_abline(intercept = ssp_shps_sda_samples$b_Intercept[1:200], 
              slope     = ssp_shps_sda_samples$b_shps[1:200],
              size = 1/3, alpha = .3, colour = "#ef8a62") + 
  geom_point(alpha = 0.6, colour = "#67a9cf", size = 2) + 
  labs(x = "SHPS Score", y = "Parameter  Value") + 
  ggtitle("Parameter: sda")


# draw all plots
setwd(here("paper_figures"))
ggsave("ssp_shps.png", 
       grid.arrange(ssp_shps_a_plot, ssp_shps_ter_plot, ssp_shps_p_plot, 
                    ssp_shps_rd_plot, ssp_shps_sda_plot, ncol = 3), 
       width = 10, 
       height = 7, 
       dpi = 1200)
setwd(here())


setwd(here("paper_figures"))
pdf("ssp_shps.pdf", width = 10, height = 7)
grid.arrange(ssp_shps_a_plot, ssp_shps_ter_plot, ssp_shps_p_plot, 
             ssp_shps_rd_plot, ssp_shps_sda_plot, ncol = 3)
dev.off()
setwd(here())


#==============================================================================
# group participants into high QIDS (>=14) and low QIDS (<8) and do group 
# comparisons on all main DVs
#==============================================================================

# sort the participants
high_qids_dstp <- all_data_dstp %>% 
  filter(qids >= 14) %>% 
  mutate(Group = "High")
low_qids_dstp <- all_data_dstp %>% 
  filter(qids < 8) %>% 
  mutate(Group = "Low")

qids_dstp <- rbind(high_qids_dstp, low_qids_dstp)
qids_dstp$Group <- as.factor(qids_dstp$Group)
 

high_qids_ssp <- all_data_ssp %>% 
  filter(qids >= 14) %>% 
  mutate(Group = "High")
low_qids_ssp <- all_data_ssp %>% 
  filter(qids < 8) %>% 
  mutate(Group = "Low")

qids_ssp <- rbind(high_qids_ssp, low_qids_ssp) 
qids_ssp$Group <- as.factor(qids_ssp$Group)



#======================================
# behavioural data
#======================================

# flanker rt
qids_flanker_rt <- brm(rt ~ Group, data = qids_dstp, 
                       cores = 4)

# flanker accuracy
qids_flanker_acc <- brm(accuracy ~ Group, data = qids_dstp, 
                        cores = 4, family = skew_normal())

# mean RT
qids_mean_rt <- brm(mean_rt ~ Group, data = qids_dstp, 
                    cores = 4, family = exgaussian())

# mean accuracy
qids_mean_acc <- brm(mean_accuracy ~ Group, data = qids_dstp, 
                     cores = 4, family = skew_normal())


## do all of the density plots

# flanker rt
qids_flanker_rt_plot <- qids_dstp %>% 
  ggplot(aes(x = rt, colour = Group)) + 
  geom_density(aes(colour = Group, fill = Group), 
               alpha = 0.8) + 
  scale_colour_manual(values = c("#3182bd", "#9ecae1")) +
  scale_fill_manual(values = c("#3182bd", "#9ecae1")) +
  labs(x = "RT Flanker Effect (ms)") 


# flanker accuracy
qids_flanker_acc_plot <- qids_dstp %>% 
  ggplot(aes(x = accuracy, colour = Group)) + 
  geom_density(aes(colour = Group, fill = Group), 
               alpha = 0.8) + 
  scale_colour_manual(values = c("#3182bd", "#9ecae1")) +
  scale_fill_manual(values = c("#3182bd", "#9ecae1")) +
  labs(x = "Accuracy Flanker Effect") 

# mean rt
qids_mean_rt_plot <- qids_dstp %>% 
  ggplot(aes(x = mean_rt, colour = Group)) + 
  geom_density(aes(colour = Group, fill = Group), 
               alpha = 0.8) + 
  scale_colour_manual(values = c("#3182bd", "#9ecae1")) +
  scale_fill_manual(values = c("#3182bd", "#9ecae1")) +
  labs(x = "Mean Response Time (ms)") 

# mean accuracy
qids_mean_acc_plot <- qids_dstp %>% 
  ggplot(aes(x = mean_accuracy, colour = Group)) + 
  geom_density(aes(colour = Group, fill = Group), 
               alpha = 0.8) + 
  scale_colour_manual(values = c("#3182bd", "#9ecae1")) +
  scale_fill_manual(values = c("#3182bd", "#9ecae1")) +
  labs(x = "Mean Accuracy") 


setwd(here("paper_figures"))
pdf("extreme_qids_behavioural_plots.pdf", width = 10, height = 8)
grid.arrange(qids_mean_rt_plot, qids_mean_acc_plot, qids_flanker_rt_plot, 
             qids_flanker_acc_plot, ncol = 2)
dev.off()
setwd(here())



#======================================
# CDFS & CAFs
#======================================

setwd(here("model_data"))
all_data_dstp = read.csv("all_data_dstp.csv")
all_data_ssp = read.csv("all_data_ssp.csv")


# get a list of the ids with high qids
high_qids_ids <- all_data_dstp %>% 
  filter(qids >= 14) %>% 
  pull(id) %>% 
  as.character()


# get a list of the ids with low qids
low_qids_ids <- all_data_dstp %>% 
  filter(qids < 8) %>% 
  pull(id) %>% 
  as.character()


# cdfs
probs <- c(0.1, 0.3, 0.5, 0.7, 0.9)

cdf_high_qids_con <- flanker_data %>% 
  filter(id %in% high_qids_ids) %>% 
  filter(accuracy == 1, congruency == "congruent") %>% 
  rename(subject = id) %>% 
  flankr::cdf(.)

cdf_high_qids_incon <- flanker_data %>% 
  filter(id %in% high_qids_ids) %>% 
  filter(accuracy == 1, congruency == "incongruent") %>% 
  rename(subject = id) %>% 
  flankr::cdf(.)

cdf_low_qids_con <- flanker_data %>% 
  filter(id %in% low_qids_ids) %>% 
  filter(accuracy == 1, congruency == "congruent") %>% 
  rename(subject = id) %>% 
  flankr::cdf(.)

cdf_low_qids_incon <- flanker_data %>% 
  filter(id %in% low_qids_ids) %>% 
  filter(accuracy == 1, congruency == "incongruent") %>% 
  rename(subject = id) %>% 
  flankr::cdf(.)

# cafs
caf_high_qids_con <- flanker_data %>% 
  filter(id %in% high_qids_ids) %>% 
  filter(congruency == "congruent") %>% 
  rename(subject = id) %>% 
  flankr::caf(.)

caf_high_qids_incon <- flanker_data %>% 
  filter(id %in% high_qids_ids) %>% 
  filter(congruency == "incongruent") %>% 
  rename(subject = id) %>% 
  flankr::caf(.)

caf_low_qids_con <- flanker_data %>% 
  filter(id %in% low_qids_ids) %>% 
  filter(congruency == "congruent") %>% 
  rename(subject = id) %>% 
  flankr::caf(.)

caf_low_qids_incon <- flanker_data %>% 
  filter(id %in% low_qids_ids) %>% 
  filter(congruency == "incongruent") %>% 
  rename(subject = id) %>% 
  flankr::caf(.)



## plotting
setwd(here("paper_figures"))

png("cdf_caf_plots.png", width = 8, height = 8, units = "in", res = 1200)
par(mfrow = c(1, 2))
do_clean_graphs()

# plot the cdfs
min_rt = min(cdf_high_qids_con, cdf_high_qids_incon, 
             cdf_low_qids_con, cdf_low_qids_incon)

max_rt = max(cdf_high_qids_con, cdf_high_qids_incon, 
             cdf_low_qids_con, cdf_low_qids_incon)


plot(cdf_high_qids_con, probs, type = "b", xlim = c(min_rt - 50, max_rt + 50), 
     ylim = c(0, 1), pch = 19, col = "#E69F00", 
     xlab = "Response Time (ms)", ylab = "Cumulative Probability", cex = 1.5, 
     lwd = 2)
lines(cdf_high_qids_incon, probs, type = "b", pch = 1, lty = 2, col = "#E69F00", 
      cex = 1.5, lwd = 2)
lines(cdf_low_qids_con, probs, type = "b", pch = 15, col = "#56B4E9", 
      cex = 1.5, lwd = 2)
lines(cdf_low_qids_incon, probs, type = "b", pch = 0, lty = 2, col = "#56B4E9", 
      cex = 1.5, lwd = 2)

# plot the cafs
min_rt = min(caf_high_qids_con[1, ], caf_high_qids_incon[1, ], 
             caf_low_qids_con[1, ], caf_low_qids_incon[1, ])
max_rt = max(caf_high_qids_con[1, ], caf_high_qids_incon[1, ], 
             caf_low_qids_con[1, ], caf_low_qids_incon[1, ])

plot(caf_high_qids_con[1, ], caf_high_qids_con[2, ], type = "b", 
     xlim = c(min_rt - 50, max_rt + 50), ylim = c(.7, 1), pch = 19, col = "#E69F00", 
     xlab = "Response Time (ms)", ylab = "Proportion Error", cex = 1.5, lwd = 2)
lines(caf_high_qids_incon[1, ], caf_high_qids_incon[2, ], type = "b", 
      pch = 1, lty = 2, col = "#E69F00", cex = 1.5, lwd = 2)
lines(caf_low_qids_con[1, ], caf_low_qids_con[2, ], type = "b", 
      pch = 15, col = "#56B4E9", cex = 1.5, lwd = 2)
lines(caf_low_qids_incon[1, ], caf_low_qids_incon[2, ], type = "b", 
      pch = 0, lty = 2, col = "#56B4E9", cex = 1.5, lwd = 2)

legend("bottomright", c("High QIDS, Congruent","High QIDS, Incongruent", 
                        "Low QIDS, Congruent", "Low QIDS, Incongruent"), 
       pch = c(19, 1, 15, 0),
       lty = c(1, 2, 1, 2), 
       lwd = 1.5,
       col = c("#E69F00", "#E69F00", "#56B4E9", "#56B4E9"), bty="n");
dev.off()


#======================================
# dstp parameters
#======================================

# A
qids_a <- brm(A ~ Group, data = qids_dstp, cores = 4, 
              family = skew_normal())
print(summary(qids_a), digits = 4)

# C
qids_c<- brm(C ~ Group, data = qids_dstp, cores = 4, 
              family = skew_normal())
print(summary(qids_c), digits = 4)

# muT
qids_muT <- brm(muT ~ Group, data = qids_dstp, cores = 4, 
                family = skew_normal())
print(summary(qids_muT), digits = 4)

# muFl
qids_muFl <- brm(muFl ~ Group, data = qids_dstp, cores = 4, 
                 family = skew_normal())
print(summary(qids_muFl), digits = 4)

# muSS
qids_muSS <- brm(muSS ~ Group, data = qids_dstp, cores = 4, 
                 family = skew_normal())
print(summary(qids_muSS), digits = 4)

# muRS2
qids_muRS2 <- brm(muRS2 ~ Group, data = qids_dstp, cores = 4, 
                  family = skew_normal())
print(summary(qids_muRS2), digits = 4)

# ter
qids_ter <- brm(ter ~ Group, data = qids_dstp, cores = 4, 
                family = skew_normal())
print(summary(qids_ter), digits = 4)


## do all of the density plots

# A
qids_a_plot <- qids_dstp %>% 
  ggplot(aes(x = A, colour = Group)) + 
  geom_density(aes(colour = Group, fill = Group), 
               alpha = 0.8) + 
  scale_colour_manual(values = c("#3182bd", "#9ecae1")) +
  scale_fill_manual(values = c("#3182bd", "#9ecae1")) +
  labs(x = "Parameter A") 

# C
qids_c_plot <- qids_dstp %>% 
  ggplot(aes(x = C, colour = Group)) + 
  geom_density(aes(colour = Group, fill = Group), 
               alpha = 0.8) + 
  scale_colour_manual(values = c("#3182bd", "#9ecae1")) +
  scale_fill_manual(values = c("#3182bd", "#9ecae1")) +
  labs(x = "Parameter C") 

# mu_t
qids_muT_plot <- qids_dstp %>% 
  ggplot(aes(x = muT, colour = Group)) + 
  geom_density(aes(colour = Group, fill = Group), 
               alpha = 0.8) + 
  scale_colour_manual(values = c("#3182bd", "#9ecae1")) +
  scale_fill_manual(values = c("#3182bd", "#9ecae1")) +
  labs(x = "Parameter muT") 

# mu_fl
qids_muFl_plot <- qids_dstp %>% 
  ggplot(aes(x = muFl, colour = Group)) + 
  geom_density(aes(colour = Group, fill = Group), 
               alpha = 0.8) + 
  scale_colour_manual(values = c("#3182bd", "#9ecae1")) +
  scale_fill_manual(values = c("#3182bd", "#9ecae1")) +
  labs(x = "Parameter muFl") 

# mu_SS
qids_muSS_plot <- qids_dstp %>% 
  ggplot(aes(x = muSS, colour = Group)) + 
  geom_density(aes(colour = Group, fill = Group), 
               alpha = 0.8) + 
  scale_colour_manual(values = c("#3182bd", "#9ecae1")) +
  scale_fill_manual(values = c("#3182bd", "#9ecae1")) +
  labs(x = "Parameter muSS") 

# mu_RS2
qids_muRS2_plot <- qids_dstp %>% 
  ggplot(aes(x = muRS2, colour = Group)) + 
  geom_density(aes(colour = Group, fill = Group), 
               alpha = 0.8) + 
  scale_colour_manual(values = c("#3182bd", "#9ecae1")) +
  scale_fill_manual(values = c("#3182bd", "#9ecae1")) +
  labs(x = "Parameter muRS2") 

# ter
qids_ter_plot <- qids_dstp %>% 
  ggplot(aes(x = ter, colour = Group)) + 
  geom_density(aes(colour = Group, fill = Group), 
               alpha = 0.8) + 
  scale_colour_manual(values = c("#3182bd", "#9ecae1")) +
  scale_fill_manual(values = c("#3182bd", "#9ecae1")) +
  labs(x = "Parameter ter") 


setwd(here("paper_figures"))
pdf("extreme_qids_dstp_plots.pdf", width = 10, height = 8)
grid.arrange(qids_a_plot, qids_c_plot, qids_muT_plot,
             qids_muFl_plot, qids_muSS_plot, qids_muRS2_plot, 
             qids_ter_plot, ncol = 3)
dev.off()
setwd(here())



#======================================
# ssp parameters
#======================================

# a
qids_a <- brm(a ~ Group, data = qids_ssp, cores = 4, 
              family = skew_normal())
print(summary(qids_a), digits = 4)

# ter
qids_ter <- brm(ter ~ Group, data = qids_ssp, cores = 4, 
              family = skew_normal())
print(summary(qids_ter), digits = 4)

# p
qids_p <- brm(p ~ Group, data = qids_ssp, cores = 4, 
                family = skew_normal())
print(summary(qids_p), digits = 4)

# rd
qids_rd <- brm(rd ~ Group, data = qids_ssp, cores = 4, 
              family = skew_normal())
print(summary(qids_rd), digits = 4)

# sda
qids_sda <- brm(sda ~ Group, data = qids_ssp, cores = 4, 
              family = skew_normal())
print(summary(qids_sda), digits = 4)


## do all of the density plots

# a
qids_a_plot <- qids_ssp %>% 
  ggplot(aes(x = a, colour = Group)) + 
  geom_density(aes(colour = Group, fill = Group), 
               alpha = 0.8) + 
  scale_colour_manual(values = c("#3182bd", "#9ecae1")) +
  scale_fill_manual(values = c("#3182bd", "#9ecae1")) +
  labs(x = "Parameter a") 

# ter
qids_ter_plot <- qids_ssp %>% 
  ggplot(aes(x = ter, colour = Group)) + 
  geom_density(aes(colour = Group, fill = Group), 
               alpha = 0.8) + 
  scale_colour_manual(values = c("#3182bd", "#9ecae1")) +
  scale_fill_manual(values = c("#3182bd", "#9ecae1")) +
  labs(x = "Parameter ter") 

# p
qids_p_plot <- qids_ssp %>% 
  ggplot(aes(x = p, colour = Group)) + 
  geom_density(aes(colour = Group, fill = Group), 
               alpha = 0.8) + 
  scale_colour_manual(values = c("#3182bd", "#9ecae1")) +
  scale_fill_manual(values = c("#3182bd", "#9ecae1")) +
  labs(x = "Parameter p") 

# rd
qids_rd_plot <- qids_ssp %>% 
  ggplot(aes(x = rd, colour = Group)) + 
  geom_density(aes(colour = Group, fill = Group), 
               alpha = 0.8) + 
  scale_colour_manual(values = c("#3182bd", "#9ecae1")) +
  scale_fill_manual(values = c("#3182bd", "#9ecae1")) +
  labs(x = "Parameter rd") 

# sda
qids_sda_plot <- qids_ssp %>% 
  ggplot(aes(x = sda, colour = Group)) + 
  geom_density(aes(colour = Group, fill = Group), 
               alpha = 0.8) + 
  scale_colour_manual(values = c("#3182bd", "#9ecae1")) +
  scale_fill_manual(values = c("#3182bd", "#9ecae1")) +
  labs(x = "Parameter sda") 


setwd(here("paper_figures"))
pdf("extreme_qids_ssp_plots.pdf", width = 10, height = 8)
grid.arrange(qids_a_plot, qids_ter_plot, qids_p_plot,
             qids_rd_plot, qids_sda_plot, ncol = 3)
dev.off()
setwd(here())






#==============================================================================
# find those with a clinical diagnosis of depression & see whether they
# perform differently from those without (i.e., similar to Dillon's design)
#==============================================================================

# get a list of the ids of participants who have a depression diagnosis
depression_diagnosis_ids <- demographics_data %>% 
  filter(Question.Key == "dem_depression") %>% 
  filter(Response == "Yes") %>% 
  pull(Participant.Public.ID) %>% 
  as.character()


write.csv(depression_diagnosis_ids, "e1_depression_diagnosis_ids.csv")

#======================================
# QIDS scores
#======================================

qids_with_dep <- all_data %>% 
  filter(id %in% depression_diagnosis_ids) %>% 
  mutate(Group = "with_depression") %>% 
  select(id, qids, Group)

qids_without_dep <- all_data %>% 
  filter(id %!in% depression_diagnosis_ids) %>% 
  mutate(Group = "without_depression") %>% 
  select(id, qids, Group)

qids_depression <- rbind(qids_with_dep, qids_without_dep)

qids_depression_plot <- qids_depression %>% 
  ggplot(aes(x = qids, colour = Group)) + 
  geom_density(aes(colour = Group, fill = Group), 
               alpha = 0.8) + 
  scale_colour_manual(values = c("#3182bd", "#9ecae1")) +
  scale_fill_manual(values = c("#3182bd", "#9ecae1")) +
  labs(x = "QIDS Scores") 


# model the group difference
qids_depression_brms <- brm(qids ~ Group, data = qids_depression, 
                            cores = 4, family = skew_normal())

# Population-Level Effects: 
#                     Estimate     Est.Error l-95% CI u-95% CI     Eff.Sample Rhat
# Intercept                  15.31      0.59    14.16    16.46       3782 1.00
# Groupwithout_depression    -5.87      0.71    -7.25    -4.46       3648 1.00
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
# sigma     5.50      0.23     5.07     5.97       3440 1.00

#======================================
# SHPS scores
#======================================

shps_with_dep <- all_data %>% 
  filter(id %in% depression_diagnosis_ids) %>% 
  mutate(Group = "with_depression") %>% 
  select(id, shps, Group)

shps_without_dep <- all_data %>% 
  filter(id %!in% depression_diagnosis_ids) %>% 
  mutate(Group = "without_depression") %>% 
  select(id, shps, Group)

shps_depression <- rbind(shps_with_dep, shps_without_dep)

shps_depression_plot <- shps_depression %>% 
  ggplot(aes(x = shps, colour = Group)) + 
  geom_density(aes(colour = Group, fill = Group), 
               alpha = 0.8) + 
  scale_colour_manual(values = c("#3182bd", "#9ecae1")) +
  scale_fill_manual(values = c("#3182bd", "#9ecae1")) +
  labs(x = "SHPS Scores")

# model the group difference
shps_depression_brms <- brm(shps ~ Group, data = shps_depression, 
                            cores = 4)

# Population-Level Effects: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
# Intercept                  31.76      0.75    30.30    33.23       3764 1.00
# Groupwithout_depression    -5.57      0.91    -7.36    -3.79       3940 1.00
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
# sigma     7.28      0.30     6.72     7.91       3659 1.00


#======================================
# plot them!
#======================================
setwd(here("paper_figures"))
pdf("depression_diagnosis_questionnaires.pdf", width = 10, height = 4)
grid.arrange(qids_depression_plot, shps_depression_plot, ncol = 2)
dev.off()
setwd(here())


#======================================
# mean rt
#======================================
rt_with_dep <- all_data %>% 
  filter(id %in% depression_diagnosis_ids) %>% 
  mutate(Group = "with_depression") %>% 
  select(id, mean_rt, Group)

rt_without_dep <- all_data %>% 
  filter(id %!in% depression_diagnosis_ids) %>% 
  mutate(Group = "without_depression") %>% 
  select(id, mean_rt, Group)

mean_rt_depression <- rbind(rt_with_dep, rt_without_dep)

mean_rt_depression_plot <- mean_rt_depression %>% 
  ggplot(aes(x = mean_rt, colour = Group)) + 
  geom_density(aes(colour = Group, fill = Group), 
               alpha = 0.8) + 
  scale_colour_manual(values = c("#3182bd", "#9ecae1")) +
  scale_fill_manual(values = c("#3182bd", "#9ecae1")) +
  labs(x = "Mean Response Time (ms)")

# model the group difference
mean_rt_depression_brms <- brm(mean_rt ~ Group, data = mean_rt_depression, 
                               cores = 4, family = exgaussian())

# Population-Level Effects: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
# Intercept                 521.30      7.84   505.76   536.68       3726 1.00
# Groupwithout_depression     2.56      8.71   -14.54    20.00       3798 1.00
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
# sigma    45.00      4.27    36.82    53.73       3378 1.00
# beta     71.08      6.53    58.80    84.73       2813 1.00

#======================================
# flanker rt
#======================================

rt_with_dep <- all_data %>% 
  filter(id %in% depression_diagnosis_ids) %>% 
  mutate(Group = "with_depression") %>% 
  select(id, rt, Group)

rt_without_dep <- all_data %>% 
  filter(id %!in% depression_diagnosis_ids) %>% 
  mutate(Group = "without_depression") %>% 
  select(id, rt, Group)

rt_depression <- rbind(rt_with_dep, rt_without_dep)

rt_depression_plot <- rt_depression %>% 
  ggplot(aes(x = rt, colour = Group)) + 
  geom_density(aes(colour = Group, fill = Group), 
               alpha = 0.8) + 
  scale_colour_manual(values = c("#3182bd", "#9ecae1")) +
  scale_fill_manual(values = c("#3182bd", "#9ecae1")) +
  labs(x = "Flanker Effect (Response Time, ms)")

# model the group difference
rt_depression_brms <- brm(rt ~ Group, data = rt_depression, 
                          cores = 4)


# Population-Level Effects: 
#                        Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
# Intercept                  64.21      3.14    58.14    70.50       3449 1.00
# Groupwithout_depression    -4.18      3.84   -11.62     3.26       3218 1.00
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
# sigma    29.44      1.21    27.18    31.84       3801 1.00



#======================================
# mean accuracy
#======================================
acc_with_dep <- all_data %>% 
  filter(id %in% depression_diagnosis_ids) %>% 
  mutate(Group = "with_depression") %>% 
  select(id, mean_accuracy, Group)

acc_without_dep <- all_data %>% 
  filter(id %!in% depression_diagnosis_ids) %>% 
  mutate(Group = "without_depression") %>% 
  select(id, mean_accuracy, Group)

mean_acc_depression <- rbind(acc_with_dep, acc_without_dep)

mean_acc_depression_plot <- mean_acc_depression %>% 
  ggplot(aes(x = mean_accuracy, colour = Group)) + 
  geom_density(aes(colour = Group, fill = Group), 
               alpha = 0.8) + 
  scale_colour_manual(values = c("#3182bd", "#9ecae1")) +
  scale_fill_manual(values = c("#3182bd", "#9ecae1")) +
  labs(x = "Mean Accuracy (%)")

# model the group difference
acc_depression_brms <- brm(mean_accuracy ~ Group, 
                           data = mean_acc_depression, 
                           cores = 4, family = skew_normal())

# Population-Level Effects: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# Intercept                 0.9583    0.0022   0.9539   0.9626       2324 1.0004
# Groupwithout_depression   0.0009    0.0020  -0.0032   0.0048       3965 0.9992
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# sigma   0.0311    0.0013   0.0285   0.0338       1714 1.0007
# alpha -12.2508    2.2034 -17.0354  -8.3738       2039 1.0003


#======================================
# flanker accuracy
#======================================

acc_with_dep <- all_data %>% 
  filter(id %in% depression_diagnosis_ids) %>% 
  mutate(Group = "with_depression") %>% 
  select(id, accuracy, Group)

acc_without_dep <- all_data %>% 
  filter(id %!in% depression_diagnosis_ids) %>% 
  mutate(Group = "without_depression") %>% 
  select(id, accuracy, Group)

acc_depression <- rbind(acc_with_dep, acc_without_dep)

acc_depression_plot <- acc_depression %>% 
  ggplot(aes(x = accuracy, colour = Group)) + 
  geom_density(aes(colour = Group, fill = Group), 
               alpha = 0.8) + 
  scale_colour_manual(values = c("#3182bd", "#9ecae1")) +
  scale_fill_manual(values = c("#3182bd", "#9ecae1")) +
  labs(x = "Flanker Effect (Accuracy, %)")

# model the group difference
acc_depression_brms <- brm(accuracy ~ Group, data = acc_depression, 
                           cores = 4, family = skew_normal())

# Population-Level Effects: 
#                    Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
# Intercept                  -0.05      0.01    -0.06    -0.04       3335 1.00
# Groupwithout_depression     0.01      0.01    -0.01     0.02       2972 1.00
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
# sigma     0.05      0.00     0.04     0.05       2526 1.00


#======================================
# plot them!
#======================================
setwd(here("paper_figures"))
pdf("depression_diagnosis_behavioural.pdf", width = 10, height = 8)
grid.arrange(mean_rt_depression_plot, mean_acc_depression_plot, 
             rt_depression_plot, acc_depression_plot, ncol = 2)
dev.off()
setwd(here())


#======================================
# CDFS & CAFs
#======================================


# get a list of the ids with depression
with_dep <- all_data %>% 
  filter(id %in% depression_diagnosis_ids) %>% 
  pull(id)


# get a list of the ids without depression
without_dep <- all_data %>% 
  filter(id %!in% depression_diagnosis_ids) %>% 
  pull(id)


# cdfs
probs <- c(0.1, 0.3, 0.5, 0.7, 0.9)

cdf_with_con <- flanker_data %>% 
  filter(id %in% with_dep) %>% 
  filter(accuracy == 1, congruency == "congruent") %>% 
  rename(subject = id) %>% 
  flankr::cdf(.)

cdf_with_incon <- flanker_data %>% 
  filter(id %in% with_dep) %>% 
  filter(accuracy == 1, congruency == "incongruent") %>% 
  rename(subject = id) %>% 
  flankr::cdf(.)

cdf_without_con <- flanker_data %>% 
  filter(id %in% without_dep) %>% 
  filter(accuracy == 1, congruency == "congruent") %>% 
  rename(subject = id) %>% 
  flankr::cdf(.)

cdf_without_incon <- flanker_data %>% 
  filter(id %in% without_dep) %>% 
  filter(accuracy == 1, congruency == "incongruent") %>% 
  rename(subject = id) %>% 
  flankr::cdf(.)

# cafs
caf_with_con <- flanker_data %>% 
  filter(id %in% with_dep) %>% 
  filter(congruency == "congruent") %>% 
  rename(subject = id) %>% 
  flankr::caf(.)

caf_with_incon <- flanker_data %>% 
  filter(id %in% with_dep) %>% 
  filter(congruency == "incongruent") %>% 
  rename(subject = id) %>% 
  flankr::caf(.)

caf_without_con <- flanker_data %>% 
  filter(id %in% without_dep) %>% 
  filter(congruency == "congruent") %>% 
  rename(subject = id) %>% 
  flankr::caf(.)

caf_without_incon <- flanker_data %>% 
  filter(id %in% without_dep) %>% 
  filter(congruency == "incongruent") %>% 
  rename(subject = id) %>% 
  flankr::caf(.)



## plotting
setwd(here("paper_figures"))

png("depression_cdf_caf_plots.png", width = 8, height = 8, units = "in", res = 1200)
par(mfrow = c(1, 2))
do_clean_graphs()

# plot the cdfs
min_rt = min(cdf_with_con, cdf_with_incon, 
             cdf_without_con, cdf_without_incon)

max_rt = max(cdf_with_con, cdf_with_incon, 
             cdf_without_con, cdf_without_incon)


plot(cdf_with_con, probs, type = "b", xlim = c(min_rt - 50, max_rt + 50), 
     ylim = c(0, 1), pch = 19, col = "#E69F00", 
     xlab = "Response Time (ms)", ylab = "Cumulative Probability", cex = 1.5, 
     lwd = 2)
lines(cdf_with_incon, probs, type = "b", pch = 1, lty = 2, col = "#E69F00", 
      cex = 1.5, lwd = 2)
lines(cdf_without_con, probs, type = "b", pch = 15, col = "#56B4E9", 
      cex = 1.5, lwd = 2)
lines(cdf_without_incon, probs, type = "b", pch = 0, lty = 2, col = "#56B4E9", 
      cex = 1.5, lwd = 2)

# plot the cafs
min_rt = min(caf_with_con[1, ], caf_with_incon[1, ], 
             caf_without_con[1, ], caf_without_incon[1, ])
max_rt = max(caf_with_con[1, ], caf_with_incon[1, ], 
             caf_without_con[1, ], caf_without_incon[1, ])

plot(caf_with_con[1, ], caf_with_con[2, ], type = "b", 
     xlim = c(min_rt - 50, max_rt + 50), ylim = c(.7, 1), pch = 19, col = "#E69F00", 
     xlab = "Response Time (ms)", ylab = "Proportion Error", cex = 1.5, lwd = 2)
lines(caf_with_incon[1, ], caf_with_incon[2, ], type = "b", 
      pch = 1, lty = 2, col = "#E69F00", cex = 1.5, lwd = 2)
lines(caf_without_con[1, ], caf_without_con[2, ], type = "b", 
      pch = 15, col = "#56B4E9", cex = 1.5, lwd = 2)
lines(caf_without_incon[1, ], caf_without_incon[2, ], type = "b", 
      pch = 0, lty = 2, col = "#56B4E9", cex = 1.5, lwd = 2)

legend("bottomright", c("W Depression, Congruent","W Depression, Incongruent", 
                        "W/O Depression, Congruent", "W/O Depression, Incongruent"), 
       pch = c(19, 1, 15, 0),
       lty = c(1, 2, 1, 2), 
       lwd = 1.5,
       col = c("#E69F00", "#E69F00", "#56B4E9", "#56B4E9"), bty="n");
dev.off()



#======================================
# DSTP parameters
#======================================
setwd(here("model_data"))
all_data_dstp = read.csv("all_data_dstp.csv")
all_data_ssp = read.csv("all_data_ssp.csv")

# sort the participants
with_dep_dstp <- all_data_dstp %>% 
  filter(id %in% with_dep) %>% 
  mutate(Group = "With")
without_dep_dstp <- all_data_dstp %>% 
  filter(id %in% without_dep) %>% 
  mutate(Group = "Without")

depression_dstp <- rbind(with_dep_dstp, without_dep_dstp)
depression_dstp$Group <- as.factor(depression_dstp$Group)


# A
depression_a <- brm(A ~ Group, data = depression_dstp, cores = 4, 
              family = skew_normal())
print(summary(depression_a), digits = 4)

# C
depression_c<- brm(C ~ Group, data = depression_dstp, cores = 4, 
             family = skew_normal())
print(summary(depression_c), digits = 4)

# muT
depression_muT <- brm(muT ~ Group, data = depression_dstp, cores = 4, 
                family = skew_normal())
print(summary(depression_muT), digits = 4)

# muFl
depression_muFl <- brm(muFl ~ Group, data = depression_dstp, cores = 4, 
                 family = skew_normal())
print(summary(depression_muFl), digits = 4)

# muSS
depression_muSS <- brm(muSS ~ Group, data = depression_dstp, cores = 4, 
                 family = skew_normal())
print(summary(depression_muSS), digits = 4)

# muRS2
depression_muRS2 <- brm(muRS2 ~ Group, data = depression_dstp, cores = 4, 
                  family = skew_normal())
print(summary(depression_muRS2), digits = 4)

# ter
depression_ter <- brm(ter ~ Group, data = depression_dstp, cores = 4, 
                family = skew_normal())
print(summary(depression_ter), digits = 4)


## do all of the density plots

# A
a_plot <- depression_dstp %>% 
  ggplot(aes(x = A, colour = Group)) + 
  geom_density(aes(colour = Group, fill = Group), 
               alpha = 0.8) + 
  scale_colour_manual(values = c("#3182bd", "#9ecae1")) +
  scale_fill_manual(values = c("#3182bd", "#9ecae1")) +
  labs(x = "Parameter A") 

# C
c_plot <- depression_dstp %>% 
  ggplot(aes(x = C, colour = Group)) + 
  geom_density(aes(colour = Group, fill = Group), 
               alpha = 0.8) + 
  scale_colour_manual(values = c("#3182bd", "#9ecae1")) +
  scale_fill_manual(values = c("#3182bd", "#9ecae1")) +
  labs(x = "Parameter C") 

# mu_t
muT_plot <- depression_dstp %>% 
  ggplot(aes(x = muT, colour = Group)) + 
  geom_density(aes(colour = Group, fill = Group), 
               alpha = 0.8) + 
  scale_colour_manual(values = c("#3182bd", "#9ecae1")) +
  scale_fill_manual(values = c("#3182bd", "#9ecae1")) +
  labs(x = "Parameter muT") 

# mu_fl
muFl_plot <- depression_dstp %>% 
  ggplot(aes(x = muFl, colour = Group)) + 
  geom_density(aes(colour = Group, fill = Group), 
               alpha = 0.8) + 
  scale_colour_manual(values = c("#3182bd", "#9ecae1")) +
  scale_fill_manual(values = c("#3182bd", "#9ecae1")) +
  labs(x = "Parameter muFl") 

# mu_SS
muSS_plot <- depression_dstp %>% 
  ggplot(aes(x = muSS, colour = Group)) + 
  geom_density(aes(colour = Group, fill = Group), 
               alpha = 0.8) + 
  scale_colour_manual(values = c("#3182bd", "#9ecae1")) +
  scale_fill_manual(values = c("#3182bd", "#9ecae1")) +
  labs(x = "Parameter muSS") 

# mu_RS2
muRS2_plot <- depression_dstp %>% 
  ggplot(aes(x = muRS2, colour = Group)) + 
  geom_density(aes(colour = Group, fill = Group), 
               alpha = 0.8) + 
  scale_colour_manual(values = c("#3182bd", "#9ecae1")) +
  scale_fill_manual(values = c("#3182bd", "#9ecae1")) +
  labs(x = "Parameter muRS2") 

# ter
ter_plot <- depression_dstp %>% 
  ggplot(aes(x = ter, colour = Group)) + 
  geom_density(aes(colour = Group, fill = Group), 
               alpha = 0.8) + 
  scale_colour_manual(values = c("#3182bd", "#9ecae1")) +
  scale_fill_manual(values = c("#3182bd", "#9ecae1")) +
  labs(x = "Parameter ter") 


setwd(here("paper_figures"))
pdf("depression_dstp_plots.pdf", width = 10, height = 8)
grid.arrange(a_plot, c_plot, muT_plot,
             muFl_plot, muSS_plot, muRS2_plot, 
             ter_plot, ncol = 3)
dev.off()
setwd(here())




#======================================
# SSP parameters
#======================================
setwd(here("model_data"))
all_data_ssp = read.csv("all_data_ssp.csv")

# sort the participants
with_dep_ssp <- all_data_ssp %>% 
  filter(id %in% with_dep) %>% 
  mutate(Group = "With")
without_dep_ssp <- all_data_ssp %>% 
  filter(id %in% without_dep) %>% 
  mutate(Group = "Without")

depression_ssp<- rbind(with_dep_ssp, without_dep_ssp)
depression_ssp$Group <- as.factor(depression_ssp$Group)


# a
depression_a <- brm(a ~ Group, data = depression_ssp, cores = 4, 
                    family = skew_normal())
print(summary(depression_a), digits = 4)

# ter
depression_ter <- brm(ter ~ Group, data = depression_ssp, cores = 4, 
                      family = skew_normal())
print(summary(depression_ter), digits = 4)

# p
depression_p <- brm(p ~ Group, data = depression_ssp, cores = 4, 
                    family = skew_normal())
print(summary(depression_p), digits = 4)

# rd
depression_rd <- brm(rd ~ Group, data = depression_ssp, cores = 4, 
                    family = skew_normal())
print(summary(depression_rd), digits = 4)


# sda
depression_sda <- brm(sda ~ Group, data = depression_ssp, cores = 4, 
                    family = skew_normal())
print(summary(depression_sda), digits = 4)


## do all of the density plots

# a
a_plot <- depression_ssp %>% 
  ggplot(aes(x = a, colour = Group)) + 
  geom_density(aes(colour = Group, fill = Group), 
               alpha = 0.8) + 
  scale_colour_manual(values = c("#3182bd", "#9ecae1")) +
  scale_fill_manual(values = c("#3182bd", "#9ecae1")) +
  labs(x = "Parameter a") 

# ter
ter_plot <- depression_ssp %>% 
  ggplot(aes(x = ter, colour = Group)) + 
  geom_density(aes(colour = Group, fill = Group), 
               alpha = 0.8) + 
  scale_colour_manual(values = c("#3182bd", "#9ecae1")) +
  scale_fill_manual(values = c("#3182bd", "#9ecae1")) +
  labs(x = "Parameter ter") 

# p
p_plot <- depression_ssp %>% 
  ggplot(aes(x = p, colour = Group)) + 
  geom_density(aes(colour = Group, fill = Group), 
               alpha = 0.8) + 
  scale_colour_manual(values = c("#3182bd", "#9ecae1")) +
  scale_fill_manual(values = c("#3182bd", "#9ecae1")) +
  labs(x = "Parameter p") 

# rd
rd_plot <- depression_ssp %>% 
  ggplot(aes(x = rd, colour = Group)) + 
  geom_density(aes(colour = Group, fill = Group), 
               alpha = 0.8) + 
  scale_colour_manual(values = c("#3182bd", "#9ecae1")) +
  scale_fill_manual(values = c("#3182bd", "#9ecae1")) +
  labs(x = "Parameter rd") 

# sda
sda_plot <- depression_ssp %>% 
  ggplot(aes(x = sda, colour = Group)) + 
  geom_density(aes(colour = Group, fill = Group), 
               alpha = 0.8) + 
  scale_colour_manual(values = c("#3182bd", "#9ecae1")) +
  scale_fill_manual(values = c("#3182bd", "#9ecae1")) +
  labs(x = "Parameter sda") 


setwd(here("paper_figures"))
pdf("depression_ssp_plots.pdf", width = 10, height = 8)
grid.arrange(a_plot, ter_plot, p_plot,
             rd_plot, sda_plot, ncol = 3)
dev.off()
setwd(here())

