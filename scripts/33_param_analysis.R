library("here")
library("tidyverse")
library("lookout")
library(brms)

params <- readRDS(here("data", "processed", "dstp_params", "dstp_params.Rds"))


dat <- params
dat$mu_ss <- ifelse(dat$mu_ss > 0.7, NA, dat$mu_ss)
dat$C <- ifelse(dat$C > 0.4, NA, dat$C)


imp <- mice::mice(
  dat, 
  maxit = 1
)

dat2 <- mice::complete(imp, 1)

# Deviation Coding. This coding system compares the mean of the dependent 
# variable for a given level to the overall mean of the dependent variable. 
dat2$exp <- factor(dat2$exp)
contrasts(dat2$exp) <- contr.sum(2)
contrasts(dat2$exp) 

lo <- lookout(dat2[, 6:7])
lo
autoplot(lo)



bf_mu_ta <- bf(mu_ta ~ exp) + skew_normal()
bf_mu_fl <- bf(mu_fl ~ exp) + skew_normal()
bf_mu_ss <- bf(mu_ss ~ exp) + skew_normal()
bf_mu_rs2 <- bf(mu_rs2 ~ exp) + skew_normal()
bf_ter <- bf(ter ~ exp) + skew_normal()


fit <- brm(
  bf_mu_ta + bf_mu_fl + bf_mu_ss + bf_mu_rs2 + bf_ter + 
    set_rescor(FALSE), 
  data = dat2, 
  chains = 4, 
  cores = 4,
  control = list(adapt_delta = 0.95),
  backend = "cmdstanr"
)

summary(fit)

bayes_R2(fit)


pp_check(fit, resp = "muss")

conditional_effects(fit, "exp", resp = "muss")





# Param A ----

# Height of response selection boundary

hist(dat2$A)


mod_A <- brm(
  bf(A ~ exp),
  family = skew_normal(),
  data = dat2, 
  cores = 4
)

summary(mod_A)
pp_check(mod_A) + xlim(c(0, 0.75))
bayes_R2(mod_A)


# Param C ----
# Height of stimulus selection boundary

hist(dat2$C)

foo <- dat2
foo$C2 <- foo$C * 10

mod_C <- brm(
  bf(C ~ exp),
  family = skew_normal(),
  data = foo, 
  cores = 4
)

summary(mod_C)
pp_check(mod_C) 
bayes_R2(mod_C)
conditional_effects(mod_C, "exp")

# Effect size.
draws <- as.data.frame(mod_C)
draws$cohens_d <- draws$b_exp2 / draws$sigma

hist(draws$cohens_d)
mean(draws$cohens_d)
# [1] -0.1288013
hdi(draws$cohens_d, credMass=0.95)
# lower       upper 
# -0.24280348 -0.02137142 



# Param mu_ta ----
# Drift rate for central target during response selection stage 1

hist(dat2$mu_ta)

dat2 %>% 
  ggplot(aes(x = exp, y = mu_fl, fill = exp)) +
  geom_violin() +
  geom_boxplot(width=.1, outlier.colour=NA) 


prior <- get_prior(
  mu_ta ~ exp
)

mod_muta <- brm(
  bf(mu_ta ~ exp),
  family = asym_laplace(),
  data = dat2, 
  cores = 4
)

summary(mod_muta)
pp_check(mod_muta) + xlim(c(0, 0.75))
bayes_R2(mod_muta)

conditional_effects(mod_muta, "exp")


# Param mu_ss ----

# Drift rate for stimulus selection

contrasts(dat2$exp) <- contr.treatment(2)

hist(dat2$mu_ss)

prior <- get_prior(
  mu_ss ~ 0 + Intercept + exp,
  data = dat2
)

prior_muss <- c(
  prior(normal(0, 0.05), class = b, coef = exp2)
)

mod_muss <- brm(
  bf(mu_ss ~ exp),
  family = asym_laplace(),
  data = dat2, 
  prior = prior_muss,
  cores = 4,
  sample_prior = TRUE
)
summary(mod_muss)

plot(hypothesis(mod_muss, "exp2 < 0"))

pp_check(
  mod_muss, heigth, type = "error_scatter_avg_vs_x", x = dat2$mu_ss) +
  geom_hline(yintercept = 0, color = "red", size = 1.5) +
  scale_x_continuous(breaks = 1:10)

post <- posterior_samples(mod_muss)

post %>% 
  ggplot(aes(x = b_Intercept, y = 0)) +
  geom_halfeyeh() +
  scale_y_continuous(NULL, breaks = NULL)

post %>% 
  ggplot(aes(x = b_Intercept, y = prior)) 

summary(mod_muss)
pp_check(mod_muss) + xlim(c(0, 0.75))
bayes_R2(mod_muss)

conditional_effects(mod_muss, "exp")
# plot(conditional_effects(mod_muss, "exp"), points = TRUE)

# Effect size.
draws <- as.data.frame(mod_muss)
draws$cohens_d <- draws$b_exp2 / draws$sigma

hist(draws$cohens_d)
mean(draws$cohens_d)

hdi(draws$cohens_d, credMass=0.95)
# lower      upper 
# -0.7914693 -0.2815002 

print(prior_summary(mod_muss, all = FALSE), show_df = FALSE)


loo(mod_muss)
loo1 <- loo(mod_muss, save_psis = TRUE)

yrep1 <- posterior_predict(mod_muss)

library(tidyverse)
library(rstan)
library(brms)
library(bayesplot)
library(loo)
library(tidybayes)


ppc_loo_pit_overlay(
  yrep = yrep1, 
  y = dat2$mu_ss, 
  lw = weights(loo1$psis_object)
) + 
  ggtitle("LOO-PIT Model")

pp_check(mod_muss, type = "stat", stat = 'median', nsamples = 100)



#--------------------------------------------------------------------
# Param mu_fl 
# Drift rate for flankers during response selection stage 1
#--------------------------------------------------------------------

hist(dat2$mu_fl)

dat2 %>% 
  ggplot(aes(x = exp, y = mu_fl, fill = exp)) +
  geom_violin() +
  geom_boxplot(width=.1, outlier.colour=NA) 

prior <- get_prior(
  mu_fl ~ 0 + Intercept + exp,
  data = dat2
)

prior_mufl <- c(
  prior(normal(0, 0.05), class = b, coef = exp2)
)


mod_mufl <- brm(
  bf(mu_fl ~ exp),
  family = skew_normal(),
  data = dat2, 
  cores = 4
)

summary(mod_mufl)
pp_check(mod_mufl) + xlim(c(0, 0.75))
bayes_R2(mod_mufl)


#--------------------------------------------------------------------
# Param mu_rs2 
# Drift rate for stage 2 of response selection
#--------------------------------------------------------------------

hist(dat2$mu_rs2)

dat2 %>% 
  ggplot(aes(x = exp, y = mu_rs2, fill = exp)) +
  geom_violin() +
  geom_boxplot(width=.1, outlier.colour=NA) 

prior <- get_prior(
  mu_rs2 ~ 0 + Intercept + exp,
  data = dat2
)

prior_mufl <- c(
  prior(normal(0, 0.05), class = b, coef = exp2)
)

mod_murs2 <- brm(
  bf(mu_rs2 ~ exp),
  family = skew_normal(),
  data = dat2, 
  cores = 4
)

summary(mod_murs2)
pp_check(mod_murs2) #+ xlim(c(0, 0.75))
bayes_R2(mod_murs2)

#--------------------------------------------------------------------
# Param ter
# Non-decision time
#--------------------------------------------------------------------

hist(dat2$ter)

dat2 %>% 
  ggplot(aes(x = exp, y = ter, fill = exp)) +
  geom_violin() +
  geom_boxplot(width=.1, outlier.colour=NA) 

prior <- get_prior(
  ter ~ 0 + Intercept + exp,
  data = dat2
)

prior_ter <- c(
  prior(normal(0, 0.05), class = b, coef = exp2),
  prior(student_t(3, 0, 1), class = sigma),
  prior(gamma(2, 0.1), class = nu)
)

mod_ter <- brm(
  bf(ter ~ exp),
  family = student(),
  data = dat2, 
  cores = 4,
  prior = prior_ter,
  sample_prior = TRUE
)

plot(hypothesis(mod_ter, "exp2 < 0"))
print(prior_summary(mod_ter, all = FALSE), show_df = FALSE)


summary(mod_ter)
pp_check(mod_ter) #+ xlim(c(0, 0.75))
bayes_R2(mod_ter)


post <- posterior_samples(mod_muss)

post %>% 
  ggplot(aes(x = b_exp2, y = 0)) +
  geom_halfeyeh() +
  scale_y_continuous(NULL, breaks = NULL)


conditional_effects(mod_ter, "exp")
# plot(conditional_effects(mod_muss, "exp"), points = TRUE)

# Effect size.
draws <- as.data.frame(mod_ter)
draws$cohens_d <- draws$b_exp2 / draws$sigma

hist(draws$cohens_d)
mean(draws$cohens_d)

HDInterval::hdi(draws$cohens_d, credMass=0.95)
# lower      upper 
# -0.7914693 -0.2815002 

loo(mod_ter)


















hist(params_surprise$mu_ss)

params %>% 
  ggplot(aes(x=C, group=exp, fill=exp)) +
  geom_density(alpha=0.5)


aa <- params


t.test(pamams$mu_ta ~ pamams$exp)


tot_df <- rbind(control_df, surpr_exp_df, nosurpr_exp_df)
tot_df$experiment <- c(rep("control", 10), rep("exp_surprise", 20), rep("exp_nosurprise", 20))

tot_df$exp <- ifelse(tot_df$experiment == "control", 0, 1)









fm <- lm(cbind(mu_ta, mu_fl, mu_ss) ~ exp, data = tot_df)
Anova(fm)
summary(fm)

fm <- lm(mu_fl ~ exp, data = tot_df)
Anova(fm)
summary(fm)

fm <- lm(mu_ss ~ exp, data = tot_df)
Anova(fm)
summary(fm)

fm <- lm(C ~ exp, data = tot_df)
Anova(fm)
summary(fm)


t.test(mu_fl ~ experiment, data = tot_df)
t.test(mu_ss ~ experiment, data = tot_df)
t.test(ter ~ experiment, data = tot_df)
t.test(A ~ experiment, data = tot_df)
t.test(C ~ experiment, data = tot_df)








# select the data of one subject
one_subj_dat <- dat %>% 
  dplyr::filter(
    subject == 2
  )






m1 <- fitDSTP(data = df_surprise)

m2 <- fitDSTP(data = df_nosurprise)

plot1 <- plotFitDSTP(modelFit = m1, data = df_surprise)
plot2 <- plotFitDSTP(modelFit = m2, data = df_nosurprise)



plot <- plotFitDSTP(modelFit = fit5, data = supr_df)




fit3 <- fitDSTP(data = dd, conditionName = "No Surprise")
# control condition (no surprise)
# (1) mu_ta (2) mu_fl (3) A (4) mu_ss (5) C (6) mu_rs2 (7) t_er
# mu = rate or drift of a given diffusion process and condition; 
# ta = target; 
# fl = flanker; 
# A = criterion for response selection; 
# SS = stimulus selection process; 
# C = criterion for stimulus selection; 
# RS2 = response selection process in Phase 2; 
# t_er = nondecisional parameter that represents time used for stimulus encoding, 
# response execution, and so forth.
# $bestParameters
# [1] 0.130 0.111 0.000 0.049 0.370 1.192 0.190
# 
# $g2
# [1] 46.63694
# 
# $bBIC
# [1] 1977.794

# surprise condition in the 2-video condition
# $bestParameters
# [1] 0.160 0.131 0.005 0.053 0.281 1.369 0.119
# 
# $g2
# [1] 24.97205
# 
# $bBIC
# [1] 1938.025

# no surprise condition in the 2-video condition
# $bestParameters
# [1] 0.215 0.142 0.096 0.131 0.300 1.373 0.096
# 
# $g2
# [1] 118.0751
# 
# $bBIC
# [1] 2045.489

# Height of response selection boundary
# Height of stimulus selection boundary
# Drift rate for central target during response selection stage 1 
# Drift rate for flankers during response selection stage 1 
# Drift rate for stimulus selection
# Drift rate for stage 2 of response selection
# Non-decision time


# $bestParameters
# [1] 0.111 0.083 0.016 0.073 0.338 1.272 0.261
#
# $g2
# [1] 39.5643
#
# $bBIC
# [1] 1959.925

# surprise NO
# $bestParameters
# [1] 0.133 0.109 0.005 0.065 0.322 1.367 0.218
#
# $g2
# [1] 48.37825
#
# $bBIC
# [1] 1985.64

plot <- plotFitDSTP(modelFit = fit2, data = dd)

#' Conclusion: The diffusion model is adequate for the control condition (no
#' surprise in any trial), but not for the experimental conditions, in which
#' in half of the trials the videos show a surprising event. The data will not
#' be analyzed with a DM approach, because the process is more complex than a
#' single decision task.


nBootstraps <- 1000

mean(tapply(mydata$rtTukey, mydata$subject_name, length))

nTrials <- 140

bestParameters <- c(0.130, 0.111, 0.000, 0.049, 0.370, 1.192, 0.190)

bootData <- matrix(0, nrow = nBootstraps, ncol = length(bestParameters))

for (i in 1:nBootstraps) {
  
  simData <- simulateDSTP(parms = bestParameters, nTrials = nTrials)
  
  fit <- fitDSTP(data = simData, 
                 parms = bestParameters,
                 nTrials = 10000,
                 multipleSubjects = FALSE)
  
  bootData[i, ] <- fit$bestParameters
}

apply(bootData, 2, sd) / sqrt(dim(bootData)[1])



library(BEST)

mean(dat2$mu_ss)
sd(dat2$mu_ss)


priors <- list(muM = 0.37, muSD = 0.8)
y1 <- dat2[dat2$exp == "control", ]$mu_ss
y2 <- dat2[dat2$exp == "surprise", ]$mu_ss

BESTout <- BESTmcmc(
  y1, y2, 
  priors=priors, 
  parallel=FALSE
)
plot(BESTout)
plotAll(BESTout)
summary(BESTout)


mean(params$mu_rs2)
sd(params$mu_rs2)

priors <- list(muM = 1.31, muSD = 0.12)
BESTout <- BESTmcmc(
  params_control$mu_rs2, params_surprise$mu_rs2, 
  priors=priors, parallel=FALSE
)
plot(BESTout)
plotAll(BESTout)


