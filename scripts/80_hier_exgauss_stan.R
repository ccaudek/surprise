#' Hierarchical Bayesian Ex-Gaussian Estimation for the Surprise Project
#'
#' This script was last modified on "Mon Jul  2 09:44:15 2018"
#' Corrado Caudek


# require("rstan")
# require("retimes")
# require("rethinking")
# 
# source("libsAndMore.R")
# source("DBDA2E-utilities.R")

if (!require("pacman")) install.packages("pacman")
library("pacman")
pacman::p_load(
  tidyverse, bayesplot, knitr, rmarkdown, forcats, rstan, lme4,
  car, effects, brms, retimes, gamlss.dist, shinystan
)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


# Stan model --------------------------------------------------------------


modelString <- "

data {
  int<lower=0> N;        // number of trials
  vector<lower=0>[N] rt; // vector of response times in seconds
  int<lower=0> J;        // number of subjects
  int<lower=0> id[N];    // vector of subjects' id (1,2,3,...)
}

parameters {
  real<lower=0> mu[J];
  real<lower=0> sigma[J];
  real<lower=0> lambda[J];
  real<lower=0> mu_m;
  real<lower=0> sigma_m;
  real<lower=0> mu_s;
  real<lower=0> sigma_s;
  real<lower=0> mu_l;
  real<lower=0> sigma_l;
}

transformed parameters {
  real<lower=0> tau;
  tau = 1 / mu_l; // transform hyperparameter mu_lambda
}

model {
  mu ~ normal(mu_m, sigma_m);
  sigma ~ normal(mu_s, sigma_s);
  lambda ~ normal(mu_l, sigma_l);
  for (n in 1:N) {
    rt[n] ~ exp_mod_normal(
    mu[id[n]], 
    sigma[id[n]], 
    lambda[id[n]]
    );
  }
}
"



# Read data ---------------------------------------------------------------


data_cntr <- read.csv("surprise_cntr.csv")
data_cntr$experiment <- "Control"
data_expr <- read.csv("surprise_expr.csv")
data_expr$experiment <- "Surprise"

surprise_df <- rbind(data_cntr, data_expr)


# Select subjects and trials ----------------------------------------------


selected_trials <- surprise_df %>% 
  dplyr::filter(
    experiment == "Control" & # "Control" "Surprise"
      is_congruent_trial == "Congruent" &
      correct == 1 &
      rtTukey < 1500 &
      trials_after_clip > 0 & trials_after_clip < 4
      #trials_after_clip == 0 
  )

# Remove rows with NA on the rt variable
thedat <- selected_trials[!is.na(selected_trials$rtTukey), ]

# subject's id from 1 to n_subjects
thedat$id <- as.numeric(factor(as.character(thedat$id)))
sort(unique(thedat$id))


# Data in the appropriate format for Stan
stan_data = list(
  rt = thedat$rtTukey / 1000, 
  N = length(thedat$rtTukey),
  J = length(unique(thedat$id)), 
  id = thedat$id
)


# Fit the model
fit_hier <- stan(
  model_code = modelString, 
  data = stan_data, 
  verbose = FALSE,
  chains = 4, 
  iter = 2000
)

# M <- 1e4
# fit_hier <- stan(model_code=modelString, data = stan_data,
#   chains = 4, iter = M / 2)

posterior <- rstan::extract(fit_hier)

print(
  fit_hier, 
  c("mu_m", "mu_s", "tau"), 
  probs = c(0.025, 0.5, 0.975), 
  digits = 3
)


#' # Without the trial immediately following the video clip

# Control condition -------------------------------------------------------

# Congruent
#      mean se_mean   sd 2.5%  50% 97.5% n_eff Rhat
# mu_m 0.46       0 0.01 0.43 0.46  0.49  4000    1
# mu_s 0.05       0 0.00 0.05 0.05  0.06  4000    1
# tau  0.08       0 0.01 0.07 0.08  0.10  4000    1

# Incongruent
#      mean se_mean   sd 2.5%  50% 97.5% n_eff Rhat
# mu_m 0.51       0 0.01 0.50 0.51  0.53  4000    1
# mu_s 0.06       0 0.00 0.06 0.06  0.06  3512    1
# tau  0.05       0 0.00 0.04 0.05  0.06  3377    1


# Surprise condition ------------------------------------------------------


# Congruent
#      mean se_mean   sd 2.5%  50% 97.5% n_eff Rhat
# mu_m 0.55       0 0.02 0.52 0.55  0.58  4000    1
# mu_s 0.09       0 0.00 0.08 0.09  0.09  4000    1
# tau  0.12       0 0.00 0.11 0.12  0.13  4000    1


# Incongruent
#      mean se_mean   sd 2.5%  50% 97.5% n_eff Rhat
# mu_m 0.57       0 0.01 0.55 0.57  0.59  4000    1
# mu_s 0.07       0 0.00 0.07 0.07  0.08  4000    1
# tau  0.10       0 0.00 0.09 0.10  0.10  4000    1



#' # Only the trial immediately following the video clip


# Control condition -------------------------------------------------------

# Congruent
#      mean se_mean   sd 2.5%  50% 97.5% n_eff Rhat
# mu_m 0.58       0 0.02 0.55 0.58  0.61  4000    1
# mu_s 0.07       0 0.01 0.06 0.07  0.08  1900    1
# tau  0.08       0 0.01 0.07 0.08  0.10  1643    1

# Incongruent
#      mean se_mean   sd 2.5%  50% 97.5% n_eff Rhat
# mu_m 0.60       0 0.01 0.57 0.60  0.62  4000    1
# mu_s 0.06       0 0.00 0.05 0.06  0.06  1343    1
# tau  0.07       0 0.01 0.05 0.07  0.08  1368    1



# Surprise condition ------------------------------------------------------

# Congruent
#      mean se_mean   sd 2.5%  50% 97.5% n_eff Rhat
# mu_m 0.60       0 0.01 0.57 0.60  0.63  4000    1
# mu_s 0.09       0 0.00 0.08 0.09  0.10  1989    1
# tau  0.11       0 0.01 0.10 0.11  0.12  1891    1

# Incongruent
#      mean se_mean   sd 2.5%  50% 97.5% n_eff Rhat
# mu_m 0.61       0 0.01 0.59 0.61  0.64  4000    1
# mu_s 0.08       0 0.00 0.08 0.08  0.09  2306    1
# tau  0.09       0 0.00 0.08 0.09  0.10  1894    1


print(
  fit_hier, c("mu"), digits = 3
)

print(
  fit_hier, c("sigma"), digits = 3
)

print(
  fit_hier, c("lambda"), digits = 3
)



# Evaluate goodness of fit ------------------------------------------------

thedat$rt <- thedat$rtTukey / 1000

# number of subjects
nsub <- length(unique(thedat$id))
nsub

# recover individual parameters
post_est <- summary(fit_hier)

mu <- post_est$summary[1:nsub, 1]
sigma <- post_est$summary[(nsub+1):(2*nsub), 1]
tau <- 1 / post_est$summary[(2*nsub+1):(3*nsub), 1]


subject_index <- 43

onesubj <- thedat[thedat$id == subject_index, ]
dim(onesubj)

xmax <- 1.5
hist(onesubj$rtTukey/1000, freq=FALSE, breaks=20, xlim=c(0.2, xmax))

x <- seq(0.2, xmax, length.out = 200)
lines(x, 
      dexgauss(x, mu[subject_index], sigma[subject_index], tau[subject_index]), 
      lty=1, 
      col="red", 
      lwd=2
)


prob_list <- seq(0.05, 0.95, length.out = 5)
prob_list

# theoretical quantiles
q_t <- qexGAUS(prob_list, mu[subject_index], sigma[subject_index], tau[subject_index])

# empirical quantiles
q_e <- quantile(onesubj$rt, probs = prob_list)

plot(q_t, q_e)
abline(0, 1)
cor(q_t, q_e)



# Loop for all estimates --------------------------------------------------


tq_list <- list()
eq_list <- list()
for (subject_index in 1:nsub) {
  onesubj <- thedat[thedat$id == subject_index, ]
  # theoretical quantiles
  tq_list[[subject_index]] <- qexGAUS(prob_list, 
                                      mu[subject_index], 
                                      sigma[subject_index], 
                                      tau[subject_index]
                                     )
  # empirical quantiles
  eq_list[[subject_index]] <- quantile(onesubj$rt, probs = prob_list)
  
}



tq_df <- as.data.frame(do.call(rbind, tq_list))
tq_df$subject <- 1:nsub
eq_df <- as.data.frame(do.call(rbind, eq_list))
eq_df$subject <- 1:nsub
new_names <- c("V1", "V2", "V3", "V4", "V5", "subject")
names(eq_df) <- new_names

long_t <- gather(tq_df, q, rt, V1:V5, factor_key=TRUE)
long_t$quantile <- "theoretical"

long_e <- gather(eq_df, q, rt, V1:V5, factor_key=TRUE)
long_e$quantile <- "empirical"

q_df <- inner_join(long_e, long_t, by = c("subject", "q"))

new_df <- data.frame(
  subject = q_df$subject,
  y = q_df$rt.x,
  x = q_df$rt.y,
  q = q_df$q
)

# scale for each subject
new_df$zy <- with(new_df, ave(y, subject, FUN = scale))
new_df$zx <- with(new_df, ave(x, subject, FUN = scale))

#' A simple Pearson's correlation coefficient method is usually employed when there 
#' are complete independent data points for both outcome variables. However, researchers 
#' often deal with correlated observations in a longitudinal setting with missing values 
#' where a simple Pearson's correlation coefficient method cannot be used. A random 
#' regression mixed model was used to estimate correlation coefficients in the present 
#' within-subjects data set.

fit <- brm(zy ~ zx + (1 + zx | subject), data = new_df)

print(fit, digits = 4)
# no intercept model
#    Estimate Est.Error    l-95% CI u-95% CI Eff.Sample   Rhat
# zx   0.9958    0.0066   0.9828   1.0085          4000 0.9997

# model with intercept 
#           Estimate Est.Error l-95% CI u-95% CI Eff.Sample   Rhat
# Intercept  -0.0001    0.0057  -0.0113   0.0112       4000 0.9995
# zx          0.9958    0.0065   0.9832   1.0090       4000 0.9992


formatC(1000*colMeans(tq_df[1:5]), digits = 1, format = "f")
# V1       V2       V3       V4       V5 
# "0.4090" "0.4827" "0.5333" "0.5972" "0.7627" 
formatC(1000*colMeans(eq_df[1:5]), digits = 4, format = "f")
# V1       V2       V3       V4       V5 
# "0.4113" "0.4831" "0.5380" "0.6043" "0.7642" 

discrepancy <- 1000*tq_df[, 1:5] - 1000*eq_df[, 1:5]

formatC(colMeans(discrepancy), digits = 4, format = "f")
# V1        V2        V3        V4        V5 
# "-0.0023" "-0.0005" "-0.0048" "-0.0071" "-0.0015" 
# standard deviation
formatC(sqrt(diag(cov(discrepancy))), digits = 4, format = "f")
# V1       V2       V3       V4       V5 
# "0.0088" "0.0080" "0.0097" "0.0124" "0.0281" 










# Compute distribution function for empirical quantiles
p <- NULL
for (i in 1:5) {
  p[i] <- pexGAUS(q[i], param)
}

# Compute theoretical quantiles
qexGAUS(p, param[1], param[2], param[3])












launch_shinystan(fit_hier)





# -------------------------------------------------------------------
# plot Hierarchical estimates

par(mfrow=c(1, 3))

# plot mu
foo <- summary(fit_hier)$summary

plot( mu , foo[1:n_subj , 1] ,
  xlim = c(0, 1.3), ylim = c(0, 1.3),
  ylab=expression(hat(mu)), xlab=expression(mu) )
abline(lm(foo[1:n_subj , 1] ~ mu), col="red")
abline( 0 , 1 , lty=2 )
cor( mu , foo[1:n_subj , 1] )

# plot sigma
plot( sigma , foo[(n_subj+1):(2*n_subj) , 1]  ,
  xlim = c(0, 1.3), ylim = c(0, 1.3),
  ylab=expression(hat(sigma)), xlab=expression(sigma) )
abline(lm(foo[(n_subj+1):(2*n_subj) , 1] ~ sigma), col="red")
abline( 0 , 1 , lty=2)

# plot tau
plot( tau , 1/foo[(2*n_subj+1):(3*n_subj) , 1] ,
  xlim = c(0, 1.3), ylim = c(0, 1.3),
  ylab=expression(hat(tau)), xlab=expression(tau) )
abline(lm(1/foo[(2*n_subj+1):(3*n_subj) , 1] ~ tau), col="red")
abline( 0 , 1 , lty=2 )














# # -------------------------------------------------------------------
# # Simulate data
# 
# set.seed(123)
# 
# n_subj <- 30
# n_rep  <- 50
# mu    <- rnorm(n_subj, 0.95, 0.10)
# sigma <- rnorm(n_subj, 0.45, 0.10)
# tau   <- rnorm(n_subj, 0.30, 0.10)
# rt <- rexgauss(n_subj * n_rep,
#   mu = mu, sigma = sigma, tau = tau, positive = TRUE)
# id <- rep(1:n_subj, n_rep)
# N <- n_subj * n_rep
# J <- n_subj
# mydata <- data.frame(id, rt)


#--------------------------------------------------------------------
# Fit model

# stan_data = list(
#   rt = mydata$rt, 
#   N = length(mydata$rt),
#   J = length(unique(mydata$id)), 
#   id = mydata$id
# )


#--------------------------------------------------------------------
# Recover parameters from pooled ML analysis

MLest <- timefit(mydata$rt)
MLest

mu_ml    <- rep(NA, n_subj)
sigma_ml <- rep(NA, n_subj)
tau_ml   <- rep(NA, n_subj)

for (i in 1:n_subj) {
  temp <- subset(mydata, id==i)
  out <- timefit(temp$rt)
  mu_ml[i] <- out@par[1]
  sigma_ml[i] <- out@par[2]
  tau_ml[i] <- out@par[3]
  rm(temp, out)
}

# -------------------------------------------------------------------
# plot ML estimates

par(mfrow=c(1,3))

# plot mu

plot( mu , mu_ml ,
  xlim = c(0, 1.3), ylim = c(0, 1.3),
  ylab=expression(hat(mu)), xlab=expression(mu) )
abline(lm(mu_ml ~ mu), col="red")
abline( 0 , 1 , lty=2 )
cor( mu , mu_ml )

# plot sigma
plot( sigma , sigma_ml  ,
  xlim = c(0, 1.3), ylim = c(0, 1.3),
  ylab=expression(hat(sigma)), xlab=expression(sigma) )
abline(lm(sigma_ml ~ sigma), col="red")
abline( 0 , 1 , lty=2)

# plot tau
plot( tau , tau_ml ,
  xlim = c(0, 1.3), ylim = c(0, 1.3),
  ylab=expression(hat(tau)), xlab=expression(tau) )
abline(lm(tau_ml ~ tau), col="red")
abline( 0 , 1 , lty=2 )


#--------------------------------------------------------------------
# Plot posteriors
# True values are indicated by a dashed red line

posteriors = extract(fit_hier, permuted = TRUE)
par(mfrow = c(1,3))
hist(posteriors$mu,
  ylab = '',
  xlab = 'Mu',
  main = 'Posterior Distribution for Mu')
abline(v = mu, col = col.alpha("black",0.2), lty = 2, lw = 2)

hist(posteriors$sigma,
  ylab = '',
  xlab = 'Sigma',
  main = 'Posterior Distribution for Sigma')
abline(v = sigma, col = col.alpha("black",0.2), lty = 2, lw = 2)

hist(posteriors$tau,
  ylab = '',
  xlab = 'Tau',
  main = 'Posterior Distribution for Tau')
abline(v = tau, col = col.alpha("black",0.2), lty = 2, lw = 2)


## TODO: add the correlation value on the title of each graph

## TODO: repeat this simulation n times and show the distribution
## of the hyperparameters (mu, sigma, tau) estimated from the
## hierarchical approach and those estimated from the complete
## pooling from the ML approach

#--------------------------------------------------------------------
#--------------------------------------------------------------------
#--------------------------------------------------------------------

codaSamples <- as.mcmc.list(posteriors)
summary(codaSamples)
