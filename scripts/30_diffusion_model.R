require("rtdists")



# Exp. 1; Wagenmakers, Ratcliff, Gomez, & McKoon (2008, JML)
data(speed_acc)

# remove excluded trials:
speed_acc <- droplevels(speed_acc[!speed_acc$censor, ])

# create numeric response variable where 1 is an error and 2 a correct response:
speed_acc$corr <- with(speed_acc, as.numeric(stim_cat == response))+1

# select data from participant 11, accuracy condition, non-word trials only
p11 <- speed_acc[speed_acc$id == 11 &
                   speed_acc$condition == "accuracy" &
                   speed_acc$stim_cat == "nonword",]
prop.table(table(p11$corr))
#          1          2
# 0.04166667 0.95833333


# surprise exp, surprise trials
p11 <- subset(all_trials, subject_name == "03BM23_F20" & is_surprise_clip == 0)
p11$corr <- p11$correct + 1
prop.table(table(p11$corr))
p11$response <- factor(p11$response)


ll_lba <- function(pars, rt, response) {
  d <- dLBA(rt = rt, response = response,
            A = pars["A"],
            b = pars["A"]+pars["b"],
            t0 = pars["t0"],
            mean_v = pars[c("v1", "v2")],
            sd_v = c(1, pars["sv"]),
            silent=TRUE)
  if (any(d == 0)) return(1e6)
  else return(-sum(log(d)))
}

start <- c(runif(3, 0.5, 3), runif(2, 0, 0.2), runif(1))

names(start) <- c("A", "v1", "v2", "b", "t0", "sv")

p11_norm <- nlminb(start,  ll_lba,  lower = c(0,  -Inf,  0,  0,  0,  0),
                   rt = p11$rt,  response = p11$corr)

p11_norm[1:3]
# $par
#          A         v1         v2          b         t0         sv
#  0.1182940 -2.7409230  1.0449963  0.4513604  0.1243441  0.2609968
#
# $objective
# [1] -211.4202
#
# $convergence
# [1] 0


ll_diffusion <- function(pars, rt, response)
{
  densities <- ddiffusion(rt, response=response,
                          a=pars["a"],
                          v=pars["v"],
                          t0=pars["t0"],
                          sz=pars["sz"],
                          st0=pars["st0"],
                          sv=pars["sv"])
  if (any(densities == 0)) return(1e6)
  return(-sum(log(densities)))
}

start <- c(runif(2, 0.5, 3), 0.1, runif(3, 0, 0.5))
names(start) <- c("a", "v", "t0", "sz", "st0", "sv")
p11_diff <- nlminb(start, ll_diffusion, lower = 0,
                   rt=p11$rt, response=p11$corr)
p11_diff[1:3]
# $par
#         a         v        t0        sz       st0        sv
# 1.3206011 3.2727202 0.3385602 0.4621645 0.2017950 1.0551706
#
# $objective
# [1] -207.5487
#
# $convergence

# quantiles:
q <- c(0.1, 0.3, 0.5, 0.7, 0.9)
## observed data:
(p11_q_c <- quantile(p11[p11$corr == 2, "rt"], probs = q))
#    10%    30%    50%    70%    90%
# 0.4900 0.5557 0.6060 0.6773 0.8231
(p11_q_e <- quantile(p11[p11$corr == 1, "rt"], probs = q))
#    10%    30%    50%    70%    90%
# 0.4908 0.5391 0.5905 0.6413 1.0653
### LBA:
# predicted error rate
(pred_prop_correct_lba <- pLBA(Inf, 2,
                               A = p11_norm$par["A"],
                               b = p11_norm$par["A"]+p11_norm$par["b"],
                               t0 = p11_norm$par["t0"],
                               mean_v = c(p11_norm$par["v1"], p11_norm$par["v2"]),
                               sd_v = c(1, p11_norm$par["sv"])))
# [1] 0.9581342
(pred_correct_lba <- qLBA(q*pred_prop_correct_lba, response = 2,
                          A = p11_norm$par["A"],
                          b = p11_norm$par["A"]+p11_norm$par["b"],
                          t0 = p11_norm$par["t0"],
                          mean_v = c(p11_norm$par["v1"], p11_norm$par["v2"]),
                          sd_v = c(1, p11_norm$par["sv"])))
# [1] 0.4871710 0.5510265 0.6081855 0.6809796 0.8301286
(pred_error_lba <- qLBA(q*(1-pred_prop_correct_lba), response = 1,
                        A = p11_norm$par["A"],
                        b = p11_norm$par["A"]+p11_norm$par["b"],
                        t0 = p11_norm$par["t0"],
                        mean_v = c(p11_norm$par["v1"], p11_norm$par["v2"]),
                        sd_v = c(1, p11_norm$par["sv"])))
# [1] 0.4684374 0.5529575 0.6273737 0.7233961 0.9314820
### diffusion:
# same result as when using Inf, but faster:
(pred_prop_correct_diffusion <- pdiffusion(rt = 20,  response = "upper",
                                      a=p11_diff$par["a"],
                                      v=p11_diff$par["v"],
                                      t0=p11_diff$par["t0"],
                                      sz=p11_diff$par["sz"],
                                      st0=p11_diff$par["st0"],
                                      sv=p11_diff$par["sv"]))
# [1] 0.964723
(pred_correct_diffusion <- qdiffusion(q*pred_prop_correct_diffusion,
                                      response = "upper",
                                      a=p11_diff$par["a"],
                                      v=p11_diff$par["v"],
                                      t0=p11_diff$par["t0"],
                                      sz=p11_diff$par["sz"],
                                      st0=p11_diff$par["st0"],
                                      sv=p11_diff$par["sv"]))
# [1] 0.4748271 0.5489903 0.6081182 0.6821927 0.8444566
(pred_error_diffusion <- qdiffusion(q*(1-pred_prop_correct_diffusion),
                                    response = "lower",
                                    a=p11_diff$par["a"],
                                    v=p11_diff$par["v"],
                                    t0=p11_diff$par["t0"],
                                    sz=p11_diff$par["sz"],
                                    st0=p11_diff$par["st0"],
                                    sv=p11_diff$par["sv"]))
# [1] 0.4776565 0.5598018 0.6305120 0.7336275 0.9770047
### plot predictions
par(mfrow=c(1,2), cex=1.2)
plot(p11_q_c, q*prop.table(table(p11$corr))[2], pch = 2, ylim=c(0, 1), xlim = c(0.4, 1.3), ylab = "Cumulative Probability", xlab = "Response Time (sec)", main = "LBA")
points(p11_q_e, q*prop.table(table(p11$corr))[1], pch = 2)
lines(pred_correct_lba, q*pred_prop_correct_lba, type = "b")
lines(pred_error_lba, q*(1-pred_prop_correct_lba), type = "b")
legend("right", legend = c("data", "predictions"), pch = c(2, 1), lty = c(0, 1))

plot(p11_q_c, q*prop.table(table(p11$corr))[2], pch = 2, ylim=c(0, 1), xlim = c(0.4, 1.3), ylab = "Cumulative Probability", xlab = "Response Time (sec)", main = "Diffusion")
points(p11_q_e, q*prop.table(table(p11$corr))[1], pch = 2)
lines(pred_correct_diffusion, q*pred_prop_correct_diffusion, type = "b")
lines(pred_error_diffusion, q*(1-pred_prop_correct_diffusion), type = "b")


rand_rts <- rdiffusion(1e5, a=p11_diff$par["a"],
                            v=p11_diff$par["v"],
                            t0=p11_diff$par["t0"],
                            sz=p11_diff$par["sz"],
                            st0=p11_diff$par["st0"],
                            sv=p11_diff$par["sv"])
plot(ecdf(rand_rts[rand_rts$response == "upper","rt"]))
normalised_pdiffusion = function(rt,...) pdiffusion(rt,...)/pdiffusion(rt=Inf,...)
curve(normalised_pdiffusion(x, response = "upper",
                            a=p11_diff$par["a"],
                            v=p11_diff$par["v"],
                            t0=p11_diff$par["t0"],
                            sz=p11_diff$par["sz"],
                            st0=p11_diff$par["st0"],
                            sv=p11_diff$par["sv"]),
      add=TRUE, col = "yellow", lty = 2)















#### START HERE!!
# ---------------------------------------------------------------------
# Diffusion model analysis
# http://singmann.org/rtdists-0-7-2-response-time-distributions-now-with-rcpp-and-faster/

# can we recover the parameters?
ll_diffusion <- function(pars, rt, response) {
  densities <- tryCatch(
     ddiffusion(rt, response = response,
                          a=pars[1],
                          v=pars[2],
                          t0=pars[3],
                          sz=pars[4],
                          st0=pars[5],
                          sv=pars[6]),
                          error = function(e) 0)
  if (any(densities == 0)) return(1e6)
  return(-sum(log(densities)))
}

# a: threshold separation.  Amount of information that is considered for a decision.
# Large values indicate a conservative decisional style. Typical range: 0.5 < a < 2.

# v: drift rate.   Average slope of the information accumulation process.
# The drift gives information about the speed and direction of the accumulation of
# information.  Large (absolute) values of drift indicate a good performance.
# If receivedinformation supports the response linked to the upper threshold the sign
# will be positive and vice versa. Typical range: -5 < v < 5.

# t0: non-decision  time  or  response  time  constant  (in  seconds).
# Lower  bound  for the duration of all non-decisional processes (encoding and
# response execution). Typical range: 0.1 < t0 < 0.5.



df <- all_trials %>%
  dplyr::filter(response != "NORESP")

df$response <- factor(ifelse(df$response == "right", "upper", "lower"))
# create numeric response variable where 1 is an error and 2 a correct response
df$corr <- df$correct + 1

# transform rt in seconds
df$rt <- df$rt / 1000
df$rt <- ifelse(df$rt  < 0.2 | df$rt > 2.0, NA, df$rt )

# remove NAs
df_clean <- df[!is.na(df$rt), ]

df_clean %>%
  group_by(id, is_surprise_clip) %>%
  summarise(
    n = n()
    )

# -------------------------------------------------------------------
# surprise trials
n_subj <- length(unique(df_clean$id))
results_surprise <- vector("list",  n_subj)

for (id_index in 1:n_subj) {
  df <- df_clean %>%
    # dplyr::filter(id == id_index & is_surprise_clip == 1)
    dplyr::filter(id == id_index & is_congruent_trial == 1)
  start <- c(runif(2,  0.5,  3),  0.1,  runif(3,  0,  0.5))
  names(start) <- c("a",  "v",  "t0",  "sz",  "st0",  "sv")
  recov <- nlminb(start,  ll_diffusion,  lower = 0,  rt=df$rt,  response=df$response)
  # round(recov$par, 3)
  results_surprise[[id_index]] <- recov$par
  message('Processing participant ', id_index, ' of ',  n_subj)
}


# not surprise trials
n_subj <- length(unique(df_clean$id))
results_not_surprise <- vector("list",  n_subj)

for (id_index in 1:n_subj) {
  df <- df_clean %>%
    # dplyr::filter(id == id_index & is_surprise_clip == 0)
    dplyr::filter(id == id_index & is_congruent_trial == 0)
  start <- c(runif(2,  0.5,  3),  0.1,  runif(3,  0,  0.5))
  names(start) <- c("a",  "v",  "t0",  "sz",  "st0",  "sv")
  recov <- nlminb(start,  ll_diffusion,  lower = 0,  rt=df$rt,  response=df$response)
  # round(recov$par, 3)
  results_not_surprise[[id_index]] <- recov$par
  message('Processing participant ', id_index, ' of ',  n_subj)
}

# transform lists in data.frames
df_surpr <- data.frame(matrix(unlist(results_surprise), nrow = n_subj, byrow = TRUE))
df_not_surpr <- data.frame(matrix(unlist(results_not_surprise), nrow = n_subj,
                                                byrow = TRUE))

foo <- data.frame(df_surpr$X6, df_not_surpr$X6)
colMeans(foo)
x_dif <- foo[, 1] - foo[, 2]
hist(x_dif)
t.test(x_dif)



# check fit
plot(p11_q_c, q*prop.table(table(p11$corr))[2], pch = 2, ylim=c(0, 1), xlim = c(0.4, 1.3), ylab = "Cumulative Probability", xlab = "Response Time (sec)", main = "Diffusion")
points(p11_q_e, q*prop.table(table(p11$corr))[1], pch = 2)
lines(pred_correct_diffusion, q*pred_prop_correct_diffusion, type = "b")
lines(pred_error_diffusion, q*(1-pred_prop_correct_diffusion), type = "b")









require("rtdists")
require("dplyr")   # for data manipulations and looping
require("tidyr")   # for data manipulations
require("purrr")   # for data manipulations
require("lattice") # for plotting and corresponding themes
require("latticeExtra")
lattice.options(default.theme = standard.theme(color = FALSE))
lattice.options(default.args = list(as.table = TRUE))
options(digits = 3) # only three decimal digits
require("binom")  # for binomial confidence intervals


# clean
df1 <- all_trials %>%
   dplyr::filter(subject_name != "05MR06M32" &
                           response != "NORESP" &
                           trials_after_clip < 5 &
                           rt > 200 & rt < 1500
                         )

df1$rt <- df1$rt / 1000
df1$response_num <- as.numeric(df1$response)

df <- df1[!is.na(df1$rt), ]




# aggregate data for first plot:
agg_rr98 <- df  %>%
  group_by(trials_after_clip, is_surprise_clip, is_congruent_trial) %>%
  summarise(
    prop = mean(correct),
    mean_rt = mean(rt),
    median_rt = mean(rt)) %>%
  ungroup()

xyplot(prop ~ trials_after_clip | is_surprise_clip,
  agg_rr98,
  group = is_congruent_trial,
  type = "b",
  auto.key = list(lines = TRUE), ylab = "Proportion of 'correct' responses")


quantiles <- c(0.1, 0.3, 0.5, 0.7, 0.9)
## aggregate data for quantile plot
quantiles_rr98 <- df  %>%
  group_by(is_surprise_clip, trials_after_clip, is_congruent_trial) %>%
  nest() %>%
  mutate(quantiles = map(data, ~ as.data.frame(t(quantile(.x$rt, probs = quantiles))))) %>%
  unnest(quantiles) %>%
  gather("quantile", "rt",`10%`:`90%`) %>%
  arrange(is_surprise_clip, trials_after_clip)

xyplot(rt ~ trials_after_clip| is_surprise_clip,
             quantiles_rr98, group = quantile, type = "b",
       auto.key = list(lines = TRUE), ylab = "RT (in seconds)", subset = is_congruent_trial == 1)

xyplot(rt ~ trials_after_clip| is_surprise_clip,
             quantiles_rr98, group = quantile, type = "b",
       auto.key = list(lines = TRUE), ylab = "RT (in seconds)", subset = is_congruent_trial == 0)






df$instruction <- df$is_surprise_clip
df$response <- df$correct
df$response <- factor(ifelse(df$response == 1, 'upper', 'lower'))

d_nested <- df %>%
  dplyr::filter(id < 3) %>%
  group_by(id, instruction) %>% # we loop across both, id and instruction
  nest()
d_nested


# objective function for diffusion with 1 a. loops over drift to assign drift rates to strength
objective_diffusion_separate <- function(pars, rt, response, drift, ...) {
  non_v_pars <- grep("^v", names(pars), invert = TRUE, value = TRUE)
  base_par <- length(non_v_pars)  # number of non-drift parameters
  densities <- vector("numeric", length(rt))
  for (i in seq_along(levels(drift))) {
    densities[drift == levels(drift)[i]] <-
      ddiffusion(rt[drift == levels(drift)[i]], response=response[drift == levels(drift)[i]],
                 a=pars["a"],
                 t0=pars["t0"],
                 sv=pars["sv"],
                 sz=if ("sz" %in% non_v_pars) pars["sz"] else 0.1,
                 z=if ("z" %in% non_v_pars) pars["z"]*pars["a"] else 0.5*pars["a"],
                 st0=if ("st0" %in% non_v_pars) pars["st0"] else 0,
                 v=pars[base_par+i])
  }
  if (any(densities == 0)) return(1e6)
  return(-sum(log(densities)))
}

# function that creates random start values, also
get_start <- function(base_par, n_drift = 5) {
  start1 <- c(
    a = runif(1, 0.5, 3),
    a_1 = runif(1, 0.5, 3),
    a_2 = runif(1, 0.5, 3),
    t0 = runif(1, 0, 0.5),
    z = runif(1, 0.4, 0.6),
    sz = runif(1, 0, 0.5),
    sv = runif(1, 0, 0.5),
    st0 = runif(1, 0, 0.5),
    d = rnorm(1, 0, 0.05)
  )
  start2 <- sort(rnorm(n_drift), decreasing = FALSE)
  names(start2) <- paste0("v_", seq_len(n_drift))
  c(start1[base_par], start2)
}

# function that tries different random start values until it works:
ensure_fit <-
  function(data, start_function, objective_function,
           base_pars, n_drift = 5, n_fits = 1,
           lower = c(rep(0, length(base_pars)), -Inf,
                     rep(-Inf,length(start_function(base_pars))-length(base_pars)))) {
    best_fit <- list(objective = 1e+06)
  for (i in seq_len(n_fits)) {
    start_ll <- 1e+06
    #browser()
    while(start_ll == 1e+06) {
      start <- start_function(base_pars, n_drift=n_drift)
      start_ll <- objective_function(start,
                                     rt = data$rt, response = data$response_num,
                                     drift = factor(data$trials_after_clip, seq_len(n_drift)),
                                     instruction = data$instruction)
    }
    cat("\nstart fitting.\n") # just for information to see if it is stuck

    fit <- nlminb(start, objective_function,
                  rt = data$rt, response = data$response_num,
                  drift = factor(data$trials_after_clip, seq_len(n_drift)),
                  instruction = data$instruction,
                  lower = lower)

    if (fit$objective < best_fit$objective) best_fit <- fit
  }
  out <- as.data.frame(t(unlist(best_fit[1:3])))
  colnames(out) <- sub("par.", "", colnames(out))
  out
}




fit_diffusion <- d_nested %>%
  mutate(fit =
           map(data,
               ~ensure_fit(data = ., start_function = get_start,
                            objective_function = objective_diffusion_separate,
                            base_pars = c("a", "t0", "sv", "sz", "z")))) %>%
  unnest(fit)
