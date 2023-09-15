# -------------------------------------------------------------------

correct_df <- df %>%
  dplyr::filter(correct == 1)

correct_df$blockf <- factor(correct_df$block)
correct_df$trials_after_clip <- factor(correct_df$trials_after_clip)

correct_df$n_after <- as.numeric(correct_df$trials_after_clip)
correct_df$congr <- as.numeric(correct_df$is_congruent_trial)

# maximum random structure
fm <- lmer(vel ~ is_first_trial * condition * is_surprise_clip * block +
  (1 + is_first_trial * condition * is_surprise_clip * block | id) +
  (1 | video_id),
  data = correct_df)
# Model failed to converge: degenerate  Hessian with 5 negative eigenvalues

fm1 <- lmer(vel ~ is_first_trial * condition * is_surprise_clip * block +
  (1 + is_first_trial * condition * block + is_surprise_clip | id) +
  (1 | video_id),
  data = correct_df)
# Model failed to converge: degenerate  Hessian with 1 negative eigenvalues

# control condition
# correct_df <- correct_df[!is.na(correct_df$n_after) & correct_df$block < 5, ]
#
# fm2 <- lmer(vel ~ poly(n_after, 3) * condition * block +
#   (1 + n_after + condition + block || id) +
#   (1 | video_id),
#   data = correct_df)

fm2 <- lmer(vel ~ poly(n_after, 3) * condition * is_surprise_clip * block +
  (1 + n_after + condition + block + is_surprise_clip | id) +
  (1 | video_id),
  data = correct_df)

# show the model
print(summary(fm2), corr = FALSE)
Anova(fm2)

mydf <- ggpredict(fm2, terms = c("n_after", "block", "condition"), ci.lvl = 0.68)
plot(mydf, ci = FALSE)

mydf <- ggpredict(fm2, terms = c("block", "condition"), ci.lvl = 0.68)
plot(mydf, ci = FALSE)

fm3 <- lmer(vel ~ is_clip_trial * condition * block + is_surprise_clip + i_trial +
  (1 + is_clip_trial + condition + block + is_surprise_clip + i_trial | id) +
  (1 | video_id),
  data = correct_df)


# As suggested by the Box-Cox test (R implementation in the car package by Fox & Weisberg, 2011), the response time latencies (RTs) were subjected to an inverse transformation.
#Analysis. A linear mixed effects model (LME) was fit to the data, using the lme4 package (Bates, Maechler, Bolker, & Walker, 2015) in the R software environment (R Core Team, 2017). We defined trial order as covariate and used its by-participant adjustments (random effects) for both the Intercept and the Slope as indicators of individuals’ average speed and change in speed. The dependent variable was reaction time latency, which had been log-transformed to facilitate statistical analysis (as suggested by the Box-Cox test). The results showed a significant effect of trial order (Estimate = -0.007, SE = 0. 003, t = -2.2). Estimates of participants’ coefficients for the Intercept and Slope were also justified (Chi-sq.(1) = 370.2, p < 0.001, AIC-difference = 368). We extracted these random intercepts and random slopes from the model for use as indicators of individual differences in average speed and change in speed.


#  We find longer RTs for *congruent* than incongruent trials. This difference decrease with
# practice (blocks). This result could be explained by saying that the surprising events
# induce a strategy leading subjects to expect some 'incongruency' in the stimuli.  when
# they don't find it, they look for a longer time. This result would imply that a surprising
# event modifies the cognitive strategy used by subject in processing the stimuli.
# It is necessary a control in which the unsurprising videos produce the normal result:
# longer RT for incongruent than for congruent trials.

fm2 <- lmer(vel ~ is_first_trial * condition * is_surprise_clip * block +
  (1 + is_first_trial * condition + block + is_surprise_clip || id) +
  (1 | video_id),
  data = correct_df)

fm2 <- lmer(vel ~ poly(n_after, 3) * condition * block + is_surprise_clip + i_trial +
  (1 + n_after + condition + block + is_surprise_clip + i_trial | id) +
  (1 | video_id),
  data = correct_df)

# model criticism trimming
trimmed_df <- correct_df[abs(scale(resid(fm2))) < 2.5, ]
# percentage of deleted trials
(1 - dim(trimmed_df)[1] / dim(correct_df)[1]) * 100

fm2t <- lmer(vel ~ poly(n_after, 3) * condition * block + is_surprise_clip + i_trial +
  (1 + poly(n_after, 3) + condition + block + is_surprise_clip + i_trial | id) +
  (1 | video_id),
  data = trimmed_df)

mydf <- ggpredict(fm3, terms = c("block", "is_clip_trial", "condition"))
plot(mydf)

Anova(fm3)

# -------------------------------------------------------------------
# GAM

correct_df$NewTimeSeries <- ifelse(correct_df$i_trial == 1, TRUE, FALSE)

mod1 <- bam(vel ~ poly(n_after, 3) * condition * block + is_surprise_clip + i_trial +
   s(i_trial, id, bs='fs', m=1, k=5) +
   s(trials_after_clip, id, bs='fs', m=1, k=2) +
   s(block, id, bs='fs', m=1, k=2) +
   s(condition, id, bs = "fs", m = 1) +
   s(block, id, bs = "fs", m = 1),
   rho = 0.3, AR.start = correct_df$NewTimeSeries,
   data = correct_df)

# model criticism trimming
trimmed_df <- correct_df[abs(scale(resid(mod1))) < 2.5, ]
# percentage of deleted trials
(1 - dim(trimmed_df)[1] / dim(correct_df)[1]) * 100

mod1t <- bam(vel ~ trials_after_clip * condition * block + is_surprise_clip + i_trial +
   s(i_trial, id, bs='fs', m=1, k=5) +
   s(trials_after_clip, id, bs='fs', m=1, k=2) +
   s(block, id, bs='fs', m=1, k=2) +
   s(condition, id, bs = "fs", m = 1) +
   s(block, id, bs = "fs", m = 1),
   rho = 0.3, AR.start = trimmed_df$NewTimeSeries,
   data = trimmed_df)



par(mfrow=c(1,2), cex=1.1)
# including random effects:
plot_smooth(mod1, view="n_after", cond=list(congr=1))
# excluding random effects:
plot_smooth(mod1, view="n_after", cond=list(congr=0))
plot_smooth(mod1, view="i_trial")

plot_data(mod1, view="i_trial", split_by="congr", cex=.5)

plot_smooth(mod1, view="n_after", plot_all="congr", rm.ranef=TRUE, ylim=c(1.4, 2))

fvisgam(mod1, view=c("n_after", "block"), cond=list(congr=0))

anova(mod1, test = "F")




# -------------------------------------------------------------------
# Summary stats on RTs.
out <- df %>%
  dplyr::filter(correct == 1 & trials_after_clip < 4) %>%
  group_by(is_surprise_clip, trials_after_clip) %>%
  summarise(
    m = mean(rt, na.rm = TRUE, trim = 0.1),
    se = sd(rt) / n(),
    n = n()
    )

out$is_surprise_clip <- factor(out$is_surprise_clip)

pd <- position_dodge(0.1) # move them .05 to the left and right
ggplot(out, aes(x = trials_after_clip, y = m, colour = is_surprise_clip)) +
    geom_errorbar(aes(ymin=m-2*se, ymax=m+2*se), width=.0, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd, size = 3) +
    #facet_wrap( ~ is_congruent_trial) +
    labs(x = "Trial position after the clip",
         y = "Reaction Times (ms)",
         title = "Effect of surprise on response latencies") +
   theme_ipsum_rc(base_family = "Helvetica", grid="XY")

# -------------------------------------------------------------------
#  PEA
out <- df %>%
  dplyr::filter(correct == 0 & trials_after_clip != 4) %>%
  group_by(is_surprise_clip, trials_after_clip) %>%
  summarise(
    m = mean(pea, na.rm = TRUE, trim = 0.1),
    se = sd(pea, na.rm = TRUE) / n(),
    n = n()
    )

out$is_surprise_clip <- factor(out$is_surprise_clip)

pd <- position_dodge(0.1) # move them .05 to the left and right
ggplot(out, aes(x = trials_after_clip, y = m, colour = is_surprise_clip)) +
    geom_errorbar(aes(ymin=m-se, ymax=m+se), width=.0, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd, size = 3) +
    #facet_wrap( ~ is_congruent_trial) +
    labs(x = "Trial position after the clip",
         y = "Robust PEA",
         title = "Post-error adjustment") +
   theme_ipsum_rc(base_family = "Helvetica", grid = "XY")

df$trials_after_clip_fac <- factor(df$trials_after_clip)

err_trials <- df %>%
  dplyr::filter(correct == 0 & trials_after_clip != 4)

fm <- lmer(pea ~ is_surprise_clip * trials_after_clip_fac +
  (1 + trials_after_clip_fac || subject_name) + (1 | video_id),
  data = err_trials)
# not enough trials!! Ignore.

# -------------------------------------------------------------------
#  RT after an error
out_err <- df %>%
  dplyr::filter(correct == 0 & trials_after_clip != 4) %>%
  group_by(is_surprise_clip, trials_after_clip) %>%
  summarise(
    m = mean(rt_follow1, na.rm = TRUE, trim = 0.1),
    se = sd(rt_follow1) / n(),
    n = n()
    )
out_err


out_cor <- df %>%
  dplyr::filter(correct == 1 & trials_after_clip != 4) %>%
  group_by(subject_name, is_surprise_clip, trials_after_clip) %>%
  summarise(
    m_c = mean(rt_follow1, na.rm = TRUE, trim = 0.1),
    se = sd(rt_follow1) / n(),
    n = n()
    )


out_err <- df %>%
  dplyr::filter(correct == 0 & trials_after_clip != 4) %>%
  group_by(subject_name, is_surprise_clip, trials_after_clip) %>%
  summarise(
    m = mean(rt_follow1, na.rm = TRUE, trim = 0.1),
    n = n()
    )

# difference with what happens with correct trials, by subject and condition!!!
# then, analysis on the difference, weighted by the number of trials.

new <- semi_join(out_err, out_cor,
   by = c("subject_name", "is_surprise_clip", "trials_after_clip"))

new$dif <- new$m - new$m_c

fm <- lmer(dif ~ is_surprise_clip * trials_after_clip +
   (1 | subject_name),
   data = new)


# Better summary stats.
df$ntrial <- factor(df$trials_after_clip)
temp <- df %>% dplyr::filter(trials_after_clip != 4 & rt > 200)
fm <- lmer(log(rt) ~ ntrial * is_surprise_clip +
   (1 + ntrial | subject_name) + (1 | video_id),
   data = temp)
summary(fm)


### Accuracy by participant.

df_acc <- df %>%
    group_by(subject_name) %>%
    summarise(
        p = mean(correct, na.rm = TRUE)
        )

# data.frame(df_acc)


# Accuracy as a function of trial type and surprise.
df %>%
    filter(rt > 200 & rt < 2000) %>%
    group_by(is_clip_trial, is_surprise_clip) %>%
    summarise(
        p = mean(correct),
        m = mean(rt, na.rm = TRUE, trim = 0.10)
        )

### Cleaned df.
cleaned_df <- df %>%
    dplyr::filter(rt > 200 & rt < 2000 &
                  correct == 1 &
                  trials_after_clip < 4)

hist(log(cleaned_df$rt))

tapply(cleaned_df$rt, list(cleaned_df$block, cleaned_df$is_congruent_trial), mean)


# Effect of surprise on response latencies
# ++++++++++++++++++++++++++++++++
# As a function of block, surprise, and
# trial number after the clip.
out <- df %>%
    dplyr::filter(rts > .200 & rts < 2.0 &
                  correct == 1) %>%
    group_by(is_surprise_clip,trials_after_clip, is_congruent_trial) %>%
    summarise(
        m = mean(log(rts), na.rm = TRUE, trim = 0.05),
        # m = mean(rt, na.rm = TRUE, trim = 0.05),
        se = sqrt(var(log(rts)) / n())
        )
out$is_surprise_clip <- factor(out$is_surprise_clip)

pd <- position_dodge(0.1) # move them .05 to the left and right
ggplot(out, aes(x = trials_after_clip, y = m, colour = is_surprise_clip)) +
    geom_errorbar(aes(ymin=m-se, ymax=m+se), width=.0, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd, size = 3) +
    facet_wrap( ~ is_congruent_trial) +
    labs(x = "Trial position after the clip",
         y = "Log reaction times",
         title = "Effect of surprise on response latencies")


# ++++++++++++++++++++++++++++++++
out <- df %>%
    dplyr::filter(rt > 200 & rt < 2000 &
                  correct == 1 &
                  trials_after_clip < 4) %>%
    group_by(subject_name, block, is_surprise_clip, trials_after_clip) %>%
    summarise(
        y = mean(log(rt), na.rm = TRUE, trim = 0.05)
        )

sur_df <- dplyr::filter(out, is_surprise_clip == 1)
nosur_df <- dplyr::filter(out, is_surprise_clip == 0)

tot_df <- sur_df
tot_df$y <- sur_df$y - nosur_df$y

out1 <- tot_df %>%
    group_by(is_surprise_clip, trials_after_clip) %>%
    summarise(
        m = mean(y, na.rm = TRUE, trim = 0.05),
        se = sqrt(var(y) / n())
        )
out1$is_surprise_clip <- NA

pd <- position_dodge(0.1) # move them .05 to the left and right
ggplot(out1, aes(x = trials_after_clip, y = m)) +
    geom_errorbar(aes(ymin=m-se, ymax=m+se), width=.0, position=pd) +
    #facet_wrap(~ block) +
    geom_line(position=pd) +
    geom_point(position=pd, size = 3) +
    labs(x = "Trial position after the clip",
         y = "Log[RT(Surprise) / RT (Nosurprise)]",
         title = "Surprise Effect") +
   theme_ipsum_rc(base_family = "Helvetica", grid="XY")



# ++++++++++++++++++++++++++++++++

data.frame(
    df %>%
        group_by(subject_name) %>%
        summarise(
            yy = mean(block)
    )
)



out <- df %>%
    dplyr::filter(rt > 200 & rt < 2000 &
                  correct == 1 &
                  trials_after_clip < 4) %>%
    group_by(subject_name, is_congruent_trial, is_surprise_clip, trials_after_clip) %>%
    summarise(
        y = mean(log(rt), na.rm = TRUE, trim = 0.05)
        )

sur_df <- dplyr::filter(out, is_surprise_clip == 1)
nosur_df <- dplyr::filter(out, is_surprise_clip == 0)

tot_df <- sur_df
tot_df$y <- sur_df$y - nosur_df$y

out1 <- tot_df %>%
    group_by(is_congruent_trial, is_surprise_clip, trials_after_clip) %>%
    summarise(
        m = mean(y, na.rm = TRUE, trim = 0.05),
        se = sqrt(var(y) / n())
        )
out1$is_surprise_clip <- NA

pd <- position_dodge(0.1) # move them .05 to the left and right
ggplot(out1, aes(x = trials_after_clip, y = m)) +
    geom_errorbar(aes(ymin=m-se, ymax=m+se), width=.0, position=pd) +
    facet_wrap(~ is_congruent_trial) +
    geom_line(position=pd) +
    geom_point(position=pd, size = 3) +
    labs(x = "Trial position after the clip",
         y = "Log[RT(Surprise) / RT (Nosurprise)]",
         title = "Surprise Effect") +
   theme_ipsum_rc(base_family = "Helvetica", grid="Y")







# Effect of surprise on accuracy
# ++++++++++++++++++++++++++++++++

out <- df %>%
    dplyr::filter(rt > 200 & rt < 2000 & trials_after_clip < 4) %>%
    group_by(is_surprise_clip, trials_after_clip) %>%
    summarise(
        m = mean(correct),
        se = sqrt(var(correct) / n())
        )
out$is_surprise_clip <- factor(out$is_surprise_clip)

# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge(0.1) # move them .05 to the left and right

ggplot(out, aes(x=trials_after_clip, y=m, colour=is_surprise_clip)) +
    #facet_wrap( ~ block) +
    geom_errorbar(aes(ymin=m-se, ymax=m+se), width=.0, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd, size = 5) +
    theme_ipsum_rc(base_family = "Helvetica", grid="XY") +
    labs(x = "Trial position after the clip",
       y = "Accuracy",
       title = "Effect of surprise on accuracy")


# Congruency effect.
out <- cleaned_df %>%
    group_by(is_surprise_clip, is_congruent_trial, block, subject_name) %>%
    summarise(m = mean(rt, na.rm = TRUE, trim = 0.05))

inc_trials <- out[out$is_congruent_trial == 0, ]
con_trials <- out[out$is_congruent_trial == 1, ]

ce_df <- con_trials
ce_df$is_congruent_trial <- NULL
ce_df$cngr_eff <- inc_trials$m - con_trials$m

out <- ce_df %>%
    group_by(block, is_surprise_clip) %>%
    summarise(
        m = mean(cngr_eff),
        se = sqrt(var(cngr_eff) / n())
        )
out$is_surprise_clip <- factor(out$is_surprise_clip)

# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge(0.1) # move them .05 to the left and right

ggplot(out, aes(x=block, y=m, color=is_surprise_clip)) +
    #facet_wrap( ~ fist_trial) +
    geom_errorbar(aes(ymin=m-se, ymax=m+se), width=.0, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd, size = 5) +
    scale_fill_manual(values = c("#ea484b", "#abd9e9")) +
    labs(x = "Blocks of trials",
       y = "Congruency effect",
       title = "Effect of surprise on the congruency effect")



# Effect of surprise on the congruency effect.
out <- cleaned_df %>%
    dplyr::filter(subject_name != "34AM03_F19") %>%
    group_by(subject_name, is_clip_trial, is_surprise_clip, is_congruent_trial,
             block) %>%
    summarise(m = mean(rt, na.rm = TRUE, trim = 0.05))

inc_trials <- out[out$is_congruent_trial == 0, ]
con_trials <- out[out$is_congruent_trial == 1, ]

ce_df <- con_trials
ce_df$is_congruent_trial <- NULL
ce_df$cngr_eff <- inc_trials$m - con_trials$m

out <- ce_df %>%
    group_by(block, is_clip_trial, is_surprise_clip) %>%
    summarise(
        m = mean(cngr_eff),
        se = sqrt(var(cngr_eff) / n())
        )
out$is_clip_trial <- factor(out$is_clip_trial)
out$is_surprise_clip <- factor(out$is_surprise_clip)

# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge(0.3) # move them .05 to the left and right

ggplot(out, aes(x=block, y=m, color=is_surprise_clip)) +
    facet_wrap( ~ is_clip_trial) +
    geom_errorbar(aes(ymin=m-se, ymax=m+se), width=.0, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd, size = 4) +
    scale_fill_manual(values = c("#ea484b", "#abd9e9")) +
    labs(x = "Blocks of trials",
       y = "Congruency effect",
       title = "Effect of surprise on congruency effect") +
    theme_ipsum_rc(grid="Y")


### Effect of surprise on the congruency effect as a function of trial number.
out <- cleaned_df %>%
    dplyr::filter(subject_name != "34AM03_F19") %>%
    group_by(subject_name, is_clip_trial, is_surprise_clip, is_congruent_trial,
             block) %>%
    summarise(m = mean(rt, na.rm = TRUE, trim = 0.05))

inc_trials <- out[out$is_congruent_trial == 0, ]
con_trials <- out[out$is_congruent_trial == 1, ]

ce_df <- con_trials
ce_df$is_congruent_trial <- NULL
ce_df$cngr_eff <- inc_trials$m - con_trials$m

out <- ce_df %>%
    group_by(block, is_clip_trial, is_surprise_clip) %>%
    summarise(
        m = mean(cngr_eff),
        se = sqrt(var(cngr_eff) / n())
        )
out$is_clip_trial <- factor(out$is_clip_trial)
out$is_surprise_clip <- factor(out$is_surprise_clip)

# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge(0.2) # move them .05 to the left and right

ggplot(out, aes(x=block, y=m, color=is_surprise_clip)) +
    facet_wrap( ~ is_clip_trial) +
    geom_errorbar(aes(ymin=m-se, ymax=m+se), width=.0, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd, size = 4) +
    scale_fill_manual(values = c("#ea484b", "#abd9e9")) +
    labs(x = "Blocks of trials",
       y = "Congruency effect",
       title = "Effect of surprise on congruency effect") +
    theme_ipsum_rc(grid="Y")





















# Effect of surprise on the congruency effect
# ++++++++++++++++++++++++++++++++

out <- df %>%
    dplyr::filter(rt > 200 & rt < 2000 & correct == 1 & trials_after_clip < 4) %>%
    group_by(isClipTrial, surpriseClip, congruency, subject_name) %>%
    summarise(m = mean(rt, na.rm = TRUE, trim = 0.10))

inc_trials <- out[out$congruency == 0, ]
con_trials <- out[out$congruency == 1, ]

ce_df <- con_trials
ce_df$concgruency <- NULL
ce_df$cngr_eff <- inc_trials$m - con_trials$m

out <- ce_df %>%
    group_by(surpriseClip, isClipTrial) %>%
    summarise(
        m = mean(cngr_eff),
        se = sqrt(var(cngr_eff) / n())
        )
out$isClipTrial <- factor(out$isClipTrial)
out$surpriseClip <- factor(out$surpriseClip)

# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge(0.25) # move them .05 to the left and right

p <- ggplot(out, aes(x=isClipTrial, y=m, color=surpriseClip)) +
    #facet_wrap( ~ isClipTrial) +
    geom_errorbar(aes(ymin=m-se, ymax=m+se), width=.0, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd, size = 3) +
    #scale_fill_manual(values = c("#ea484b", "#abd9e9")) +
    labs(x = "Blocks of trials",
       y = "Congruency effect",
       title = "Effect of surprise on the congruency effect") +
   scale_fill_ipsum(name="surpriseClip") +
   guides(fill = guide_legend(reverse=TRUE)) +
   theme_ipsum_rc(grid="Y")



out <- df %>%
    dplyr::filter(rt > 200 & rt < 2000 & correct == 1 & trials_after_clip < 4) %>%
    group_by(isClipTrial, surpriseClip, congruency, subject_name) %>%
    summarise(m = mean(rt, na.rm = TRUE, trim = 0.10))

inc_trials <- out[out$congruency == 0, ]
con_trials <- out[out$congruency == 1, ]

ce_df <- con_trials
ce_df$concgruency <- NULL
ce_df$cngr_eff <- inc_trials$m - con_trials$m

out <- ce_df %>%
    group_by(surpriseClip, isClipTrial) %>%
    summarise(
        m = mean(cngr_eff),
        se = sqrt(var(cngr_eff) / n())
        )
out$isClipTrial <- factor(out$isClipTrial)
out$surpriseClip <- factor(out$surpriseClip)

# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge(0.25) # move them .05 to the left and right

ggplot(out, aes(x=isClipTrial, y=m, color=surpriseClip)) +
    geom_errorbar(aes(ymin=m-se, ymax=m+se), width=.0, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd, size = 3) +
    scale_fill_manual(values = c("#ea484b", "#abd9e9")) +
    labs(x = "Does the trial follows the video clip?",
       y = "Congruency effect",
       title = "Effect of surprise on the congruency effect") +
    theme_apa()


# DONE =========================================================================





# Add subjects' names.
mydata$file_name <- file_names

mydata$id <- as.factor(substr(mydata$file_name, 0, 10))

if (cond == "POSTCUE") {
    mydata$cue <- as.factor(substr(mydata$file_name, 14, 20))
} else {
    mydata$cue <- as.factor(substr(mydata$file_name, 14, 19))
}

mydata$gender <- as.factor(substr(mydata$file_name, 8, 8))
mydata$age <- as.numeric(substr(mydata$file_name, 9, 10))

# rename variables
mydata$block <- mydata$V1
mydata$V1 <- NULL

mydata$stim_set <- mydata$V2
mydata$V2 <- NULL

# each participant only sees one set out of 3 sets! In each stimulus set there
# are 3 different identities: frames 1--20: identity 1, frames 21--40: identity
# 2, frames 41--60: identity 3.  In setA, identity 1 is happy (?), identity 2 is
# angry (?), identity 3 is neutral (?)  In setB, identity 1 is ...  In setC,
# identity 1 is ...

mydata$trial_in_block <- mydata$V3
mydata$V3 <- NULL
# trial within each block

mydata$contrast <- mydata$V4
mydata$V4 <- NULL
# contrast of target image

mydata$im1 <- mydata$V5
mydata$im2 <- mydata$V6
mydata$im3 <- mydata$V7
mydata$V5 <- NULL
mydata$V6 <- NULL
mydata$V7 <- NULL
# the three images shown in each trial: im1: one frame randomly chosen from those
# belonging to the identity 1 im2: one frame randomly chosen from those belonging
# to the identity 2 im3: one frame randomly chosen from those belonging to the
# identity 3

mydata$im_trg_char <- mydata$V8
mydata$im_t <- mydata$V9
mydata$V8 <- NULL
mydata$V9 <- NULL
# image target: V8 in the form <img_49>; V9 in the form <49>

mydata$resp <- mydata$V10
mydata$V10 <- NULL
# the image selected by the participant

mydata$im_diff <- mydata$V11
mydata$V11 <- NULL
# the response error, coded as response in the current trial - target in the
# current trial

mydata$on_target_resp <- mydata$V12
mydata$V12 <- NULL
# whether the image selected by the participant is exactly the target image
# (on_target_resp == 1) or not (on_target_resp == 0) Better would be to estimate
# the JND and to code as 1 those responses that are within one JND!!

mydata$trg_id <- mydata$V13
mydata$resp_id <- mydata$V14
mydata$V13 <- NULL
mydata$V14 <- NULL
# target's orientation (V13) and response orientation (V14)

mydata$id_cor <- mydata$V15
mydata$V15 <- NULL
# whether the participant has selected the correct orientation (id_cor == 1)

mydata$rt <- mydata$V16
mydata$V16 <- NULL
# reaction times between the presentation of the stimulus and the participant's
# response (??)

mydata$start_img <- mydata$V17
mydata$V17 <- NULL
# the random image that was presented before the participant starts the
# exploration of the 3 continua in order to select his/her response

mydata$target_orientation <- mydata$V18
mydata$V18 <- NULL

mydata$target_sf <- mydata$V19  # target's spatial frequency
mydata$V19 <- NULL

mydata$response_sf <- mydata$V20  # response's spatial frequency
mydata$V20 <- NULL

mydata$response_orientation <- mydata$V21
mydata$V21 <- NULL

mydata$day <- mydata$V22
mydata$V22 <- NULL

mydata$hour <- mydata$V23
mydata$V23 <- NULL


# Summary statistics ----

mydata$sf_cor <- ifelse(abs(mydata$im_diff) < 3, 1, 0)

gabor <- with(mydata, data.frame(x = id_cor, y = sf_cor))
head(gabor, 5)

t_abs <- with(gabor, table(x, y))
t_abs

n <- with(gabor, sum(table(x, y)))
n

t_rel <- with(gabor, table(x, y) / n)
t_rel
sum(t_rel)

# Distribuzioni marginali.
px <- rowSums(t_rel)
px

ex <- 0 * 0.2095047 + 1 * 0.7904953
ex

py <- colSums(t_rel)
py

ey <- 0 * 0.6184739 + 1 * 0.3815261
ey

t_rel

s_xy <- (0 - ex) * (0 - ey) * 0.20571174 +
        (1 - ex) * (0 - ey) * 0.41276216 +
        (0 - ex) * (1 - ey) * 0.00379295 +
        (1 - ex) * (1 - ey) * 0.37773315
s_xy

with(gabor, cov(x, y))


# SD for the two variables.

vx <- with(gabor, mean(x^2) - mean(x)^2)
std_x <- sqrt(vx)

vy <- with(gabor, mean(y^2) - mean(y)^2)
std_y <- sqrt(vy)

s_xy / (std_x * std_y)

with(gabor, cor(x, y))




# Write CSV file in R.
tem <- with(mydata, data.frame(x=))
write.csv(mydata, file = "gabor.csv", row.names = FALSE)









length(unique(mydata$id))

mydata %>% group_by(block) %>%
    summarise(m = mean(id_cor, na.rm = TRUE))

mydata %>%
  filter(id_cor == 1) %>%
  group_by(block) %>%
  summarise(m = sd(im_diff))



ggplot(data = rDat, aes(abs(im_diff))) +
    geom_histogram() +
    facet_wrap(~contrast) +
    labs(title = "Histogram for Spatial Frequency Errors as a Function of Contrast") +
    labs(x = "Spatial Frequency Difference", y = "Count")

ggplot(data = rDat, aes(x = target_sf, y = response_sf)) +
    geom_point(shape = 1) +
    facet_wrap(~contrast) +
    labs(title = "Histogram for Age") +
    labs(x = "Age", y = "Count")
