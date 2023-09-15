#' Surprise Experiment
#' Exploratory analyses for the control condition.
#'
#' File: 20_descr_stats_cntr.R
#'
#' This script was last modified on "Sat May 12 16:25:46 2018"



bysub_acc <- mydata %>% 
  dplyr::filter(!is.na(rtTukey)) %>% 
  group_by(subject_name) %>% 
  summarise(
    avg_acc = mean(correct, na.rm = TRUE),
    n = n()
  )
data.frame(bysub_acc)

hist(bysub_acc$avg_acc)


mydata <- mydata[!is.na(mydata$rtTukey), ]

mydata$rt_bins <- cut(
  mydata$rtTukey,
  c(seq(200, 1200, length.out = 10), 2000),
  include.lowest = TRUE,
  labels = c("q1", "q2", "q3", "q4", "q5", "q6", "q7", 
             "q8", "q9", "q10")
)

# mydata$rt_bins <- cut(
#   mydata$rtTukey,
#   quantile(
#     mydata$rtTukey,
#     seq(0, 1, length.out = 11)
#   ),
#   include.lowest = TRUE,
#   labels = c(
#     "q1", "q2", "q3", "q4", "q5", "q6", "q7",
#     "q8", "q9", "q10"
#   )
# )

ntot <- nrow(mydata)

tapply(mydata$correct, mydata$is_congruent_trial, length)

acc_df <- mydata %>% 
  group_by(is_congruent_trial, rt_bins) %>% 
  summarise(
    acc = mean(correct, na.rm = TRUE),
    n = n() / ntot * 2
  )


ggplot(acc_df, aes(x=rt_bins, y=acc, group=is_congruent_trial)) +
  geom_line(aes(color=is_congruent_trial)) +
  geom_point(aes(color=is_congruent_trial))




# Effect of the experiment on the RTs -------------------------------------

# Are RTs longer, on average, when the video clips signal a 
# surprising event?

hist(final$rtTukey)
rtdata <- complete.cases(final$rtTukey)
retimes::mexgauss(rtdata)

# Parameter of the ex-gaussian distribution in the two conditions
# control condition
#         mu      sigma        tau 
# 0.93562532 0.04419088 0.05892118 

# surprise condition
#         mu      sigma        tau 
# 0.89897785 0.06645159 0.08860212 
# The difference is in the parameters sigma and tau.
# Better do this analysis with a Bayesian hierarchic model!

plot_df <- final %>%
  group_by(block, is_congruent_trial, trials_after_clip) %>%
  summarise(
    y = mean(rtTukey, na.rm = TRUE, trim = 0.05),
    stderr = sqrt(var(rtTukey, na.rm = TRUE) / n())
  )

plot_df$lower <- plot_df$y - plot_df$stderr
plot_df$upper <- plot_df$y + plot_df$stderr

pd <- position_dodge(0.2)

ggplot(plot_df, aes(x = trials_after_clip, y = y, color = is_congruent_trial)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), lwd = 1.05, position = pd) +
  geom_line(position = pd, lwd = 1.05) +
  geom_point(position = pd) +
  ylim(500, 700) +
  facet_wrap(~ block) +
  xlab("Trials after the video clip") +
  ylab("Congruency Effect (ms)") +
  # ggtitle(title)
  geom_hline(yintercept = 0, lty = 2) +
  theme(legend.position = "bottom")







fit1 <- brm(formula = rtTukey ~ is_surprise_trial + (1 + is_surprise_trial | subject_name) ,
            data = data, family = lognormal(),
            warmup = 10000, iter = 20000, chains = 4, cores = 8,
            control = list(adapt_delta = 0.9999))






# Congruency effect for each subject and condition ------------------------

bysubraw <- final %>%
  group_by(
    subject_name, block, trials_after_clip, 
    is_congruent_trial
  ) %>%
  summarise(
    mrt = mean(rtTukey, trim = 0.1, na.rm = TRUE)
  )

bysub_con <- bysubraw %>% dplyr::filter(
  is_congruent_trial == "Congruent"
)
bysub_inc <- bysubraw %>% dplyr::filter(
  is_congruent_trial == "Incongruent"
)

bysub_df <- merge(
  bysub_con,
  bysub_inc,
  by = c("subject_name", "block", 
         "trials_after_clip")
)

bysub_df <- bysub_df %>% 
  mutate(
    interference_eff = mrt.y - mrt.x
  )



plot_df <- bysub_df %>%
    group_by(block, trials_after_clip) %>%
    summarise(
      y = mean(interference_eff, trim = .1, na.rm = TRUE),
      stderr = sqrt(var(interference_eff, na.rm = TRUE) / n()),
      n()
    )
plot_df$lower <- plot_df$y - plot_df$stderr
plot_df$upper <- plot_df$y + plot_df$stderr

pd <- position_dodge(0.2)

ggplot(plot_df, aes(x=trials_after_clip, y=y)) +
  geom_pointrange(aes(ymin=lower, ymax=upper), lwd=1.05, position=pd) +
  geom_line(position=pd, lwd=1.05) +
  geom_point(position=pd, size=5) +
  facet_wrap(~ block) +
  xlab("Trials after the video clip") +
  ylab("Congruency Effect (ms)") +
  # ggtitle(title)
  geom_hline(yintercept=0, lty=2) +
  theme(legend.position="bottom")


plot_df <- bysub_df %>%
  group_by(trials_after_clip) %>%
  summarise(
    y = mean(interference_eff, trim = .1, na.rm = TRUE),
    stderr = sqrt(var(interference_eff, na.rm = TRUE) / n()),
    n()
  )
plot_df$lower <- plot_df$y - plot_df$stderr
plot_df$upper <- plot_df$y + plot_df$stderr

pd <- position_dodge(0.2)

ggplot(plot_df, aes(x=trials_after_clip, y=y)) +
  geom_pointrange(aes(ymin=lower, ymax=upper), lwd=1.05, position=pd) +
  geom_line(position=pd, lwd=1.05) +
  geom_point(position=pd, size=5) +
  xlab("Trials after the video clip") +
  ylab("Congruency Effect (ms)") +
  # ggtitle(title)
  geom_hline(yintercept=0, lty=2) +
  theme(legend.position="bottom")



mydf <- bysub_df %>%
  group_by(subject_name, trials_after_clip) %>%
  summarise(
    y = mean(interference_eff, trim = .1, na.rm = TRUE)
  )


m1 <- brm(formula = y ~ trials_after_clip + 
            (1 + trials_after_clip | subject_name),
          data = mydf, 
          #family = lognormal(),
          warmup = 1000, 
          iter = 2000, 
          chains = 4)

summary(m1)

plot(m1, pars = c("trials_after_clip")) 

plot(marginal_effects(m1, effects = "trials_after_clip"))





bysub_df %>%
  group_by(block, is_clip_trial) %>%
  summarise(
    y = mean(interference_eff, na.rm = TRUE)
  )

length(unique(final$subject_name))


dd <- merge(bysub_df, temp, by = "subject_name")

bysub_df$acc <- temp$acc
bysub_df$med_rt <- temp$med_rt
bysub_df

mean(dd$med_rt) + c(-1, 1) * 2.5 * sd(dd$med_rt)
























# Stimulus sequence
with(
  df,
  data.frame(
    trial_within_blk, is_surprise_clip, is_clip_trial, trials_after_clip
  )[1:50, ]
)


temp <- mydata %>%
  dplyr::filter(
    block < 5
  )

ggerrorplot(temp,
            x = "block", 
            y = "correct",
            desc_stat = "mean_se",
            color = "is_congruent_trial", 
            palette = "jco",
            facet.by = "is_surprise_clip",
            ylab = "Accuracy",
            xlab = "Block",
            position = position_dodge(0.3) # Adjust the space between bars
)



ggerrorplot(temp,
            x = "block", 
            y = "rtTukey",
            desc_stat = "mean_se",
            color = "is_congruent_trial", 
            palette = "jco",
            facet.by = "is_surprise_clip",
            ylab = "Reaction Times (ms)",
            xlab = "Block",
            position = position_dodge(0.3) # Adjust the space between bars
)


ggerrorplot(temp,
            x = "is_clip_trial", 
            y = "rtTukey",
            desc_stat = "mean_se",
            color = "is_congruent_trial", 
            palette = "jco",
            facet.by = "is_surprise_clip",
            ylab = "Reaction Times (ms)",
            xlab = "Block",
            position = position_dodge(0.3) # Adjust the space between bars
)





# Histogram of accuracy
out <- mydata %>%
  group_by(subject_name, is_surprise_clip) %>%
  summarise(
    acc = mean(correct, na.rm = TRUE)
  )
# data.frame(out)

p <- ggdensity(out,
  x = "acc",
  add = "mean",
  rug = FALSE,
  color = "is_surprise_clip",
  fill = "is_surprise_clip",
  palette = c("#00AFBB", "#E7B800", "#FC4E07")
)

p2 <- ggpar(p,
  title = "Density Plot for Accuracy",
  subtitle = "Flanker task after video clip",
  caption = "Surprise project",
  xlab = "Accuracy",
  ylab = "Density",
  legend.title = "Surprise Video"
)
p2


# Histogram of median RTs as a function of surprise
out <- mydata %>%
  dplyr::filter(correct == 1) %>%
  group_by(subject_name, is_surprise_clip) %>%
  summarise(
    median_rt = median(rt)
  )
# data.frame(out)

p <- ggdensity(out,
  x = "median_rt",
  add = "mean", rug = FALSE,
  color = "is_surprise_clip", fill = "is_surprise_clip",
  palette = c("#00AFBB", "#E7B800", "#FC4E07")
)

p2 <- ggpar(p,
  title = "Density Plot for Response Latencies",
  subtitle = "Flanker task after video clip",
  caption = "Surprise project",
  xlab = "Median Reaction Times (ms)",
  ylab = "Density",
  legend.title = "Surprise Video"
)
p2


# Histogram of median RTs as a function of number of trials after the clip

temp <- mydata
temp$f_n_after <- factor(temp$trials_after_clip)
temp <- temp[temp$trials_after_clip != 4, ]

ggplot() + 
  geom_density(data=temp, aes(x=lrt, group=f_n_after, fill=f_n_after),alpha=0.5, adjust=2) + 
  xlab("Log reaction time") +
  ylab("Density")


# -------------------------------------------------------------------
# Log transformation of RTs works best
par(mfrow = c(1, 3))
qqPlot(mydata$trt, ylab = "Response Latencies (s)", main = "Reaction Times")
qqPlot(log(mydata$trt), ylab = "Log Response Latencies (s)", main = "Log Reaction Times")
qqPlot(1 / mydata$trt, ylab = "1 / RT", main = "Reciprocal of Reaction Times")
par(mfrow = c(1, 1))


ggerrorplot(final,
  x = "trials_after_clip", 
  y = "lrt",
  desc_stat = "mean_se",
  color = "is_congruent_trial", 
  palette = "jco",
  # facet.by = "block",
  ylab = "Log RTs",
  xlab = "Trial after the Video Clip",
  position = position_dodge(0.3) # Adjust the space between bars
)



ggerrorplot(four_trials_df,
  x = "trials_after_clip", 
  y = "lrt",
  desc_stat = "mean_se",
  color = "is_congruent_trial", 
  palette = "jco",
  facet.by = "is_surprise_clip",
  ylab = "Log RTs",
  xlab = "Trial after the Video Clip",
  position = position_dodge(0.3) # Adjust the space between bars
)


b1 <- final %>%
  dplyr::filter(block == 1)

ggerrorplot(b1,
  x = "trials_after_clip",
  y = "lrt",
  desc_stat = "mean_se",
  color = "is_congruent_trial", 
  palette = "jco",
  facet.by = "is_surprise_clip",
  ylab = "Log RTs",
  xlab = "Trials after the clip",
  position = position_dodge(0.3) # Adjust the space between bars
)


# -------------------------------------------------------------------
#  table of congruency effect and surprise
out <- final %>%
  dplyr::filter(correct == 1) %>%
  # group_by(block, is_surprise_clip, is_congruent_trial) %>%
  group_by(is_first_trial, is_congruent_trial) %>%
  summarise(
    m = mean(rt, na.rm = TRUE, trim = 0.1),
    se = sqrt(var(rt, na.rm = TRUE) / n()),
    n = n()
  )
data.frame(out)


# -------------------------------------------------------------------
# Plot reaction times
temp <- final %>%
  filter(trials_after_clip == 1)

ggerrorplot(final,
  x = "block",
  y = "lrt",
  desc_stat = "mean_se",
  color = "is_congruent_trial",
  palette = "jco",
  facet.by = "is_surprise_clip",
  ylab = "Log Reaction Times",
  xlab = "Block of Trials",
  position = position_dodge(0.3) # Adjust the space between bars
)


ggerrorplot(final,
            x = "block",
            y = "lrt",
            desc_stat = "mean_se",
            color = "is_surprise_clip",
            palette = "jco",
            facet.by = "is_congruent_trial",
            ylab = "Log Reaction Times",
            xlab = "Block of Trials",
            position = position_dodge(0.3) # Adjust the space between bars
)


bysub_df <- final %>%
  group_by(id, block, is_congruent_trial, is_surprise_clip) %>%
  summarise(
    mrt = mean(lrt, na.rm = TRUE), 
    n = n()
  )

bysub_df$mrt <- ifelse(bysub_df$n < 5, NA, bysub_df$mrt)

con_df <- dplyr::filter(bysub_df, is_congruent_trial == "Congruent")
con_df$is_congruent_trial <- NULL
inc_df <- dplyr::filter(bysub_df, is_congruent_trial == "Incongruent")
inc_df$is_congruent_trial <- NULL

con_df$n <- NULL
inc_df$n <- NULL

int_eff_df <- merge(con_df, inc_df, 
                    by = c("id", "block", "is_surprise_clip"))

int_eff_df$interf_eff <- int_eff_df$mrt.y - int_eff_df$mrt.x


ggerrorplot(int_eff_df,
            x = "block", 
            y = "interf_eff",
            desc_stat = "mean_se",
            color = "is_surprise_clip", 
            palette = "jco",
            ylab = "Log Reaction Times",
            xlab = "Block",
            position = position_dodge(0.3) # Adjust the space between bars
)






b1 <- final %>%
  filter(trials_after_clip == 1)

ggerrorplot(b1,
  x = "block", 
  y = "lrt",
  desc_stat = "mean_se",
  color = "is_congruent_trial", 
  palette = "jco",
  facet.by = c("is_surprise_clip"),
  ylab = "Log Reaction Times",
  xlab = "Block",
  position = position_dodge(0.3) # Adjust the space between bars
)


ggerrorplot(final,
  x = "is_surprise_clip", 
  y = "rt",
  desc_stat = "mean_se",
  palette = "jco",
  # facet.by = "is_surprise_clip",
  ylab = "Reaction Time (ms)",
  xlab = "Condition",
  position = position_dodge(0.3) # Adjust the space between bars
)


ggerrorplot(mydata,
  x = "is_surprise_clip", 
  y = "correct",
  desc_stat = "mean_se",
  palette = "jco",
  # facet.by = "is_surprise_clip",
  ylab = "Accuracy",
  xlab = "Surprise",
  position = position_dodge(0.3) # Adjust the space between bars
)


# -------------------------------------------------------------------
# Plot accuracy
all_trials_df <- df

all_trials_df <- all_trials_df %>%
  mutate(is_surprise_clip = fct_recode(is_surprise_clip,
    "No Surprise" = "No",
    "Surprise" = "Yes"
  ))

all_trials_df <- rename(all_trials_df, Condition = is_congruent_trial)
all_trials_df <- all_trials_df %>%
  mutate(Condition = fct_recode(Condition,
    "Incongruent" = "No",
    "Congruent" = "Yes"
  ))

ggerrorplot(all_trials_df,
  x = "block", y = "correct",
  desc_stat = "mean_se",
  color = "Condition", palette = "jco",
  facet.by = "is_surprise_clip",
  ylab = "Accuracy",
  xlab = "Block of Trials",
  ylim = c(0.75, 1.0),
  position = position_dodge(0.3) # Adjust the space between bars
)

ggerrorplot(all_trials_df,
  x = "is_surprise_clip", y = "correct",
  desc_stat = "mean_se",
  color = "condition", palette = "jco",
  facet.by = "block",
  ylab = "Accuracy",
  xlab = "Video Clip",
  ylim = c(0.75, 1.0),
  position = position_dodge(0.3) # Adjust the space between bars
)

all_trials_df <- all_trials_df[!is.na(all_trials_df$is_clip_prev_1) &
  all_trials_df$block < 5, ]

ggerrorplot(all_trials_df,
  x = "is_clip_prev_1", y = "correct",
  desc_stat = "mean_se",
  color = "condition", palette = "jco",
  facet.by = c("block"),
  ylab = "Accuracy",
  xlab = "First Trial After the Video Clip",
  ylim = c(0.75, 1.0),
  position = position_dodge(0.3) # Adjust the space between bars
)

# -------------------------------------------------------------------

names_df <- c(
  "id", "vel", "block", "i_trial", "video_id", "is_surprise_clip",
  "is_congruent_trial", "trials_after_clip"
)

dd <- dplyr::select(correct_df, names_df)

block_1_df <- dd %>%
  dplyr::filter(block == 4)

subject_means <- block_1_df %>%
  group_by(id, is_congruent_trial) %>%
  summarize(vel = mean(vel, na.rm = TRUE))
subject_means

# barplot <- ggplot(subject_means, aes(x = is_congruent_trial, y = vel)) +
#     stat_summary(
#     geom = "bar",
#     fun.y = "mean",
#     col = "black",
#     fill = "gray70"
#     )
# barplot

subject_means_wide <- spread(subject_means,
  key = is_congruent_trial,
  value = vel,
  sep = "_"
)
subject_means_wide

lims <- c(min(correct_df$vel, na.rm = TRUE), max(correct_df$vel, na.rm = TRUE))
wsplot <-
  ggplot(subject_means_wide, aes(
    x = is_congruent_trial_No,
    y = is_congruent_trial_Yes
  )) +
  geom_point() +
  geom_abline() +
  scale_x_continuous("Incongruent", limits = lims) +
  scale_y_continuous("Congruent", limits = lims) +
  theme(aspect.ratio = 1)
wsplot

# -------------------------------------------------------------------

message("\nDescriptive step: done!")
