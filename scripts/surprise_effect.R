library("BEST")


bysub <- four_trials_df %>%
  dplyr::filter(
    trials_after_clip == 0) %>%
  group_by(id, is_surprise_clip) %>%
  summarise(
    y = mean(rtTukey, trim = 0.1, na.rm = TRUE)
  )

bysub_sur <- bysub %>%
  dplyr::filter(is_surprise_clip == "Surprise")

bysub_not <- bysub %>%
  dplyr::filter(is_surprise_clip == "No Surprise")

yd <- bysub_sur$y - bysub_not$y

new <- data.frame(bysub_sur, yd)
new$is_surprise_clip <- NULL
new$y <- NULL

write.csv(new, "surprise_data.csv")


# -------------------------------------

dat <- read.csv("surprise_data.csv")

head(dat)

hist(dat$yd)
boxplot(dat$yd)

t.test(dat$yd)


# Bayesian analysis
priors <- list(muM = 0, muSD = 100)
BESTout <- BESTmcmc(yd, priors = priors)
plot(BESTout)
