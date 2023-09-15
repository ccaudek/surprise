

# simulated Wiener process

library("RWiener")


alpha <- 1.5
tau <- 0.25
beta <- 0.5
delta <- 0.2


# rwiener(n, alpha, tau, beta, delta)
rt_df <- rwiener(150, alpha, tau, beta, delta)
summary(rt_df)
hist(rt_df$q)



fp <- wdm(rt_df)
fp$coefficients



## generate random data
dat <- rbind(cbind(rwiener(100, 2.0, 0.3, .25, 2), group="A"),
             cbind(rwiener(100, 1.5, 0.3, .50, 4), group="B"))

## fit wdm
wdm1 <- wdm(dat)

## extract parameters
coef(wdm1)

## further models
wdm2 <- wdm(dat, beta=.5)
wdm3 <- wdm(dat, alpha=wdm1$coefficients[1], tau=wdm1$coefficients[2],
            beta=wdm1$coefficients[3], xvar="group")


wdm3 <- wdm(dat, beta=wdm1$coefficients[3], xvar="group")

wdm3 <- wdm(dat, xvar="group")
coef(wdm3)







selected_trials <- surprise_df %>% 
  dplyr::filter(
    experiment == "Control" & # "Control" "Surprise"
      trials_after_clip > 0 &trials_after_clip < 4
  )


# have RT_i < RT_c!!
bad <- c(20, 21, 22, 24, 25, 26, 30, 32, 41)

selected_trials <- selected_trials[!selected_trials$id %in% bad, ]


# Remove rows with NA on the rt variable
thedat <- selected_trials[!is.na(selected_trials$rtTukey), ]

# subject's id from 1 to n_subjects
thedat$id <- as.numeric(factor(as.character(thedat$id)))
sort(unique(thedat$id))


# temp <- thedat
# 
# # select one subject
# temp <- thedat %>% 
#   dplyr::filter(id == 5)
# 
# 
# # Put data in proper format wdm()
# temp$resp <- factor(ifelse(temp$correct == 1, "upper", "lower"))
# temp$q    <- temp$rtTukey / 1000
# temp$group <- factor(ifelse(temp$is_congruent_trial == "Congruent", "C", "I"))
# 
# temp <- as.data.frame(temp)
# 
# onesubj <- temp %>% 
#   dplyr::select(q, resp, group) %>%
#   ungroup()
# 
# wdm1 <- wdm(onesubj, beta=wdm1$coefficients[3], xvar="group")
# coef(wdm1)




n_sub <- length(unique(thedat$id))
n_sub

res <- list()

for (i_sub in 1:n_sub) {
  
  # select one subject
  temp <- thedat %>% 
    dplyr::filter(id == i_sub)
  
  # Put data in proper format wdm()
  temp$resp <- factor(ifelse(temp$correct == 1, "upper", "lower"))
  temp$q    <- temp$rtTukey / 1000
  temp$group <- factor(ifelse(temp$is_congruent_trial == "Congruent", "C", "I"))
  
  temp <- as.data.frame(temp)
  
  onesubj <- temp %>% 
    dplyr::select(q, resp, group) %>%
    ungroup()
  
  #wdm1 <- wdm(onesubj, beta=wdm1$coefficients[3], xvar="group")
  wdm1 <- wdm(onesubj, xvar="group")
  
  # save results
  res[[i_sub]] <- coef(wdm1)
  print(i_sub)
  
}


param_names <- c("C:alpha", "C:tau", "C:beta", "C:delta", "I:alpha", "I:tau", 
                 "I:beta", "I:delta")

param_df <- do.call(rbind.data.frame, res)
colnames(param_df) <- param_names


with(param_df, t.test(`C:delta` - `I:delta`))
with(param_df, t.test(`I:alpha` - `C:alpha`))



y <- param_df$`I:alpha` - param_df$`C:alpha`
hist(y)

t.test(y)

bysub <- thedat %>% 
  group_by(id, is_congruent_trial) %>% 
  summarise(
    y = mean(rtTukey, na.rm = TRUE)
  )

c_df <- bysub %>% 
  dplyr::filter(is_congruent_trial == "Congruent")

i_df <- bysub %>% 
  dplyr::filter(is_congruent_trial == "Incongruent")

new <- data.frame(param_df$`I:alpha`, param_df$`C:alpha`, congr_eff=i_df$y - c_df$y)
new <- data.frame(param_df$`I:delta`, param_df$`C:delta`, congr_eff=i_df$y - c_df$y)






# One subject analysis ----------------------------------------------------

thedat <- selected_trials

thedat$resp <- factor(ifelse(thedat$correct == 1, "upper", "lower"))
thedat$q    <- thedat$rtTukey / 1000
thedat$group <- factor(as.character(ifelse(thedat$is_congruent_trial == "Congruent", "C", "I")))

# NO 30, 32 (few trials), 41 (few trials)
# 20, 21, 22, 24, 25, 26, 
# select one subject
temp <- thedat %>%
  dplyr::filter(id == 13) %>%  # 23, 27, 28 only correct responses
  dplyr::select(q, resp, group) %>%
  ungroup()

onesubj <- as.data.frame(temp)
onesubj <- onesubj[!is.na(onesubj$q), ]

onesubj %>% 
  # dplyr::filter(resp == "upper") %>% 
  group_by(group) %>% 
  summarise(
    mrt = mean(q, na.rm = TRUE),
    acc = mean(as.numeric(resp) - 1),
    n = n()
  )






wdm1 <- wdm(
  onesubj, 
  #beta=wdm1$coefficients[3], 
  xvar="group")

wdm1$coefficients



