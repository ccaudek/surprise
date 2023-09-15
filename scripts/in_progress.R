
dat_s <- dt_cor |> 
  dplyr::filter(experiment == "surprise")
dat_s$blk <- factor(dat_s$block)

mod1 <- brm(
  rt ~ is_congruent_trial * blk +
    (1 + is_congruent_trial + blk | subj_id) 
    # (1 | movie_id2),
  family = exgaussian(),
  data = dat_s,
  # algorithm = "meanfield",
  backend = "cmdstanr"
)

pp_check(mod1)
print(mod1, 4)


pp_check(m1)
summary(m1)

conditional_effects(m1, "blk:is_congruent_trial")
conditional_effects(m1, "is_surprise_clip")

conditions <- make_conditions(m1, "is_surprise_clip")
conditional_effects(m1, "blk:is_congruent_trial", conditions=conditions)


