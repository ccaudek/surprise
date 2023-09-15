library("tidyverse")
library("rstan")
library("ezStan")

rstan_options(auto_write = TRUE)


# Prep the data for Stan ----

a <- ds

a$lrt <- log(a$rt)

#sort by subject
a %>% 
  arrange(id) -> a

#make sure everything's a factor with sensible level orders
a$session = factor(a$block)
a$warning = factor(a$is_surprise_clip)
a$cuing = factor(a$is_clip_trial)
a$flankers = factor(a$is_congruent_trial)

# a$sstroop = factor(a$sstroop,levels=c('incongruent','congruent'))

#define a function for nicer contrasts
halfhelmert_contrasts = function(...) contr.helmert(...)*.5

# specify the contrast matrix
W = ezStan::get_contrast_matrix(
  data = a, 
  contrast_kind = halfhelmert_contrasts, 
  formula = ~ session*warning*cuing*flankers
)
ncol(W)
#View(head(W))

# compute the unique entries in the model matrix
temp = as.data.frame(W)
temp = tidyr::unite_(data = temp, col = 'combined', from = names(temp))
temp_unique = unique(temp)
W_unique = W[row.names(W) %in% row.names(temp_unique), ]

#double check that contrasts are orthogonal
table(cor(W_unique)) #should be 0s & 1s only

# for each row in W, get its index W_unique
indexWer = match(temp$combined,temp_unique$combined)
indexWrt = match(temp$combined[!a$error], temp_unique$combined)

#get subject integers & indices of subjects into flattened W_unique
Ser = as.numeric(factor(a$id))
indexWerSubj = indexWer + (Ser-1)*nrow(W_unique)

#same for RT, but using accurate trials only
Srt = as.numeric(factor(a$id[!a$error]))
indexWrtSubj = indexWrt + (Srt-1)*nrow(W_unique)

#collapse to just between subjects variables
a %>%
  dplyr::group_by(
    id, experiment
  ) %>%
  dplyr::summarize() ->
  b

#generate between-subjects matrix
B = ezStan::get_contrast_matrix(
  data = b, 
  contrast_kind = halfhelmert_contrasts, 
  formula = ~ experiment
)
nrow(B)
ncol(B)
#View(head(B))

#package in list for Stan
data_for_stan = list(
  #nTrials: num trials
  nER = nrow(a),
  #ER: Error outcomes
  ER = as.numeric(a$error),
  #nRT: number of RT trials
  nRT = nrow(a[!a$error, ]),
  #RT: RT outcomes
  RT = scale(a$lrt[!a$error])[, 1], #note scaling
  #nS: num subjects
  nS = length(unique(a$id)),
  #S: trial-by-trial subject labels
  Ser = as.numeric(factor(a$id)),
  #Sacc: trial-by-trial subject labels
  Srt = as.numeric(factor(a$id[!a$error])),
  #nWerr: num within predictors
  nW = ncol(W_unique),
  rowsW = nrow(W_unique),
  #W: unique within predictors
  W = W_unique,
  indexWerSubj = indexWerSubj,
  indexWrtSubj = indexWrtSubj,
  nB = ncol(B),
  B = B
)



# Compile & sample the model ----
mod = rstan::stan_model(
  file = here("stan", "er_rt_faster_robust.stan") 
)

ezStan::kill_stan()
ezStan::clean_stan()

warmup = 1e3
ezStan::start_stan(
  mod = mod, 
  data = data_for_stan, 
  iter = warmup*2, 
  warmup = warmup, 
  cores = 4, 
  seed_start = 1, 
  args = "init=0, include=F, pars=c('cor_helper','multi_normal_helper')"
)
ezStan::watch_stan()
beepr::beep()

post = ezStan::collect_stan()
#ezStan::clean_stan()


filename = 'rdas/post.rda'
save(post,W,B,file=filename,compress=T)
beepr::beep()

