# Script name: 00_prelims.R
# Project: Surprise
# Script purpose: load packages
# @author: Corrado Caudek <corrado.caudek@unifi.it>
# Date Created: Tue Jul 12 09:35:06 2022
# Last Modified Date: Tue Sep 20 06:38:16 2022
# 
# Notes: 

# Load the necessary libraries
if (!require("pacman")) {
  install.packages("pacman")
}
library("pacman")
pacman::p_load(
  "tidyverse", "knitr", "rmarkdown", "forcats", "gdata",
  "rstan", "brms", "ggeffects", "ggthemes",
  "bayesplot", "ggmcmc", "mgcv", "emmeans",
  "ggthemes", "bayesplot", "grid", "brms",
  "tidybayes", "loo", "tidyr", "cmdstanr"
)

options(show.signif.stars = FALSE)
options(max.print = 999999)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
reuse_models <- TRUE
theme_set(bayesplot::theme_default(base_family = "sans", base_size = 14))

