library(tidyverse)

sim_data <- read_csv("sim/sim-data.csv", col_types = cols())

data(example_par_tab, package = "serosolver")
par_tab <- example_par_tab %>% filter(names != "phi")

write_csv(par_tab, "sim/fit-par-tab.csv")

prior_version <- 2

sim_agmap <- read_csv("sim/sim-agmap.csv", col_types = cols())

strain_isolation_times <- unique(sim_agmap$inf_times)

# NOTE(sen) Need to load plyr globally because someone doesn't know how to
# import functions in packages properly
library(plyr)
model_func <- serosolver::create_posterior_func(
  par_tab = par_tab,
  titre_dat = sim_data %>% as.data.frame(),
  strain_isolation_times = strain_isolation_times,
  antigenic_map = sim_agmap %>% as.data.frame(),
  version = prior_version
)
detach("package:plyr")

start_prob <- -Inf
while (!is.finite(start_prob)) {
  ## Generating starting antibody kinetics parameters
  start_tab <- serosolver::generate_start_tab(par_tab)

  ## Generate starting infection history
  start_inf <- serosolver::setup_infection_histories_titre(
    sim_data, strain_isolation_times,
    space = 5 * 4, titre_cutoff = 4
  )
  # titre_cutoff-specifies how high the titre must be to imply an infection how
  # many epochs (a particular period of time in history or a person's life) must
  # separate proposed infections. Change this to 20-representing 5 years - 20
  # quarters
  start_prob <- sum(model_func(start_tab$values, start_inf)[[1]])
}

library(plyr)
fit <- serosolver::run_MCMC(
  par_tab = start_tab,
  titre_dat = sim_data %>% as.data.frame(),
  antigenic_map = sim_agmap %>% as.data.frame(),
  strain_isolation_times = strain_isolation_times,
  start_inf_hist = start_inf,
  mcmc_pars = c(
    "iterations" = 5000,
    "popt" = 0.234,
    "opt_freq" = 2000,
    "thin" = 1,
    "adaptive_period" = 5000,
    "save_block" = 1000,
    "thin_hist" = 10,
    "hist_sample_prob" = 1,
    "inf_propn" = 1,
    "move_size" = 3 * 4,
    "hist_opt" = 0,
    "hist_switch_prob" = 0.8,
    "year_swap_propn" = 1
  ),
  filename = "sim/fit",
  CREATE_POSTERIOR_FUNC = serosolver::create_posterior_func,
  version = prior_version
)
detach("package:plyr")
